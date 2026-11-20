"""Single-point energy comparison: QM9 geometry vs RDKit ETKDG + MMFF geometry, scored by several methods.

For each sampled QM9 molecule, both structures are scored as given (no re-optimization), and
delta = E(RDKit) - E(QM9) in kcal/mol is reported per method (positive: QM9 geometry is lower).
DLPNO-CCSD(T) is the reference; every other method's error is measured against it.

Usage:
    python eval_singlepoint.py --n 10 --seed 0 [--parallel 4] [--nprocs 16] [--tag 10mol] [--exclude tables/singlepoint_10mol.csv]
    python eval_singlepoint.py --combine tables/singlepoint_10mol.csv tables/singlepoint_20mol.csv --tag 30mol

Writes tables/singlepoint_<tag>.csv (per molecule x method) and tables/singlepoint_<tag>_summary.csv.
Structures and run directories go to temp/singlepoint_<tag>/ and temp/{xtb,orca}_runs/.
"""
import argparse
import csv
import random
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import mean, stdev as _stdev

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "energy"))
import compare_conformer as cc  # noqa: E402
from molio import read_structures  # noqa: E402
from orca_energy import HARTREE_TO_KCAL, run_orca  # noqa: E402
from xtb_energy import run_xtb  # noqa: E402

ORCA_METHODS = ["r2scan3c", "qm9", "b3lyp-d3", "dlpno"]
METHOD_LABELS = {
    "mmff": "MMFF94", "xtb": "GFN2-xTB",
    "r2scan3c": "r2SCAN-3c", "qm9": "B3LYP/6-31G(2df,p)", "b3lyp-d3": "B3LYP-D3(BJ)/def2-TZVP",
    "dlpno": "DLPNO-CCSD(T)/cc-pVTZ",
}


def stdev(x):
    x = list(x)
    return _stdev(x) if len(x) > 1 else 0.0


def sample_molecules(n, seed, exclude=()):
    """Random characterized QM9 molecules without formal charges (those come from misperceived SDF bonds)."""
    bad = cc.uncharacterized_ids()
    suppl = Chem.SDMolSupplier(str(cc.DATA / "gdb9.sdf"), removeHs=False)
    rng = random.Random(seed)
    picked = []
    while len(picked) < n:
        mol = suppl[rng.randrange(len(suppl))]
        if mol is None:
            continue
        mol_id = mol.GetProp("_Name")
        if mol_id in bad or mol_id in exclude or mol_id in [m.GetProp("_Name") for m in picked]:
            continue
        if any(a.GetFormalCharge() for a in mol.GetAtoms()):
            continue
        picked.append(mol)
    return picked


def mmff_energy(mol):
    props = AllChem.MMFFGetMoleculeProperties(mol)
    return AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--parallel", type=int, default=4, help="ORCA jobs run at once")
    ap.add_argument("--nprocs", type=int, default=16, help="cores per ORCA job")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="per-molecule CSVs / failure lists from earlier runs; their molecules are not drawn again")
    ap.add_argument("--combine", nargs="+", default=None,
                    help="only merge these per-molecule CSVs into tables/singlepoint_<tag>*.csv; no calculations")
    args = ap.parse_args()
    (ROOT / "tables").mkdir(exist_ok=True)
    if args.combine:
        return combine(args.combine, args.tag or "combined")
    exclude = {"gdb_52425"}
    for p in args.exclude:
        exclude |= set(re.findall(r"gdb_\d+", Path(p).read_text()))
    tag = args.tag or f"{args.n}mol_seed{args.seed}"
    work = ROOT / "temp" / f"singlepoint_{tag}"
    work.mkdir(parents=True, exist_ok=True)

    # 1) Sample molecules and generate RDKit conformers; molecules RDKit cannot embed are recorded
    #    as failures and replaced by the next random draw, so n molecules are always scored
    mols, failures = [], []
    candidates = iter(sample_molecules(args.n * 3, args.seed, exclude=exclude))
    while len(mols) < args.n:
        ref = next(candidates)
        mol_id = ref.GetProp("_Name")
        smiles = Chem.MolToSmiles(Chem.RemoveHs(ref))
        try:
            gen = cc.etkdg_mmff(smiles)
        except AssertionError:
            failures.append(f"{mol_id} {smiles}")
            print(f"{mol_id}  {smiles}  RDKit embedding failed, skipped", flush=True)
            continue
        ref_f, gen_f = work / f"{mol_id}_qm9.sdf", work / f"{mol_id}_rdkit.sdf"
        Chem.MolToMolFile(ref, str(ref_f))
        Chem.MolToMolFile(gen, str(gen_f))
        mols.append({"mol_id": mol_id, "smiles": smiles, "n_heavy": ref.GetNumHeavyAtoms(),
                     "ref": ref, "gen": gen, "files": {"qm9": ref_f, "rdkit": gen_f}})
        print(f"{mol_id}  {smiles}", flush=True)

    rows = []  # (mol_id, method, structure) -> energy (kcal/mol) and time (s)
    energy, timing = {}, {}

    # 2) Fast methods
    for m in mols:
        for kind, mol in (("qm9", m["ref"]), ("rdkit", m["gen"])):
            t0 = time.perf_counter()
            energy[(m["mol_id"], "mmff", kind)] = mmff_energy(Chem.Mol(mol))
            timing[(m["mol_id"], "mmff", kind)] = time.perf_counter() - t0
            atoms = next(read_structures(m["files"][kind]))[1]
            e, secs = run_xtb(f"{m['mol_id']}_{kind}", atoms, False, 1, 0, 0, None)
            energy[(m["mol_id"], "xtb", kind)], timing[(m["mol_id"], "xtb", kind)] = e * HARTREE_TO_KCAL, secs
    print("fast methods done", flush=True)

    # 3) ORCA methods, several jobs at once; cheap methods are queued first
    jobs = [(m, meth, kind) for meth in ORCA_METHODS for m in mols for kind in ("qm9", "rdkit")]

    def run(job):
        m, meth, kind = job
        atoms = next(read_structures(m["files"][kind]))[1]
        e, secs = run_orca(f"{m['mol_id']}_{kind}", atoms, meth, False, args.nprocs, 0, 1, 2000)
        print(f"  {meth:9s} {m['mol_id']}_{kind}: {secs:.0f} s", flush=True)
        return job, e, secs

    with ThreadPoolExecutor(args.parallel) as pool:
        for (m, meth, kind), e, secs in pool.map(run, jobs):
            energy[(m["mol_id"], meth, kind)] = e * HARTREE_TO_KCAL
            timing[(m["mol_id"], meth, kind)] = secs

    # 4) Per-molecule table
    methods = ["mmff", "xtb"] + ORCA_METHODS
    for m in mols:
        rmsd = cc.rdMolAlign.GetBestRMS(Chem.RemoveHs(m["gen"]), Chem.RemoveHs(m["ref"]))
        ref_delta = energy[(m["mol_id"], "dlpno", "rdkit")] - energy[(m["mol_id"], "dlpno", "qm9")]
        for meth in methods:
            d = energy[(m["mol_id"], meth, "rdkit")] - energy[(m["mol_id"], meth, "qm9")]
            ts = [timing.get((m["mol_id"], meth, k)) for k in ("qm9", "rdkit")]
            rows.append({
                "mol_id": m["mol_id"], "smiles": m["smiles"], "n_heavy": m["n_heavy"],
                "heavy_rmsd_A": round(rmsd, 3),
                "method": METHOD_LABELS[meth], "delta_kcal_mol": round(d, 3),
                "error_vs_dlpno": round(d - ref_delta, 3),
                "sp_time_s": round(mean(ts), 4) if None not in ts else "",
            })
    out = ROOT / "tables" / f"singlepoint_{tag}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    if failures:
        (ROOT / "tables" / f"singlepoint_{tag}_rdkit_failures.txt").write_text("\n".join(failures) + "\n")
    print(f"\nRDKit embedding failures: {len(failures)}/{len(failures) + len(mols)} {failures}")
    out_s = write_summary(rows, tag)
    print(f"wrote {out} and {out_s}")


def read_rows(path):
    rows = list(csv.DictReader(open(path)))
    for x in rows:
        x["delta_kcal_mol"], x["error_vs_dlpno"] = float(x["delta_kcal_mol"]), float(x["error_vs_dlpno"])
        x["sp_time_s"] = float(x["sp_time_s"]) if x["sp_time_s"] else ""
    return rows


def write_summary(rows, tag):
    """Mean/std over molecules, per method -> tables/singlepoint_<tag>_summary.csv."""
    ref_sign = {x["mol_id"]: x["delta_kcal_mol"] > 0 for x in rows if x["method"] == METHOD_LABELS["dlpno"]}
    summary = []
    for label in dict.fromkeys(x["method"] for x in rows):
        r = [x for x in rows if x["method"] == label]
        d = [x["delta_kcal_mol"] for x in r]
        abs_err = [abs(x["error_vs_dlpno"]) for x in r]
        t = [x["sp_time_s"] for x in r if x["sp_time_s"] != ""]
        summary.append({
            "method": label, "n": len(r),
            "delta_mean": round(mean(d), 2), "delta_std": round(stdev(d), 2),
            "abs_err_vs_dlpno_mean": round(mean(abs_err), 2), "abs_err_vs_dlpno_std": round(stdev(abs_err), 2),
            "same_sign_as_dlpno": f"{sum((x['delta_kcal_mol'] > 0) == ref_sign[x['mol_id']] for x in r)}/{len(r)}",
            "sp_time_mean_s": round(mean(t), 4) if t else "", "sp_time_std_s": round(stdev(t), 4) if len(t) > 1 else "",
        })
    out_s = ROOT / "tables" / f"singlepoint_{tag}_summary.csv"
    with open(out_s, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary[0].keys())
        w.writeheader()
        w.writerows(summary)
    for s in summary:
        print(s)
    return out_s


def combine(paths, tag):
    """Merge per-molecule tables from several runs (must not share molecules) and summarize them together."""
    rows = [x for p in paths for x in read_rows(p)]
    per_run = [set(x["mol_id"] for x in read_rows(p)) for p in paths]
    assert sum(map(len, per_run)) == len(set().union(*per_run)), "runs share molecules"
    out = ROOT / "tables" / f"singlepoint_{tag}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    print(f"combined {len(set().union(*per_run))} molecules from {len(paths)} runs -> {out}")
    write_summary(rows, tag)


if __name__ == "__main__":
    main()
