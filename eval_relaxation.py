"""How far are given 3D structures from an energy minimum? Relax them with GFN2-xTB and r2SCAN-3c.

For every molecule in a per-molecule table from eval_singlepoint.py, the QM9 structure (and, with
--structures qm9 rdkit, the RDKit ETKDG + MMFF one) is optimized from the given geometry, and we report:
    relax_kcal_mol   E(given) - E(optimized), the energy released by relaxing to the nearest minimum
    rmsd_heavy_A / rmsd_all_A   RMSD between given and optimized geometry (same atom order, Kabsch-aligned)
    grad_*           gradient at the given geometry (Eh/bohr): xTB norm, ORCA RMS and max (internal coords)
QM9 structures are B3LYP/6-31G(2df,p) minima, so at another level they relax only by that method's
own differences; this measures how close they are to a minimum of xTB / r2SCAN-3c.

Usage:
    python eval_relaxation.py tables/singlepoint_30mol.csv --tag 30mol [--structures qm9 rdkit] [--parallel 6] [--nprocs 16]

Writes tables/relaxation_<tag>.csv and tables/relaxation_<tag>_summary.csv.
Structures are read from temp/singlepoint_*/; runs go to temp/{xtb,orca}_runs/ as <mol_id>_<kind>_opt.
"""
import argparse
import csv
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import mean, median, stdev as _stdev
import sys

import numpy as np

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "energy"))
from molio import read_structures  # noqa: E402
from orca_energy import HARTREE_TO_KCAL, RUNS as ORCA_RUNS, run_orca  # noqa: E402
from xtb_energy import RUNS as XTB_RUNS, run_xtb  # noqa: E402

def find_structure(mol_id, kind):
    hits = sorted(ROOT.glob(f"temp/singlepoint_*/{mol_id}_{kind}.sdf"))
    assert hits, f"no structure file for {mol_id}_{kind}"
    return next(read_structures(hits[0]))[1]


def stdev(x):
    x = list(x)
    return _stdev(x) if len(x) > 1 else 0.0


def rmsd(atoms_a, atoms_b, heavy_only):
    """Kabsch RMSD between two geometries of the same atoms in the same order."""
    keep = [i for i, (s, _) in enumerate(atoms_a) if not heavy_only or s != "H"]
    a = np.array([atoms_a[i][1] for i in keep])
    b = np.array([atoms_b[i][1] for i in keep])
    a, b = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((a @ rot - b) ** 2).sum(1).mean()))


def xtb_relax(mol_id, kind, atoms):
    name = f"{mol_id}_{kind}_opt"
    e_opt, secs = run_xtb(name, atoms, True, 1, 0, 0, None)
    log = (XTB_RUNS / name / "job.out").read_text()
    # The first total energy / gradient norm printed are for the starting geometry
    e0 = float(re.search(r":: total energy\s+(-?\d+\.\d+)", log).group(1))
    g0 = float(re.search(r":: gradient norm\s+(\d+\.\d+)", log).group(1))
    opt = next(read_structures(XTB_RUNS / name / f"{name}_opt.xyz"))[1]
    return {"e0": e0, "e_opt": e_opt, "grad_norm": g0, "grad_rms": "", "grad_max": "", "opt": opt, "secs": secs}


def orca_relax(mol_id, kind, atoms, nprocs, retries=2):
    name = f"{mol_id}_{kind}_opt"
    wd = ORCA_RUNS / "r2scan3c" / name
    done = (wd / f"{name}_opt.xyz").exists() and "ORCA TERMINATED NORMALLY" in (wd / "job.out").read_text()
    if done:  # reuse a finished run from an earlier, interrupted invocation
        d, h, mi, s, ms = map(int, re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec",
                                             (wd / "job.out").read_text()).groups())
        secs = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000
    for attempt in range(retries + 1):
        if done:
            break
        try:
            _, secs = run_orca(name, atoms, "r2scan3c", True, nprocs, 0, 1, 2000)
            break
        except RuntimeError:
            # mpirun occasionally segfaults on this shared machine; a rerun usually succeeds
            if attempt == retries:
                raise
            print(f"  r2scan3c {name}: ORCA failed, retrying", flush=True)
    log = (wd / "job.out").read_text()
    e0 = float(re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", log).group(1))
    g_rms = float(re.search(r"^\s+RMS gradient\s+(\d+\.\d+)\s+\d+\.\d+\s+(?:YES|NO)", log, re.M).group(1))
    g_max = float(re.search(r"^\s+MAX gradient\s+(\d+\.\d+)\s+\d+\.\d+\s+(?:YES|NO)", log, re.M).group(1))
    opt = next(read_structures(wd / f"{name}_opt.xyz"))[1]
    e_opt = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", log)[-1])
    return {"e0": e0, "e_opt": e_opt, "grad_norm": "", "grad_rms": g_rms, "grad_max": g_max, "opt": opt, "secs": secs}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="per-molecule CSV from eval_singlepoint.py (its molecules are used)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--parallel", type=int, default=6, help="ORCA jobs run at once")
    ap.add_argument("--nprocs", type=int, default=16, help="cores per ORCA job")
    ap.add_argument("--structures", nargs="+", choices=["qm9", "rdkit"], default=["qm9"])
    args = ap.parse_args()
    (ROOT / "tables").mkdir(exist_ok=True)

    mols = {}
    for x in csv.DictReader(open(args.table)):
        mols.setdefault(x["mol_id"], x["smiles"])
    structs = {(m, k): find_structure(m, k) for m in mols for k in args.structures}

    results = {}
    for (m, k), atoms in structs.items():
        results[(m, "GFN2-xTB", k)] = xtb_relax(m, k, atoms)
    print("xTB done", flush=True)

    def run(key):
        r = orca_relax(*key, structs[key], args.nprocs)
        print(f"  r2scan3c {key[0]}_{key[1]}: {r['secs']:.0f} s", flush=True)
        return key, r

    with ThreadPoolExecutor(args.parallel) as pool:
        for (m, k), r in pool.map(run, structs):
            results[(m, "r2SCAN-3c", k)] = r

    rows = []
    for (m, method, k), r in results.items():
        given = structs[(m, k)]
        rows.append({
            "mol_id": m, "smiles": mols[m], "structure": k, "method": method,
            "relax_kcal_mol": round((r["e0"] - r["e_opt"]) * HARTREE_TO_KCAL, 3),
            "rmsd_heavy_A": round(rmsd(given, r["opt"], True), 4),
            "rmsd_all_A": round(rmsd(given, r["opt"], False), 4),
            "grad_norm": r["grad_norm"], "grad_rms": r["grad_rms"], "grad_max": r["grad_max"],
            "opt_time_s": round(r["secs"], 2),
        })
    out = ROOT / "tables" / f"relaxation_{args.tag}.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    summary = []
    for method in ("GFN2-xTB", "r2SCAN-3c"):
        for k in args.structures:
            r = [x for x in rows if x["method"] == method and x["structure"] == k]
            relax = [x["relax_kcal_mol"] for x in r]
            summary.append({
                "method": method, "structure": k, "n": len(r),
                "relax_mean": round(mean(relax), 2), "relax_std": round(stdev(relax), 2),
                "relax_median": round(median(relax), 2), "relax_max": round(max(relax), 2),
                "rmsd_heavy_mean": round(mean(x["rmsd_heavy_A"] for x in r), 3),
                "rmsd_all_mean": round(mean(x["rmsd_all_A"] for x in r), 3),
                "rmsd_heavy_max": round(max(x["rmsd_heavy_A"] for x in r), 3),
            })
    out_s = ROOT / "tables" / f"relaxation_{args.tag}_summary.csv"
    with open(out_s, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary[0].keys())
        w.writeheader()
        w.writerows(summary)
    print(f"wrote {out} and {out_s}")
    for s in summary:
        print(s)


if __name__ == "__main__":
    main()
