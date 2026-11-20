"""GFN2-xTB energies: fast semi-empirical QM, suited to ranking conformers at dataset scale.

Usage:
    python energy/xtb_energy.py examples/gdb_52425_qm9.sdf examples/gdb_52425_etkdg_mmff.sdf [--opt] [--nprocs 4]

Single point by default; --opt relaxes each structure first (seconds per QM9-sized molecule).
Runs go to temp/xtb_runs/<structure>/, and an optimized geometry is written there as
<structure>_opt.xyz. Gas phase; pass --solvent water etc. for ALPB implicit solvation.
"""
import argparse
import os
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from molio import read_structures, write_xyz  # noqa: E402

XTB = Path(os.environ.get("XTB_BIN", "/shared/nas2/thaotn2/anaconda3/envs/xtb/bin/xtb"))
RUNS = Path(__file__).parent.parent / "temp" / "xtb_runs"
HARTREE_TO_KCAL = 627.5095


def run_xtb(name, atoms, opt, nprocs, charge, uhf, solvent):
    wd = RUNS / name
    wd.mkdir(parents=True, exist_ok=True)
    write_xyz(atoms, wd / "input.xyz", comment=name)
    cmd = [str(XTB), "input.xyz", "--gfn", "2", "--chrg", str(charge), "--uhf", str(uhf), "-P", str(nprocs)]
    if opt:
        cmd += ["--opt", "tight"]
    if solvent:
        cmd += ["--alpb", solvent]
    env = dict(os.environ, OMP_NUM_THREADS=str(nprocs), OMP_STACKSIZE="1G")
    t0 = time.time()
    with open(wd / "job.out", "w") as out:
        proc = subprocess.run(cmd, cwd=wd, env=env, stdout=out, stderr=subprocess.STDOUT)
    log = (wd / "job.out").read_text()
    if proc.returncode != 0 or "normal termination of xtb" not in log:
        raise RuntimeError(f"xtb failed for {name}; see {wd / 'job.out'}")
    energy = float(re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", log)[-1])
    if opt:
        (wd / "xtbopt.xyz").replace(wd / f"{name}_opt.xyz")
    return energy, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help=".sdf/.mol/.xyz files (multi-structure files are fine)")
    ap.add_argument("--opt", action="store_true", help="geometry-optimize before reporting the energy")
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--uhf", type=int, default=0, help="number of unpaired electrons")
    ap.add_argument("--solvent", default=None, help="ALPB solvent name, e.g. water; gas phase if omitted")
    args = ap.parse_args()

    rows = []
    for f in args.files:
        structures = list(read_structures(f))
        for i, (_, atoms) in enumerate(structures):
            label = Path(f).stem if len(structures) == 1 else f"{Path(f).stem}_{i + 1}"
            e, secs = run_xtb(label, atoms, args.opt, args.nprocs, args.charge, args.uhf, args.solvent)
            print(f"  done {label}: {e:.6f} Ha ({secs:.1f} s)", flush=True)
            rows.append((label, "".join(sorted(s for s, _ in atoms)), e))

    # Relative energies only make sense between conformers of the same molecule (same atoms)
    e_min = {}
    for _, formula, e in rows:
        e_min[formula] = min(e, e_min.get(formula, e))
    width = max(len(r[0]) for r in rows)
    kind = "optimized" if args.opt else "single point"
    print(f"\nGFN2-xTB ({args.solvent or 'gas phase'}), {kind}")
    print(f"{'structure':{width}s}  E (Hartree)       rel. to lowest same-formula (kcal/mol)")
    for label, formula, e in rows:
        print(f"{label:{width}s}  {e:16.6f}  {(e - e_min[formula]) * HARTREE_TO_KCAL:8.2f}")


if __name__ == "__main__":
    main()
