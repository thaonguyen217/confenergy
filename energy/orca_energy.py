"""Quantum-chemistry energies with ORCA.

Usage:
    python energy/orca_energy.py examples/gdb_52425_qm9.sdf examples/gdb_52425_etkdg_mmff.sdf \
        [--method qm9|r2scan3c|b3lyp-d3|dlpno] [--opt] [--nprocs 16]

Methods:
    qm9       B3LYP/6-31G(2df,p), QM9's own level (default). Use to compare with QM9 labels.
              No dispersion correction and a small basis, so not the best for conformer energies.
    r2scan3c  r2SCAN-3c composite DFT (dispersion included). Cheap; good for conformer energies.
    b3lyp-d3  B3LYP-D3(BJ)/def2-TZVP. Slower, more accurate.
    dlpno     DLPNO-CCSD(T)/cc-pVTZ, near coupled-cluster quality. Single point only: run it on
              geometries already optimized with a DFT method.

Single point by default; --opt relaxes each structure first (needed when the input geometry is
not a minimum at that level, e.g. an MMFF conformer). Runs go to temp/orca_runs/<method>/<structure>/,
and an optimized geometry is written there as <structure>_opt.xyz.

For qm9: B3LYP/G is ORCA's Gaussian-compatible B3LYP (VWN3), which is what QM9 used. ORCA always
uses spherical d/f functions while Gaussian's 6-31G-family default is Cartesian d (6D), so absolute
energies are 0.5-5 mHa per heavy atom above QM9's u0 - zpve; compare relative energies computed here.
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

ORCA_HOME = Path(os.environ.get("ORCA_HOME", "/shared/nas2/thaotn2/orca-6.1.0-f.0_linux_x86-64"))
MPI_HOME = Path(os.environ.get("MPI_HOME", "/shared/nas2/thaotn2/anaconda3/envs/shm_env"))  # OpenMPI 4.1.6, as ORCA 6.1 requires
RUNS = Path(__file__).parent.parent / "temp" / "orca_runs"
HARTREE_TO_KCAL = 627.5095
# name -> (label, simple-input keywords, extra input block)
METHODS = {
    # ORCA has no built-in 6-31G(2df,p): use 6-31G(2df) for heavy atoms and 6-31G(d,p)'s hydrogen (one p shell)
    "qm9": ("B3LYP/6-31G(2df,p)", "B3LYP/G 6-31G(2df) TightSCF", '%basis\n  NewGTO H "6-31G(d,p)" end\nend\n'),
    "r2scan3c": ("r2SCAN-3c", "r2SCAN-3c TightSCF", ""),
    "b3lyp-d3": ("B3LYP-D3(BJ)/def2-TZVP", "B3LYP D3BJ def2-TZVP def2/J RIJCOSX TightSCF", ""),
    "dlpno": ("DLPNO-CCSD(T)/cc-pVTZ", "DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C TightPNO TightSCF", ""),
}


def run_orca(name, atoms, method, opt, nprocs, charge, mult, maxcore):
    _, simple, block = METHODS[method]
    wd = RUNS / method / name
    wd.mkdir(parents=True, exist_ok=True)
    write_xyz(atoms, wd / "input.xyz", comment=name)
    keywords = f"! {simple}" + (" TightOpt" if opt else "")
    (wd / "job.inp").write_text(
        f"{keywords}\n{block}%pal nprocs {nprocs} end\n%maxcore {maxcore}\n* xyzfile {charge} {mult} input.xyz\n"
    )
    env = dict(os.environ)
    env["PATH"] = f"{ORCA_HOME}/bin:{ORCA_HOME}:{MPI_HOME}/bin:{env.get('PATH', '')}"
    env["LD_LIBRARY_PATH"] = f"{ORCA_HOME}/lib:{MPI_HOME}/lib:{env.get('LD_LIBRARY_PATH', '')}"
    t0 = time.time()
    with open(wd / "job.out", "w") as out:
        # ORCA must be called by its full path for parallel runs
        proc = subprocess.run([str(ORCA_HOME / "bin" / "orca"), "job.inp"], cwd=wd, env=env,
                              stdout=out, stderr=subprocess.STDOUT)
    log = (wd / "job.out").read_text()
    if proc.returncode != 0 or "ORCA TERMINATED NORMALLY" not in log:
        raise RuntimeError(f"ORCA failed for {name}; see {wd / 'job.out'}")
    energy = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", log)[-1])
    if opt:
        (wd / "job.xyz").replace(wd / f"{name}_opt.xyz")
    return energy, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help=".sdf/.mol/.xyz files (multi-structure files are fine)")
    ap.add_argument("--method", choices=METHODS, default="qm9")
    ap.add_argument("--opt", action="store_true", help="geometry-optimize before reporting the energy")
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--maxcore", type=int, default=2000, help="MB of memory per process")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    args = ap.parse_args()
    if args.opt and args.method == "dlpno":
        ap.error("dlpno has no analytic gradients here; optimize with a DFT method, then run a dlpno single point")

    rows = []
    for f in args.files:
        structures = list(read_structures(f))
        for i, (_, atoms) in enumerate(structures):
            label = Path(f).stem if len(structures) == 1 else f"{Path(f).stem}_{i + 1}"
            e, secs = run_orca(label, atoms, args.method, args.opt, args.nprocs, args.charge, args.mult, args.maxcore)
            print(f"  done {label}: {e:.6f} Ha ({secs:.0f} s)", flush=True)
            rows.append((label, "".join(sorted(s for s, _ in atoms)), e))

    # Relative energies only make sense between conformers of the same molecule (same atoms)
    e_min = {}
    for _, formula, e in rows:
        e_min[formula] = min(e, e_min.get(formula, e))
    width = max(len(r[0]) for r in rows)
    kind = "optimized" if args.opt else "single point"
    print(f"\n{METHODS[args.method][0]}, {kind}")
    print(f"{'structure':{width}s}  E (Hartree)       rel. to lowest same-formula (kcal/mol)")
    for label, formula, e in rows:
        print(f"{label:{width}s}  {e:16.6f}  {(e - e_min[formula]) * HARTREE_TO_KCAL:8.2f}")


if __name__ == "__main__":
    main()
