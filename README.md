# confenergy

Tools to compute energies of 3D molecular structures (conformers) with GFN2-xTB and ORCA.

| Method | Command | Cost per QM9-sized molecule |
|---|---|---|
| GFN2-xTB | `python energy/xtb_energy.py <files>` | ~0.03 s, 1 core |
| r2SCAN-3c | `python energy/orca_energy.py <files> --method r2scan3c` | ~11 s, 16 cores |
| B3LYP/6-31G(2df,p) (QM9's level) | `... --method qm9` | ~19 s, 16 cores |
| B3LYP-D3(BJ)/def2-TZVP | `... --method b3lyp-d3` | ~24 s, 16 cores |
| DLPNO-CCSD(T)/cc-pVTZ | `... --method dlpno` | ~15 min, 16 cores |

- Input: `.sdf`, `.mol` or `.xyz` files (multi-structure files are fine).
- Output: absolute energies, plus energies relative to the lowest structure with the same formula.
- Single points by default. Add `--opt` (xTB and ORCA DFT methods only, not DLPNO) to optimize first.
- Other flags: `--nprocs`, `--charge`, `--mult` (ORCA) / `--uhf` and `--solvent` (xTB).
- Run directories go to `temp/xtb_runs/` and `temp/orca_runs/<method>/`.
- Which method to trust: `notes/energy_oracles.md`.

Example:

```bash
python energy/xtb_energy.py examples/gdb_52425_qm9.sdf examples/gdb_52425_etkdg_mmff.sdf
python energy/orca_energy.py examples/gdb_52425_qm9.sdf examples/gdb_52425_etkdg_mmff.sdf --method r2scan3c --nprocs 16
```

## Requirements

- [ORCA](https://www.faccts.de/orca/) 6.1 with OpenMPI 4.1.x (see below).
- [xTB](https://github.com/grimme-lab/xtb) 6.7 (see below).
- Python 3 standard library only.

## Installing ORCA

ORCA is free for academic use but cannot be redistributed, so each user downloads it after registering.

1. Register at the [ORCA Forum](https://orcaforum.kofo.mpg.de) and accept the academic license (approval can take a day).
2. In the forum's **Downloads** section, get the Linux x86-64 build of ORCA 6.1 (the archive linked against OpenMPI 4.1) and unpack it:
   ```bash
   tar -xf orca_6_1_*_linux_x86-64*.tar.xz -C ~/software    # adjust to the downloaded file name
   ```
3. ORCA's parallel runs need OpenMPI 4.1 (not 5.x). With conda:
   ```bash
   conda create -n orca_mpi -c conda-forge openmpi=4.1
   ```
4. Point the scripts at both:
   ```bash
   export ORCA_HOME=~/software/<unpacked ORCA folder>     # folder containing bin/orca
   export MPI_HOME=$(conda info --base)/envs/orca_mpi    # folder containing bin/mpirun
   ```
   `orca_energy.py` adds both to `PATH` and `LD_LIBRARY_PATH` and calls ORCA by its full path, as ORCA requires for parallel runs.
5. Check it:
   ```bash
   python energy/orca_energy.py examples/gdb_52425_qm9.sdf --method r2scan3c --nprocs 4
   ```

Use `--nprocs` up to the number of free cores, with about 2 GB of memory per process (`--maxcore`, in MB, default 2000).

## Installing xTB

```bash
conda create -n xtb -c conda-forge xtb=6.7
export XTB_BIN=$(conda info --base)/envs/xtb/bin/xtb
```
