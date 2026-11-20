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

- [ORCA](https://www.faccts.de/orca/) 6.1 with OpenMPI 4.1.x. Set `ORCA_HOME` (ORCA folder) and `MPI_HOME` (folder containing `bin/mpirun`).
- [xTB](https://github.com/grimme-lab/xtb) 6.7, e.g. `conda install -c conda-forge xtb`. Set `XTB_BIN` to the `xtb` executable.
- Python 3 standard library only.
