# Energy oracles: DLPNO-CCSD(T), r2SCAN-3c, GFN2-xTB

Which of two structures of the same molecule is lower in energy, scored as given (single point).
Source: single points on QM9 vs RDKit ETKDG + MMFF structures of 30 random QM9 molecules (2026-10-23).

| Role | Method | Error vs DLPNO, kcal/mol (mean / median / max) | Same lower structure as DLPNO | Cost per structure |
|---|---|---|---|---|
| Reference | DLPNO-CCSD(T)/cc-pVTZ | — | — | ~15 min (16 cores) |
| Better, slower oracle | r2SCAN-3c | 0.70 / 0.50 / 3.8 | 30/30 | ~11 s (16 cores) |
| Quick oracle | GFN2-xTB | 1.45 / 1.02 / 8.2 | 30/30 | ~0.03 s (1 core) |

- **DLPNO**: not exact (~0.5–1 kcal/mol uncertainty). Use for spot checks only.
- **r2SCAN-3c**: about 80× cheaper than DLPNO. Use for close calls.
- **xTB**: unreliable for gaps under ~2 kcal/mol (wrong on gdb_52425: DLPNO +0.26, xTB −0.40). Use for screening.
- **Workflow**: screen with xTB → re-score with r2SCAN-3c if the xTB gap is under 2 kcal/mol → spot-check with DLPNO.
- **Caveats**: electronic energies only, no zero-point or thermal corrections. Tested only on 30 small, neutral QM9 molecules. MMFF and the QM9-trained TensorNet are not usable as oracles.

Run: see `README.md`.
