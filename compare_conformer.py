"""Compare a QM9 DFT conformer with an RDKit ETKDG + MMFF conformer for the same molecule.

Usage: python compare_conformer.py [mol_id]   (e.g. gdb_100000; picks one at random if omitted)
Writes temp/out/<mol_id>_qm9.sdf, <mol_id>_etkdg_mmff.sdf, <mol_id>_overlay.json and an interactive
temp/out/<mol_id>_overlay.html (3Dmol.js) showing both conformers in the same frame.
"""
import json
import os
import random
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors, rdMolTransforms

RDLogger.DisableLog("rdApp.*")

DATA = Path(os.environ.get("QM9_DIR", Path(__file__).parent / "data" / "qm9"))  # gdb9.sdf, uncharacterized.txt
OUT = Path(__file__).parent / "temp" / "out"
SEED = 42


def uncharacterized_ids():
    # uncharacterized.txt: header lines, then rows whose first column is the integer index
    ids = set()
    for line in (DATA / "uncharacterized.txt").read_text().splitlines():
        parts = line.split()
        if parts and parts[0].isdigit():
            ids.add(f"gdb_{parts[0]}")
    return ids


def load_qm9(mol_id=None):
    bad = uncharacterized_ids()
    suppl = Chem.SDMolSupplier(str(DATA / "gdb9.sdf"), removeHs=False, sanitize=True)
    if mol_id is not None:
        idx = int(mol_id.split("_")[1]) - 1
        mol = suppl[idx]
        assert mol is not None and mol.GetProp("_Name") == mol_id, f"could not load {mol_id}"
        return mol
    # Random pick: a characterized 9-heavy-atom molecule with a ring and a rotatable bond,
    # so the comparison is more interesting than a rigid cage.
    rng = random.Random(SEED)
    while True:
        mol = suppl[rng.randrange(len(suppl))]
        if mol is None or mol.GetProp("_Name") in bad:
            continue
        if (mol.GetNumHeavyAtoms() == 9 and mol.GetRingInfo().NumRings() >= 1
                and rdMolDescriptors.CalcNumRotatableBonds(mol) >= 1):
            return mol


def etkdg_mmff(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        # Strained cages often fail the default eigenvalue-based start; random coordinates are RDKit's standard fallback
        params.useRandomCoords = True
        assert AllChem.EmbedMolecule(mol, params) == 0, "ETKDG embedding failed"
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    Chem.SanitizeMol(mol)  # MMFF setup kekulizes; restore aromatic flags so it matches the QM9 graph
    return mol


def main():
    ref = load_qm9(sys.argv[1] if len(sys.argv) > 1 else None)
    mol_id = ref.GetProp("_Name")
    smiles = Chem.MolToSmiles(Chem.RemoveHs(ref))
    gen = etkdg_mmff(smiles)

    # Symmetry-aware best alignment of the generated conformer onto the QM9 one
    rmsd_heavy = rdMolAlign.GetBestRMS(Chem.RemoveHs(gen), Chem.RemoveHs(ref))
    rmsd_all, transform, atom_map = rdMolAlign.GetBestAlignmentTransform(gen, ref)
    rdMolTransforms.TransformConformer(gen.GetConformer(), transform)

    # Per-atom displacement, in QM9 atom order, using the symmetry-matched atom mapping
    ref_pos, gen_pos = ref.GetConformer().GetPositions(), gen.GetConformer().GetPositions()
    probe_of = {r: p for p, r in atom_map}
    displacements = [
        {"elem": ref.GetAtomWithIdx(r).GetSymbol(),
         "d": round(float(np.linalg.norm(ref_pos[r] - gen_pos[probe_of[r]])), 4)}
        for r in range(ref.GetNumAtoms())
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    Chem.MolToMolFile(ref, str(OUT / f"{mol_id}_qm9.sdf"))
    gen.SetProp("_Name", f"{mol_id}_etkdg_mmff")
    Chem.MolToMolFile(gen, str(OUT / f"{mol_id}_etkdg_mmff.sdf"))
    result = {
        "mol_id": mol_id,
        "smiles": smiles,
        "rmsd_heavy": round(rmsd_heavy, 3),
        "rmsd_all": round(rmsd_all, 3),
        "qm9_molblock": Chem.MolToMolBlock(ref),
        "gen_molblock": Chem.MolToMolBlock(gen),
        "displacements": displacements,
    }
    (OUT / f"{mol_id}_overlay.json").write_text(json.dumps(result))
    template = (Path(__file__).parent / "overlay_template.html").read_text()
    (OUT / f"{mol_id}_overlay.html").write_text(template.replace("/*DATA*/null", json.dumps(result)).replace("__MOL_ID__", mol_id))
    print(f"{mol_id}  {smiles}  heavy-atom RMSD {rmsd_heavy:.3f} Å  all-atom RMSD {rmsd_all:.3f} Å")


if __name__ == "__main__":
    main()
