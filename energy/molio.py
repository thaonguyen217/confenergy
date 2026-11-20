"""Minimal SDF/XYZ reader shared by the energy tools.

No RDKit dependency: the tools only need element symbols and Cartesian coordinates (Å).
"""
from pathlib import Path

SYMBOLS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
           "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
           "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
           "Ga", "Ge", "As", "Se", "Br", "Kr"]
Z_OF = {s: i + 1 for i, s in enumerate(SYMBOLS)}


def _parse_molblock(lines):
    name = lines[0].strip()
    n_atoms = int(lines[3][:3])
    atoms = []
    for line in lines[4:4 + n_atoms]:
        x, y, z, sym = float(line[0:10]), float(line[10:20]), float(line[20:30]), line[31:34].strip()
        atoms.append((sym, (x, y, z)))
    return name, atoms


def read_structures(path):
    """Yield (name, [(symbol, (x, y, z)), ...]) for every structure in an .sdf/.mol or .xyz file."""
    path = Path(path)
    text = path.read_text()
    if path.suffix.lower() == ".xyz":
        lines = text.splitlines()
        i = 0
        while i < len(lines) and lines[i].strip():
            n = int(lines[i])
            name = lines[i + 1].strip() or path.stem
            atoms = [(p[0], tuple(map(float, p[1:4])))
                     for p in (l.split() for l in lines[i + 2:i + 2 + n])]
            yield name, atoms
            i += 2 + n
        return
    for i, block in enumerate(text.split("$$$$")):
        # Drop only the newline that follows $$$$: the molblock's first line (its name) may be blank
        lines = (block.removeprefix("\n") if i else block).splitlines()
        if len(lines) >= 4 and "V2000" in lines[3]:
            name, atoms = _parse_molblock(lines)
            yield name or path.stem, atoms


def write_xyz(atoms, path, comment=""):
    lines = [str(len(atoms)), comment]
    lines += [f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}" for s, (x, y, z) in atoms]
    Path(path).write_text("\n".join(lines) + "\n")
