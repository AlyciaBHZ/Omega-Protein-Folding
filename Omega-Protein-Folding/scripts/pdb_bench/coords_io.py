from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser


def _longest_chain_ca_coords(model) -> Tuple[str, np.ndarray]:
    best_chain = None
    best_coords = None
    for chain in model:
        coords = []
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].get_coord())
        if len(coords) >= 2:
            arr = np.asarray(coords, dtype=np.float64)
            if best_coords is None or len(arr) > len(best_coords):
                best_coords = arr
                best_chain = chain.id
    if best_chain is None or best_coords is None:
        raise ValueError("No CA trace found")
    return str(best_chain), best_coords


def load_ca_coords(path: Path) -> Tuple[str, np.ndarray]:
    """
    Load a CA trace from .pdb or .cif/.mmcif, returning (chain_id, coords[N,3])
    for the longest chain in the first model.
    """
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".pdb":
        parser = PDBParser(QUIET=True)
        s = parser.get_structure(p.stem, str(p))
    elif suf in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure(p.stem, str(p))
    else:
        raise ValueError(f"Unsupported structure suffix: {p.suffix} (expected .pdb/.cif/.mmcif)")

    model = next(iter(s.get_models()))
    return _longest_chain_ca_coords(model)

