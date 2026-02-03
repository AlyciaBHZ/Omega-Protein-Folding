#!/usr/bin/env python3
from __future__ import annotations

"""
Extract FASTA sequences for QUARK / I-TASSER benchmark targets from the local
Decoys 'R' Us native PDB cache.

Design goals:
- Reuse the repo's convention of selecting the longest CA-trace chain.
- Be robust to missing/unknown residues (emit 'X').
- Emit auditable metadata: chosen chain, length, source path.

Outputs (committable):
- docs/runs/quark_itasser_homology_ablation/targets_6.fasta
- docs/runs/quark_itasser_homology_ablation/targets.md
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from Bio.PDB.PDBParser import PDBParser


AA3_TO_1: Dict[str, str] = {
    # standard
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # common variants
    "MSE": "M",  # selenomethionine
    "SEC": "U",
    "PYL": "O",
}


TARGETS_6: Tuple[str, ...] = ("1R69", "2CRO", "4PTI", "1CTF", "1DTK", "1SHF-A")


@dataclass(frozen=True)
class TargetSeq:
    target_id: str
    dataset: str
    native_path: Path
    chain_id: str
    seq: str


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _iter_manifest_rows(path: Path) -> Iterable[dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            yield row


def _native_paths_from_manifests(root: Path) -> Dict[str, Tuple[str, Path]]:
    """
    Return map target_id -> (dataset, native_path).
    Uses the decoy manifests already tracked in-repo.
    """
    out: Dict[str, Tuple[str, Path]] = {}
    man_paths = [
        root / "docs" / "runs" / "decoys_4state_auric" / "decoy_manifest.csv",
        root / "docs" / "runs" / "decoys_lmds_auric" / "decoy_manifest.csv",
    ]
    for mp in man_paths:
        if not mp.exists():
            continue
        for row in _iter_manifest_rows(mp):
            tid = str(row.get("target_id", "")).strip()
            ds = str(row.get("dataset", "(unknown)")).strip()
            npth = str(row.get("native_path", "")).strip()
            if tid and npth and tid not in out:
                out[tid] = (ds, (root / npth).resolve())
    return out


def _extract_longest_ca_chain_seq(pdb_path: Path) -> Tuple[str, str]:
    """
    Return (chain_id, seq) where seq is built from residues that contain CA
    (mirrors coords_io longest-CA-chain choice).
    """
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(s.get_models()))

    best_chain_id = None
    best_seq: List[str] = []

    for chain in model:
        seq: List[str] = []
        for res in chain:
            # Skip hetero and waters
            if res.id[0].strip():
                continue
            if "CA" not in res:
                continue
            aa3 = str(res.get_resname()).upper().strip()
            seq.append(AA3_TO_1.get(aa3, "X"))
        if len(seq) > len(best_seq):
            best_seq = seq
            best_chain_id = str(chain.id)

    if best_chain_id is None or not best_seq:
        raise ValueError(f"No CA-trace residues found in {pdb_path}")
    return best_chain_id, "".join(best_seq)


def _wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def main() -> None:
    root = _repo_root()
    native_map = _native_paths_from_manifests(root)

    missing = [t for t in TARGETS_6 if t.lower() not in {k.lower() for k in native_map.keys()}]
    if missing:
        raise SystemExit(f"Missing native paths for targets in manifests: {missing}")

    # resolve by case-insensitive match
    def _resolve_tid(t: str) -> str:
        for k in native_map.keys():
            if k.lower() == t.lower():
                return k
        return t

    seqs: List[TargetSeq] = []
    for t in TARGETS_6:
        tid = _resolve_tid(t)
        ds, npth = native_map[tid]
        if not npth.exists():
            raise SystemExit(f"Native PDB not found: target={tid} path={npth}")
        chain_id, seq = _extract_longest_ca_chain_seq(npth)
        # Some PDBs may have blank/space chain IDs. Normalize to 'A' for server submission headers.
        chain_id_norm = chain_id.strip() if str(chain_id).strip() else "A"
        seqs.append(
            TargetSeq(target_id=tid, dataset=ds, native_path=npth, chain_id=chain_id_norm, seq=seq)
        )

    run_dir = root / "docs" / "runs" / "quark_itasser_homology_ablation"
    _ensure_dir(run_dir)

    out_fasta = run_dir / "targets_6.fasta"
    out_md = run_dir / "targets.md"

    # FASTA
    fasta_lines: List[str] = []
    for s in seqs:
        fasta_lines.append(f">{s.target_id}|chain={s.chain_id}|len={len(s.seq)}|dataset={s.dataset}")
        fasta_lines.append(_wrap_fasta(s.seq))
    out_fasta.write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")

    # Targets markdown (with placeholders for job IDs)
    md: List[str] = []
    md.append("# QUARK / I-TASSER homology-ablation targets (N=6)")
    md.append("")
    md.append("This folder supports the final benchmark described in the plan: QUARK vs I-TASSER,")
    md.append("and each with a homology/fragment exclusion ablation to assess template dependence.")
    md.append("")
    md.append("## Target list")
    md.append("")
    md.append("| target_id | dataset | chain | length | native_pdb | QUARK_A | QUARK_B | ITASSER_A | ITASSER_B | notes |")
    md.append("| --- | --- | --- | ---: | --- | --- | --- | --- | --- | --- |")
    for s in seqs:
        native_rel = s.native_path.relative_to(root).as_posix()
        md.append(
            f"| {s.target_id} | {s.dataset} | {s.chain_id} | {len(s.seq)} | `{native_rel}` |  |  |  |  |  |"
        )
    md.append("")
    md.append("## Inputs")
    md.append("")
    md.append(f"- FASTA (copy/paste into servers): `{out_fasta.relative_to(root).as_posix()}`")
    md.append("")
    md.append("## Submission protocol reminder")
    md.append("")
    md.append("- QUARK-A: default fragments")
    md.append("- QUARK-B: exclude fragments from proteins with >30% sequence identity")
    md.append("- I-TASSER-A: default")
    md.append("- I-TASSER-B: exclude homologous templates (benchmark option)")
    md.append("")
    out_md.write_text("\n".join(md) + "\n", encoding="utf-8")

    print(f"[ok] wrote: {out_fasta}")
    print(f"[ok] wrote: {out_md}")


if __name__ == "__main__":
    main()

