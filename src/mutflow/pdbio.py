from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HID": "H",
    "HIE": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
}


class PDBError(ValueError):
    """Raised when a PDB cannot provide an unambiguous protein residue map."""


@dataclass(frozen=True)
class ResidueIdentity:
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    one_letter: str

    @property
    def key(self) -> tuple[str, int, str]:
        return self.chain, self.residue_number, self.insertion_code


def read_protein_residues(path: Path) -> dict[tuple[str, int, str], ResidueIdentity]:
    residues: dict[tuple[str, int, str], ResidueIdentity] = {}
    atom_records = 0
    try:
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise PDBError(f"cannot read PDB: {exc}") from exc

    for line_number, line in enumerate(lines, start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_records += 1
        if len(line) < 27:
            raise PDBError(f"short atom record at line {line_number}")
        residue_name = line[17:20].strip().upper()
        if residue_name not in AA3_TO_1:
            continue
        chain = line[21:22].strip()
        insertion_code = line[26:27].strip().upper()
        try:
            residue_number = int(line[22:26])
        except ValueError as exc:
            raise PDBError(f"invalid residue number at line {line_number}") from exc
        identity = ResidueIdentity(
            chain=chain,
            residue_number=residue_number,
            insertion_code=insertion_code,
            residue_name=residue_name,
            one_letter=AA3_TO_1[residue_name],
        )
        existing = residues.get(identity.key)
        if existing and existing.one_letter != identity.one_letter:
            raise PDBError(
                f"conflicting residue names for chain {chain!r} residue "
                f"{residue_number}{insertion_code}: {existing.residue_name}/{residue_name}"
            )
        residues[identity.key] = identity

    if atom_records == 0:
        raise PDBError("no ATOM/HETATM records found")
    if not residues:
        raise PDBError("no standard protein residues found")
    return residues
