from __future__ import annotations

import re
from dataclasses import dataclass, replace


STANDARD_AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
STANDARD_AMINO_ACID_SET = frozenset(STANDARD_AMINO_ACIDS)

_MUTATION_RE = re.compile(
    r"^(?:(?P<chain>[^:_\s]+):)?"
    r"(?P<wt>[A-Z])(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)(?P<target>[A-Z])$"
)


@dataclass(frozen=True, order=True)
class Mutation:
    chain: str | None
    wt: str
    residue_number: int
    insertion_code: str
    target: str

    @property
    def residue_key(self) -> tuple[str | None, int, str]:
        return self.chain, self.residue_number, self.insertion_code

    @property
    def canonical(self) -> str:
        chain = f"{self.chain}:" if self.chain else ""
        return f"{chain}{self.wt}{self.residue_number}{self.insertion_code}{self.target}"

    def with_chain(self, chain: str) -> "Mutation":
        return replace(self, chain=chain)


def parse_mutation(token: str) -> Mutation:
    text = token.strip().upper()
    match = _MUTATION_RE.fullmatch(text)
    if not match:
        raise ValueError(f"invalid mutation token: {token!r}")
    wt = match.group("wt")
    target = match.group("target")
    if wt not in STANDARD_AMINO_ACID_SET or target not in STANDARD_AMINO_ACID_SET:
        raise ValueError(f"mutation uses a non-standard amino-acid code: {token!r}")
    if wt == target:
        raise ValueError(f"mutation target equals WT: {token!r}")
    return Mutation(
        chain=match.group("chain"),
        wt=wt,
        residue_number=int(match.group("resnum")),
        insertion_code=match.group("icode").upper(),
        target=target,
    )


def parse_variant(text: str) -> tuple[Mutation, ...]:
    tokens = [token for token in text.strip().split("_") if token]
    if not tokens:
        raise ValueError("variant has no mutation tokens")
    mutations = tuple(parse_mutation(token) for token in tokens)
    keys = [mutation.residue_key for mutation in mutations]
    if len(keys) != len(set(keys)):
        raise ValueError(f"variant mutates the same residue more than once: {text!r}")
    return mutations


def canonical_variant(mutations: tuple[Mutation, ...]) -> str:
    return "_".join(mutation.canonical for mutation in mutations)


def variant_id_from_name(name: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")
    if not normalized:
        raise ValueError(f"cannot normalize variant identifier from {name!r}")
    return normalized
