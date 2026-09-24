from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ProteinInput:
    gene: str
    category: str = "Proteins"
    label: str = ""
    disorder_text: str = ""


@dataclass(frozen=True)
class PtmSite:
    site: int
    residue: str
    family: str
    raw_types: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class ProteinRecord:
    category: str
    query: str
    gene: str
    systematic: str
    sgdid: str
    uniprot: str
    label: str
    length: int
    user_disorder_text: str
    user_disorder: tuple[tuple[int, int], ...]
    raw_disorder: tuple[tuple[int, int], ...]
    disorder: tuple[tuple[int, int], ...]
    ptms: tuple[PtmSite, ...]
    qc_note: str
    qc_status: str

