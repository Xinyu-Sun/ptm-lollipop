from __future__ import annotations

import csv
from pathlib import Path

from .models import ProteinInput


def _sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8-sig").splitlines()[:5]
    joined = "\n".join(sample)
    if "\t" in joined:
        return "\t"
    return ","


def read_protein_inputs(path: str | Path) -> list[ProteinInput]:
    input_path = Path(path)
    delimiter = _sniff_delimiter(input_path)
    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or "gene" not in {f.strip() for f in reader.fieldnames}:
            raise ValueError("Input must include a 'gene' column.")

        rows: list[ProteinInput] = []
        for i, raw in enumerate(reader, start=2):
            row = {(k or "").strip(): (v or "").strip() for k, v in raw.items()}
            gene = row.get("gene", "")
            if not gene:
                raise ValueError(f"Missing gene value on input line {i}.")
            rows.append(
                ProteinInput(
                    gene=gene,
                    category=row.get("category") or "Proteins",
                    label=row.get("label") or gene,
                    disorder_text=row.get("disorder") or "",
                )
            )
    return rows

