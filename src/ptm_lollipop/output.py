from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .intervals import fmt_ranges
from .models import ProteinRecord


def write_tables(records: list[ProteinRecord], outdir: str | Path) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "disorder_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "category",
                "gene",
                "systematic",
                "length",
                "user_disorder",
                "sgd_mobidblite_raw",
                "sgd_mobidblite_collapsed",
                "qc_status",
                "qc_note",
            ],
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "category": record.category,
                    "gene": record.gene,
                    "systematic": record.systematic,
                    "length": record.length,
                    "user_disorder": record.user_disorder_text or "blank",
                    "sgd_mobidblite_raw": fmt_ranges(record.raw_disorder),
                    "sgd_mobidblite_collapsed": fmt_ranges(record.disorder),
                    "qc_status": record.qc_status,
                    "qc_note": record.qc_note,
                }
            )

    with (out / "ptm_sites.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["category", "gene", "systematic", "site", "residue", "ptm_family", "raw_sgd_types"],
        )
        writer.writeheader()
        for record in records:
            for ptm in record.ptms:
                writer.writerow(
                    {
                        "category": record.category,
                        "gene": record.gene,
                        "systematic": record.systematic,
                        "site": ptm.site,
                        "residue": ptm.residue,
                        "ptm_family": ptm.family,
                        "raw_sgd_types": "; ".join(ptm.raw_types),
                    }
                )

    serializable = [asdict(record) for record in records]
    (out / "records.json").write_text(json.dumps(serializable, indent=2), encoding="utf-8")

