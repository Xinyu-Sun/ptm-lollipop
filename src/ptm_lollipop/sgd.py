from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from .intervals import compare_disorder, merge_intervals, parse_ranges
from .models import ProteinInput, ProteinRecord, PtmSite

SGD_BASE = "https://www.yeastgenome.org/backend"


def normalize_ptm_type(value: str) -> str:
    text = (value or "").lower()
    if "phosph" in text:
        return "phosphorylation"
    if "sumoy" in text or "sumo" in text:
        return "sumoylation"
    if "ubiquit" in text:
        return "ubiquitination"
    if "acetyl" in text:
        return "acetylation"
    if "methyl" in text:
        return "methylation"
    if "succin" in text:
        return "succinylation"
    if "glutathion" in text:
        return "glutathionylation"
    return "other"


class SgdClient:
    def __init__(self, cache_dir: str | Path, base_url: str = SGD_BASE, delay_s: float = 0.12):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.base_url = base_url.rstrip("/")
        self.delay_s = delay_s

    def json(self, identifier: str, suffix: str | None = None):
        name = identifier.upper()
        cache_name = f"{name}.{suffix or 'locus'}.json"
        path = self.cache_dir / cache_name
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))

        url = f"{self.base_url}/locus/{identifier}"
        if suffix:
            url += f"/{suffix}"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=20) as response:
                    data = response.read().decode("utf-8")
                path.write_text(data, encoding="utf-8")
                time.sleep(self.delay_s)
                return json.loads(data)
            except urllib.error.HTTPError as err:
                if err.code == 429 and attempt < 3:
                    time.sleep(2 + attempt)
                    continue
                raise

    def protein_record(self, item: ProteinInput) -> ProteinRecord:
        locus = self.json(item.gene)
        sequence = self.json(item.gene, "sequence_details")
        domains = self.json(item.gene, "protein_domain_details")
        ptm_rows = self.json(item.gene, "posttranslational_details")

        protein_seq = _choose_s288c_protein(sequence)
        residues = (protein_seq.get("residues") or "").rstrip("*")
        length = int(protein_seq.get("protein_length") or len(residues))

        raw_disorder: list[tuple[int, int]] = []
        for row in domains:
            source = row.get("source") or {}
            domain = row.get("domain") or {}
            source_name = source.get("format_name") or source.get("display_name")
            if source_name == "MobiDBLite" or domain.get("display_name") == "MobiDBLite":
                raw_disorder.append((int(row["start"]), int(row["end"])))
        disorder = merge_intervals(raw_disorder)

        by_key: dict[tuple[int, str, str], dict[str, object]] = {}
        for ptm in ptm_rows:
            site = ptm.get("site_index")
            if not site:
                continue
            residue = ptm.get("site_residue") or ""
            family = normalize_ptm_type(ptm.get("type") or "")
            key = (int(site), residue, family)
            by_key.setdefault(key, {"site": int(site), "residue": residue, "family": family, "raw_types": set()})
            by_key[key]["raw_types"].add(ptm.get("type") or "")

        ptms = tuple(
            PtmSite(
                site=value["site"],
                residue=value["residue"],
                family=value["family"],
                raw_types=tuple(sorted(value["raw_types"])),
            )
            for value in sorted(by_key.values(), key=lambda x: (x["site"], x["family"]))
        )

        qc_note, qc_status = compare_disorder(item.disorder_text, disorder, length)
        return ProteinRecord(
            category=item.category,
            query=item.gene,
            gene=locus.get("display_name") or item.gene,
            systematic=locus.get("format_name") or "",
            sgdid=locus.get("sgdid") or "",
            uniprot=locus.get("uniprot_id") or "",
            label=item.label or item.gene,
            length=length,
            user_disorder_text=item.disorder_text,
            user_disorder=parse_ranges(item.disorder_text, length),
            raw_disorder=tuple(raw_disorder),
            disorder=disorder,
            ptms=ptms,
            qc_note=qc_note,
            qc_status=qc_status,
        )


def _choose_s288c_protein(sequence_payload: dict) -> dict:
    proteins = sequence_payload.get("protein") or []
    for row in proteins:
        strain = row.get("strain") or {}
        if strain.get("format_name") == "S288C" or strain.get("display_name") == "S288C":
            return row
    if proteins:
        return proteins[0]
    raise ValueError("SGD sequence_details response did not include a protein sequence.")

