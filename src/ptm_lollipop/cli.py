from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .input import read_protein_inputs
from .output import write_tables
from .plot import plot_all_categories
from .sgd import SgdClient


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate PTM lollipop plots from SGD annotations.")
    parser.add_argument("input", help="Input TSV/CSV with at least a gene column.")
    parser.add_argument("--outdir", default="outputs", help="Output directory. Default: outputs")
    parser.add_argument("--cache-dir", default="cache", help="SGD response cache directory. Default: cache")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution. Default: 300")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    inputs = read_protein_inputs(args.input)
    client = SgdClient(args.cache_dir)

    records = []
    for index, item in enumerate(inputs, start=1):
        print(f"[{index}/{len(inputs)}] {item.gene}", file=sys.stderr, flush=True)
        records.append(client.protein_record(item))

    outdir = Path(args.outdir)
    write_tables(records, outdir)
    written = plot_all_categories(records, outdir, dpi=args.dpi)

    print("Wrote:")
    table_paths = [outdir / "ptm_sites.tsv", outdir / "sgd_disorder.tsv", outdir / "records.json"]
    if any(record.user_disorder_text for record in records):
        table_paths.append(outdir / "disorder_qc.tsv")
    for path in written + table_paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
