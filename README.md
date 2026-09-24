# PTM Lollipop

Reusable PTM lollipop plots for yeast proteins using Saccharomyces Genome Database (SGD) annotations.

The repo turns a simple protein list into publication-ready linear protein schematics:

- protein backbones scaled by amino-acid length
- PTM sites pulled from SGD and deduplicated by site, residue, and PTM family
- optional disorder intervals from SGD MobiDB-lite
- optional user-supplied disorder coordinates checked against SGD
- category-level PNG/PDF panels and QC tables

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

ptm-lollipop examples/alan_table_s1_minimal.tsv --outdir outputs
```

The main input is a TSV or CSV with these columns:

```text
category	gene	label	disorder
Axial Proteins	CDC14	Cdc14 E	383-516
Axial Proteins	FOB1	Fob1	1-44; 449-566
```

Required columns:

- `gene`: SGD gene name or identifier used in SGD `/backend/locus/{gene}`.

Optional columns:

- `category`: used to split output panels. Defaults to `Proteins`.
- `label`: display label. Defaults to `gene`.
- `disorder`: user-supplied disorder ranges for QC, such as `1-44; 449-566`. The plot itself uses SGD MobiDB-lite by default.

## Outputs

For each category, the CLI writes:

- `{category}_ptm_lollipop.pdf`
- `{category}_ptm_lollipop.png`

It also writes:

- `ptm_sites.tsv`: deduplicated PTM annotations used for plotting
- `disorder_qc.tsv`: user disorder coordinates compared with SGD MobiDB-lite
- `records.json`: normalized data used by the plotting layer

## Data Source

The default fetcher uses SGD backend endpoints:

- `/backend/locus/{gene}`
- `/backend/locus/{gene}/sequence_details`
- `/backend/locus/{gene}/protein_domain_details`
- `/backend/locus/{gene}/posttranslational_details`

Responses are cached in `cache/` by default. Keep the cache for reproducibility during a figure build, but do not commit bulky cache/output files unless you intentionally want a frozen snapshot.

## Development

```bash
pip install -e ".[dev]"
pytest
```

The tests use local fixtures and do not hit SGD.

