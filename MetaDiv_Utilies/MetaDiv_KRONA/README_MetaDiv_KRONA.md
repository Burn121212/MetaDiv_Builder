# MetaDiv_KRONA

**MetaDiv_KRONA** is a standalone utility for converting MetaDiv-compatible CSV tables into interactive Krona-style hierarchical charts.

It is designed primarily for visualizing taxonomic or ecological subsets extracted from a MetaDiv `FINAL_DB`, but it can also be used with other compatible CSV tables containing taxonomic ranks and abundance data.

## Output

The notebook generates:

- Krona-compatible tab-delimited `.txt` tables
- Self-contained interactive `.html` visualizations

The HTML files are generated directly from Python and can be opened locally in any modern web browser.

## Requirements

No additional Krona-specific software is required to generate the HTML charts.

MetaDiv_KRONA does **not** require:

- KronaTools
- `ktImportText`
- Excel
- Docker
- An internet connection
