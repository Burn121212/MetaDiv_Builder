
# MetaDiv_KRONA

**MetaDiv_KRONA** is a standalone utility for generating interactive Krona-style hierarchical charts from MetaDiv-compatible CSV tables.

It is designed for visualizing taxonomic and ecological subsets extracted from a MetaDiv `FINAL_DB`, but also supports other compatible CSV tables.

## Main Functions

- **Run MetaDiv_KRONA:** Generates interactive taxonomic charts from global or individual sample abundances.
- **KRONA by sample metadata:** Groups samples by metadata categories (e.g., state, ecoregion, vegetation) to visualize taxonomic composition.
- **KRONA richness:** Visualizes species richness per genus across metadata groups, based on unique species rather than read abundance.

## Integration with MetaDiv Subsets

CSV files generated using **MetaDiv Subsets** can be copied directly into `input_KRONA/`. Running **MetaDiv_KRONA** automatically generates interactive HTML charts for each subset, enabling rapid taxonomic exploration.

## Output

The notebook generates:

- Krona-compatible tab-delimited `.txt` tables
- Self-contained interactive `.html` visualizations
- Supporting `.csv` tables

HTML files can be opened locally in any modern web browser.

## Requirements

No additional Krona-specific software is required.

- Docker
- An internet connection
