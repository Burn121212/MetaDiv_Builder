# MetaDiv_KRONA

**MetaDiv_KRONA** is a standalone utility for generating interactive, offline Krona charts from MetaDiv-compatible CSV tables, including subsets exported by **MetaDiv Subsets**.

## Main functions

- **Abundance:** Visualizes global abundance, individual samples, or both.
- **Metadata groups:** Aggregates sample abundances by metadata categories.
- **Richness:** Counts unique species per genus within metadata groups. Global richness is the union of species across mapped samples.

## Python CLI for Windows, Linux and macOS

Requires **Python 3.9+**, `pandas`, and `numpy`.

No Jupyter, Docker, Excel, KronaTools, or `ktImportText` is required. Charts can be generated offline after installing dependencies.

### 1. Install dependencies

Using an existing MetaDiv environment:

```bash
conda activate metadiv
python -m pip install pandas numpy
```

Alternatively, create a dedicated environment:

```bash
conda create -n metadiv_krona python=3.11 pandas numpy
conda activate metadiv_krona
```

### 2. Prepare input files

Run the following commands from the folder containing `MetaDiv_KRONA.py`.

Create `input_KRONA/` and copy your MetaDiv-compatible CSV tables or subsets into it.

Input tables require taxonomic columns and numeric sample abundance columns. Richness requires `genus`, `species`, and sample abundances.

For grouping and richness, create `input_KRONA/metadata_gropus.txt`:

```text
sampleID,estado
SAMPLE_001,Oaxaca
SAMPLE_002,Oaxaca
SAMPLE_003,Puebla
```

The first column contains sample IDs matching the abundance column headers. The second contains the grouping category.

Tab, comma, and semicolon delimiters are supported. The filename `metadata_gropus.txt` preserves the notebook's original spelling.

### 3. Display options

```bash
python MetaDiv_KRONA.py --help
```

### 4. Run examples

Process all CSVs, generating individual sample datasets plus global abundance:

```bash
python MetaDiv_KRONA.py
```

Process one subset using global abundance:

```bash
python MetaDiv_KRONA.py --input-file subset_fungi.csv --mode global --taxonomic-cutoff genus
```

Generate abundance and metadata-grouped charts:

```bash
python MetaDiv_KRONA.py --analysis groups --mode both
```

Generate abundance, metadata groups, and richness:

```bash
python MetaDiv_KRONA.py --analysis all --mode both --taxonomic-cutoff species
```

Richness also works with a genus cutoff by recovering species identities from the original CSV:

```bash
python MetaDiv_KRONA.py --analysis richness --taxonomic-cutoff genus
```

Use custom folders and metadata:

```bash
python MetaDiv_KRONA.py --input-dir subsets --output-dir results_KRONA --metadata metadata/estados.csv --analysis all
```

Require every sample to be present in metadata:

```bash
python MetaDiv_KRONA.py --analysis all --strict-metadata
```

### Run from the MetaDiv base directory

If the script is stored in `MetaDiv_Utilies/MetaDiv_KRONA/`:

```bash
python MetaDiv_Utilies/MetaDiv_KRONA/MetaDiv_KRONA.py --input-dir MetaDiv_Utilies/MetaDiv_KRONA/input_KRONA --output-dir MetaDiv_Utilies/MetaDiv_KRONA/output_KRONA --analysis all
```

Adjust folder names to your repository. Paths are relative to your current working directory. Quote paths containing spaces.

## Options

| Option | Default / behavior |
| --- | --- |
| `--input-dir PATH` | `input_KRONA` |
| `--output-dir PATH` | `output_KRONA` |
| `--input-file NAME.csv` | Process one CSV; omitted = all CSVs in the input directory |
| `--mode MODE` | `global`, `samples`, or `both`; default: `both` |
| `--taxonomic-cutoff RANK` | Default: `genus`; `auto` uses the deepest available rank |
| `--analysis ANALYSIS` | `abundance`, `groups`, `richness`, or `all`; default: `abundance` |
| `--metadata PATH` | Default: `input-dir/metadata_gropus.txt` |
| `--strict-metadata` | Stop if any sample is missing from metadata |
| `--no-homogenize-taxonomy` | Disable lineage homogenization |
| `--include-primary-lifestyle` | Append primary lifestyle when available |
| `--include-secondary-lifestyle` | Append secondary lifestyle when available |
| `--include-sppn` | Append SPPN when available |
| `--max-krona-paths N` | Limit abundance paths; default: unlimited |
| `--no-grouped-global` | Omit global from grouped abundance |
| `--no-richness-global` | Omit global from richness |
| `--include-unclassified-species` | Include unclassified species labels in richness |
| `--genus-only-richness` | Omit ranks above genus in richness |

`--analysis groups` runs abundance followed by grouping.

`--analysis richness` and `--analysis all` run all three stages.

Grouping and richness require `--mode samples` or `--mode both`.

## Output

Results are saved under `output_KRONA/<input CSV stem>/`:

- Interactive, self-contained `.html` charts
- Tab-delimited `.txt` tables
- Supporting `.csv` tables
- Mapping and richness reports

Open HTML files locally in a browser. Use the Dataset menu to switch between samples or groups.

## Notes

- Metadata acts as the master list for grouped analyses. Unmapped samples are excluded from grouped abundance and richness and reported.
- Global grouped abundance sums mapped samples; global richness counts their unique species once.
- Commands rerun prerequisite stages and may overwrite matching outputs. Use another output directory to preserve previous runs.
- Leave `--max-krona-paths` unset for complete richness.
- Richness counts unique species labels in the input.
- Avoid unrelated numeric metadata columns in abundance CSVs because they may be detected as samples.
