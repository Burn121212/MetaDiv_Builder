# MetaDiv_metadata

MetaDiv_metadata assigns reference metadata to a `sample_metadata_to_fill.csv` obtained with MetaDiv Bulder MAIN SCRIPT in the For_R/ diretory, using a the sample_IDs of a reference sample file `reference_metadata.csv`

## Input

Place both files inside:

source_metadata/

- sample_metadata_to_fill.csv
- reference_metadata.csv

Edit the first notebook cell with the exact file names.

## Run

Run the second cell.

The script automatically matches samples using the first column of each table.

## Output

Results are saved in:

merged_metadata/

- sample_metadata.csv
- unmatched_samples.csv
- duplicated_reference_IDs.csv

`duplicated_reference_IDs.csv` is created only if duplicated reference IDs are detected.
