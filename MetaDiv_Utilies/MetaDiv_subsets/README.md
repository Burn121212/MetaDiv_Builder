# MetaDiv_subset

MetaDiv_subset creates analytical subsets from a MetaDiv Final Database.

## Modes

- `headers` → select biodiversity units using `query_lists/headers.txt`
- `sites` → select samples/sites using `query_lists/sites.txt`
- `taxonomy` → select biodiversity units by taxonomic rank

## Use

1. Edit the first notebook cell:
   - `MODE`
   - `SUBSET_MODE`
   - `SUBSET_NAME`
   - `FINAL_DB_NAME`
   - taxonomy settings if needed

2. Run the second cell.

The script automatically finds:

MetaDiv/output/<MODE>/FINAL_DB/<FINAL_DB_NAME>

## Output

Each subset contains:

final_db.csv  
summary.txt  

For_R/
- FOR_R_INFO.txt
- abundance.csv
- sample_metadata_to_fill.csv
- taxonomy.csv
- sequences.fasta
