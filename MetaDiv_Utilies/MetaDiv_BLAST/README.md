# MetaDiv_BLAST
MetaDiv_BLAST compares query DNA sequences against a MetaDiv Final Database using NCBI BLAST.

## Requirements
- Docker Desktop
- Python with:
  - pandas
  - pathlib
    
## Install the NCBI BLAST Docker image with:
docker pull ncbi/blast  ncbi/blast:latest

# Input
 1. Select the MetaDiv mode
Edit the first cell of the notebook:
mode = "16S"
Examples:
mode = "16S"
mode = "ITS"
mode = "CO1"
 2. Select the Final Database
Set the exact Final Database file name:
final_database_name = "Final_Database_species_only_all_prokaryotes_p08_sppn08.csv"
The database must be located in:
MetaDiv/
└── output/
    └── <mode>/
        └── FINAL_DB/
            └── Final_Database_*.csv
The Final Database must contain the columns:
SPPN
sequence
SPPN is used as the BLAST reference sequence ID.
 3. Add the query FASTA
Place exactly one query FASTA file in:
MetaDiv_Utilities/
└── MetaDiv_BLAST/
    └── query/
        └── query.fasta
    
## Accepted extensions:
- .fasta
- .fa
- .fna

## Run
Run the notebook cells in order.
The script automatically:
1. Detects the MetaDiv project folder.
2. Loads the selected Final Database.
3. Creates a reference FASTA from the SPPN and sequence columns.
4. Builds the local BLAST database.
5. Runs blastn.
6. Selects the best hit for each query.
7. Adds the MetaDiv annotation to the best hits.

# Output
Results are saved in:
MetaDiv_Utilities/
└── MetaDiv_BLAST/
    └── blast_results/
Output files:
blast_results.tsv
blast_all_hits.csv
blast_best_hits.csv
blast_best_hits_annotated.csv
blast_best_hits_annotated.csv contains the best BLAST match for each query together with the corresponding information from the MetaDiv Final Database.
