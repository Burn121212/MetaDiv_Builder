# MetaDiv sintax

VSEARCH/SINTAX taxonomic classification utility for MetaDiv Builder

## CLI script for Linux/MacOS

Alternative to the jupyter notebook. It is possible to install vsearch in a dedicated anaconda environment.

### Quickstart

Starting in the MetaDiv base directory:

1. Activate metadiv environment

```bash
conda activate metadiv
```

2. Display options and examples

```bash
python MetaDiv_Utilies/MetaDiv_sintax/MetaDiv_sintax.py --help
```

3. Run example assuming database is availabe (see instructions in MetaDiv_Builder/databases/README.md)

```bash
python MetaDiv_Utilies/MetaDiv_sintax/MetaDiv_sintax.py \
--mode ITS --all-fastas --reference-db SINTAX_EUKARYOME_ITS_v2.0.fasta \
 --runtime conda --conda-env-name MetaDiv_sintax --auto-create-conda-env
```

### Options

```bash
  -h, --help            show this help message and exit
  -m MODE, --mode MODE  Genetic marker. Options: ITS, 16S, COI. CO1 is accepted as an alias for COI. (default: None)
  -f FASTA_FILES [FASTA_FILES ...], --fasta-files FASTA_FILES [FASTA_FILES ...]
                        One or more FASTA filenames to classify. Paths are relative to the selected input directory. Subfolders are accepted.
                        (default: None)
  --all-fastas          Classify every supported FASTA file found recursively in the selected input directory. (default: False)
  -d REFERENCE_DB, --reference-db REFERENCE_DB
                        Reference SINTAX database filename, usually stored inside databases/taxonomic_reference_db/. (default: None)
  --project-root PROJECT_ROOT
                        MetaDiv Builder project root. If omitted, the script searches upward from the current directory. (default: None)
  --input-dir INPUT_DIR
                        Custom input directory. If omitted, input/<MODE>/ is used. (default: None)
  --reference-db-dir REFERENCE_DB_DIR
                        Custom taxonomic reference database directory. (default: None)
  --sintax-cutoff SINTAX_CUTOFF
                        VSEARCH --sintax_cutoff value. (default: 0.8)
  --strand {plus,both}  VSEARCH --strand option. (default: both)
  -t THREADS, --threads THREADS
                        Number of VSEARCH threads. (default: 16)
  --overwrite-existing  Overwrite existing <sample>_sintax.txt files. (default: True)
  --skip-existing       Skip FASTA files that already have a <sample>_sintax.txt output. (default: True)
  --runtime {docker,conda}
                        Runtime used to execute VSEARCH. (default: docker)
  --docker-image DOCKER_IMAGE
                        Docker image used when --runtime docker. (default: pipecraft/vsearch:2.30.4-pc1.2.0)
  --auto-pull-image     Automatically pull the Docker image if it is not available locally. (default: True)
  --no-auto-pull-image  Do not pull a missing Docker image automatically. (default: True)
  --conda-env-name CONDA_ENV_NAME
                        Conda environment used when --runtime conda. (default: MetaDiv_sintax)
  --auto-create-conda-env
                        When --runtime conda, create or repair the Conda environment with VSEARCH if it is missing or VSEARCH cannot be
                        executed. (default: False)
  --vsearch-package VSEARCH_PACKAGE
                        Conda package spec used with --auto-create-conda-env. Examples: vsearch, vsearch=2.30.0. (default: vsearch)
  --report-file REPORT_FILE
                        Processing report path. Default: ./MetaDiv_sintax_report.txt. (default: None)
```
