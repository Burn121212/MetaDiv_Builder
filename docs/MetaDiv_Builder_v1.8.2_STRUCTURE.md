# MetaDiv Builder v1.8.2 — Directory Structure

This document reflects the **actual directory and file structure** included in the `MetaDiv_Builder-1.8.2.zip` release. It replaces the outdated v1.7.3 tree in `docs/MetaDiv Builder structure.txt`.

## Repository tree

```text
MetaDiv_Builder-1.8.2/
├── databases/
│   ├── ecological_reference_db/
│   │   ├── FAPROTAX.txt
│   │   └── Fungal_Traits_DB.txt
│   └── taxonomic_reference_db/
│       └── README.md
├── docs/
│   ├── Diagram_Structure_MetaDiv_Builder.png
│   ├── METADIV BUILDER LOGO V1.png
│   ├── MetaDiv Builder Platform.docx
│   └── MetaDiv Builder structure.txt
├── input/
│   ├── 16S/
│   │   ├── cyano_litter_abundance.csv
│   │   ├── cyano_litter_sequences.fasta
│   │   ├── cyano_litter_taxonomy.sintax.txt
│   │   ├── cyano_soil_abundance.csv
│   │   ├── cyano_soil_sequences.fasta
│   │   ├── cyano_soil_taxonomy.sintax.txt
│   │   └── readme.txt
│   ├── CO1/
│   │   ├── co1atlas_abundance.csv
│   │   ├── co1atlas_sequences.fasta
│   │   ├── co1atlas_taxonomy.sintax.txt
│   │   ├── co1tuxtlas_abundance.csv
│   │   ├── co1tuxtlas_sequences.fasta
│   │   ├── co1tuxtlas_taxonomy.sintax.txt
│   │   └── readme.txt
│   ├── ITS/
│   │   ├── ATLASMXA_abundance.csv
│   │   ├── ATLASMXA_sequences.fasta
│   │   ├── ATLASMXA_taxonomy.sintax.txt
│   │   ├── ATLASMXBC_abundance.csv
│   │   ├── ATLASMXBC_sequences.fasta
│   │   └── ATLASMXBC_taxonomy.sintax.txt
│   └── REAMDE.md
├── LICENSE
├── MetaDiv_FunEco/
│   ├── MetaDiv_FunEco_SCRIPT.ipynb
│   └── README.md
├── MetaDiv_R_scripts/
│   ├── 00_MetaDiv_create_phyloseq.R
│   ├── 01_MetaDiv_explorer_alpha_portable.R
│   ├── 02_MetaDiv_ordination_NMDS_Hellinger_taxglom_portable.R
│   ├── 03_MetaDiv_iNEXT_Hill_Numbers_portable.R
│   └── check_environment.R
├── MetaDiv_Utilies/
│   ├── MetaDiv_BLAST/
│   │   ├── query/
│   │   │   └── example_seq.fasta
│   │   ├── MetaDiv_BLAST.ipynb
│   │   └── README.md
│   ├── MetaDiv_KRONA/
│   │   ├── input_KRONA/
│   │   │   ├── metadata_gropus_example.txt
│   │   │   └── README.md
│   │   ├── MetaDiv_KRONA.ipynb
│   │   ├── MetaDiv_KRONA.py
│   │   └── README.md
│   ├── MetaDiv_metadata/
│   │   ├── source_metadata/
│   │   │   └── README.md
│   │   ├── MetaDiv_metadata.ipynb
│   │   └── README.md
│   ├── MetaDiv_sintax/
│   │   ├── MetaDiv_sintax.ipynb
│   │   ├── MetaDiv_sintax.py
│   │   └── README.md
│   ├── MetaDiv_subsets/
│   │   ├── query_lists/
│   │   │   ├── headers.txt
│   │   │   └── sites.txt
│   │   ├── MetaDiv_subset.ipynb
│   │   └── README.md
│   └── README.md
├── MetaDiv_ValStats/
│   ├── MetaDiv_ValStats_SCRIPT.ipynb
│   └── README.md
├── MetaDiv_Builder_v1_8_2_MAIN_SCRIPT.ipynb
├── MetaDiv_Builder_v1_8_2_MAIN_SCRIPT.py
├── metadiv_environment.yml
└── README.md
```

## Directory overview

| Path | Contents |
| --- | --- |
| `MetaDiv_Builder_v1_8_2_MAIN_SCRIPT.py` / `.ipynb` | Main Python workflow, in script and notebook formats. |
| `input/` | Example metabarcoding inputs organized by molecular marker (`ITS`, `16S`, `CO1`). |
| `databases/ecological_reference_db/` | Ecological reference tables (`FAPROTAX`, `Fungal_Traits_DB`). |
| `databases/taxonomic_reference_db/` | Taxonomic reference database location; the ZIP includes a README but no `.udb` files. |
| `MetaDiv_FunEco/` | Functional ecology notebook and documentation. |
| `MetaDiv_ValStats/` | Validation/statistics notebook and documentation. |
| `MetaDiv_R_scripts/` | R analysis scripts. |
| `MetaDiv_Utilies/` | BLAST, Krona, metadata, SINTAX and subset utilities (folder spelling preserved as supplied). |
| `docs/` | Documentation, images and the previous structure text file. |
| `metadiv_environment.yml` | Environment configuration. |

## Notes

- The tree records **files physically present in the supplied ZIP**, not hypothetical files generated during execution.
- An `output/` directory is **not present in this ZIP**. Output folders may be created at runtime; do not include them as if already committed.
- The supplied archive has an `input/REAMDE.md` filename and a `MetaDiv_Utilies/` directory name. These names have been retained exactly, rather than silently renamed.
- Update or replace `docs/MetaDiv Builder structure.txt` with this document, or link to it from the repository `README.md`.
