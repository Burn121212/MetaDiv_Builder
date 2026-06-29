<p align="center">
  <img src="https://raw.githubusercontent.com/Burn121212/DB-Builder/main/docs/METADIV%20BUILDER%20LOGO%20V1.png" alt="MetaDiv Builder Logo" width="300"/>
</p>

# MetaDiv Builder

> **An integrated platform for metabarcoding data integration, ecological annotation, and biodiversity analysis**

MetaDiv Builder is an open-source platform designed to transform heterogeneous metabarcoding outputs into standardized biodiversity databases ready for ecological analyses. The platform integrates abundance tables, representative sequences, and SINTAX taxonomic classifications while preserving taxonomic traceability through persistent identifiers (SPPN).

Unlike traditional post-processing tools that focus on a single molecular marker or sequencing pipeline, MetaDiv Builder supports multiple metabarcoding markers (**ITS**, **16S rRNA**, and **CO1**) and provides a unified workflow for taxonomic harmonization, ecological annotation, database construction, visualization, and downstream analyses.

The platform currently integrates ecological reference databases such as **FungalTraits** (fungi) and **FAPROTAX** (prokaryotes), exports **phyloseq-compatible** datasets for R, generates **Krona** visualizations, and includes utilities for **BLAST** searches and taxonomic subset extraction.

---

# Table of Contents

- [Features](#features)
- [Workflow](#workflow)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Structure](#input-structure)
- [Directory Structure](#directory-structure)
- [Output Structure](#output-structure)
- [Taxonomic Collapse Strategies](#taxonomic-collapse-strategies)
- [Ecological Annotation](#ecological-annotation)
- [Integration with R](#integration-with-r)
- [License](#license)

---

# Features

✅ Supports multiple molecular markers

- ITS
- 16S 
- CO1

✅ Automatic dataset discovery

✅ Taxonomic harmonization

✅ Persistent taxonomic identifiers (SPPN)

✅ Taxonomic collapsing

- Species
- Genus
- Lowest confident rank

✅ Ecological annotation

- FungalTraits
- FAPROTAX

✅ Phyloseq-ready export

✅ Krona visualization

✅ Local BLAST utilities

✅ Taxonomic subset extraction

✅ Modular architecture for future expansion

---

# Workflow

MetaDiv Builder is organized into five independent but interconnected modules.

## 1. Dataset Integration

Automatically detects complete metabarcoding datasets and integrates:

- abundance tables
- representative sequences
- SINTAX taxonomy

Outputs:

- Concatenated datasets
- Fungi-only datasets
- Bacteria-only datasets
- Metazoa-only datasets
- All prokaryotes datasets
- All eukaryotes datasets
- Quality-control summaries

---

## 2. Taxonomic Harmonization

Performs:

- Taxonomy parsing
- Confidence filtering
- Taxonomic collapsing
- Representative sequence selection
- SPPN identifier assignment

Supported collapse strategies:

- `species_only`
- `genus`
- `lowest_rank`

---

## 3. Ecological Annotation

Adds ecological information using specialized databases.

Current supported databases:

- **FungalTraits** (ITS)
- **FAPROTAX** (16S)

The ecological annotation module is independent from the taxonomic engine, allowing future ecological databases to be incorporated without modifying the core workflow.

---

## 4. Output Generation

Creates standardized datasets for downstream analyses.

Outputs include:

- Final integrated database
- Functional ecology tables
- Phyloseq-ready files
- Krona tables
- Representative FASTA files

---

## 5. Accessory Modules

Additional utilities include:

- Local BLAST searches
- Taxonomic subset extraction
- Future analytical extensions

---

# Installation

## Windows

### 1. Install Miniconda

Download and install Miniconda:

https://www.anaconda.com/download

---

### 2. Clone the repository

```bash
git clone https://github.com/Burn121212/MetaDiv-Builder.git

cd MetaDiv-Builder
```

---

### 3. Create the MetaDiv Builder environment

```bash
conda env create -f metadiv_builder_win.yml
```

---

### 4. Activate the environment

```bash
conda activate metadiv
```

> **Note:** If your environment uses a different name, replace `metadiv` with the name specified inside `metadiv_builder_win.yml`.

---

### 5. Launch Jupyter Notebook

```bash
jupyter notebook
```

Open

```
MetaDiv_Builder_v1_7.ipynb
```

and run all notebook cells.

---

# macOS / Linux

## 1. Install Miniconda

Download and install Miniconda:

https://www.anaconda.com/download

---

## 2. Clone the repository

```bash
git clone https://github.com/Burn121212/MetaDiv-Builder.git

cd MetaDiv-Builder
```

---

## 3. Create the MetaDiv Builder environment

```bash
conda env create -f metadiv_builder_linux_ios.yml
```

---

## 4. Activate the environment

```bash
conda activate metadiv
```

> **Note:** If your environment uses a different name, replace `metadiv` with the name specified inside `metadiv_builder_linux_ios.yml`.

---

## 5. Launch Jupyter Notebook

```bash
jupyter notebook
```

Open

```
MetaDiv_Builder_v1_7.ipynb
```

and run all notebook cells.

---

Each dataset must contain:

```
DatasetName_abundance.csv

DatasetName_taxonomy.sintax

DatasetName_sequences.fasta
```

Configure the parameters in the first section of the notebook and execute all cells.

MetaDiv Builder automatically detects complete datasets and generates the corresponding outputs.

---

# Input Structure

Each molecular marker has its own folder.

```
input/

├── ITS/
│   ├── Dataset_abundance.csv
│   ├── Dataset_taxonomy.sintax
│   └── Dataset_sequences.fasta
│
├── 16S/
│   ├── Dataset_abundance.csv
│   ├── Dataset_taxonomy.sintax
│   └── Dataset_sequences.fasta
│
└── CO1/
    ├── Dataset_abundance.csv
    ├── Dataset_taxonomy.sintax
    └── Dataset_sequences.fasta
```

---

# Directory Structure

```
MetaDiv Builder/

├── input/
├── databases/
├── output/
├── R/
├── BLAST/
├── subsets/
└── docs/
```

---

# Output Structure

Independent output folders are generated for each molecular marker.

```
output/

├── ITS/
├── 16S/
└── CO1/
```

Each output contains:

- Final database
- Functional ecology
- For_R
- Krona
- Concatenated tables
- Log files

---

# Taxonomic Collapse Strategies

## `species_only`

Collapse only when the species assignment satisfies the selected confidence threshold.

---

## `genus`

Collapse at the genus level.

---

## `lowest_rank`

Automatically collapses at the deepest taxonomic rank satisfying the selected confidence threshold.

---

# Ecological Annotation

| Marker | Ecological Database |
|----------|--------------------|
| ITS | FungalTraits |
| 16S | FAPROTAX |
| CO1 | Planned |

---

# Integration with R

MetaDiv Builder exports standardized files directly compatible with the official R analysis scripts.

Generated files include:

- abundance table
- taxonomy table
- representative sequences
- ecological annotation tables
- sample metadata template

These files can be imported directly into:

- phyloseq
- vegan
- hillR
- picante
- ggplot2

allowing fully reproducible biodiversity analyses.

---

# Citation

If you use MetaDiv Builder in your research, please cite:

**Águila B., Romero-Guiterrez M. F. et al. **

*MetaDiv Builder: An integrated platform for metabarcoding data integration, ecological annotation, and biodiversity analysis.*

(Manuscript in preparation.)

---

# License

**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**

---

## Authors

**Bernardo Águila, Miguel F. Romero-Guiterrez**  
Instituto de Biología  
Universidad Nacional Autónoma de México (UNAM)

Contributions and suggestions are welcome.
