# MetaDiv Utilities

**MetaDiv Utilities** is a collection of companion tools designed to support
different steps of the MetaDiv Builder workflow.

These utilities provide additional functions for taxonomic assignment,
sequence comparison, metadata integration, interactive taxonomic visualization,
and database subsetting.

## Utilities

### MetaDiv_sintax

Performs taxonomic classification of representative sequences using the
**SINTAX** algorithm implemented in VSEARCH.

The utility can be used as a preprocessing step before MetaDiv Builder to
generate SINTAX-formatted taxonomic assignments from representative FASTA
sequences.

It supports different metabarcoding markers, including **ITS, 16S rRNA,
and COI**, using the corresponding reference databases.

---

### MetaDiv_BLAST

Provides tools for local sequence comparison using **BLAST**.

This utility allows sequences from MetaDiv databases to be compared against
local reference databases and can be used to explore or verify taxonomic
assignments based on sequence similarity.

---

### MetaDiv_KRONA

Generates **Krona-compatible taxonomic abundance tables** from MetaDiv
databases.

Taxonomic abundances can be summarized at different levels, including
phylum, class, order, family, genus, and species, allowing the results to be
visualized as interactive Krona plots.

The utility generates the input tables required by Krona; visualization can
then be performed using Krona-compatible tools.

---

### MetaDiv_metadata

Integrates external sample metadata with MetaDiv sample identifiers.

MetaDiv Builder generates basic sample metadata containing sample/site
identifiers. This utility matches those identifiers with an external metadata
table and appends additional ecological or methodological information, such as
**sequencing technology, country, biome, or other user-provided variables**.

The resulting metadata table can then be used in downstream ecological
analyses, including phyloseq workflows.

---

### MetaDiv_subsets

Creates subsets of MetaDiv biodiversity databases.

The utility allows users to extract selected portions of a database according
to criteria such as **taxonomic groups, sample lists, or user-defined
queries**, facilitating targeted ecological analyses without modifying the
original MetaDiv database.

---
