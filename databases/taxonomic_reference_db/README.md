# Taxonomic reference databases for MetaDiv

Best practices for MetaDiv_Builer pre-processing include annotating all fasta files using the same refference database. We provide tested versions of the following databases through Google Drive:

| Database | Marker Gene | Google Drive code | Reference |
| -------------| ------ | ------------- | ----------|
|EUKARYOME ITS v2.0|ITS|1TU7WiKNN6cHG8_dn0GpDT6bb-YbBisF5|1|
|ITGDB|16S rRNA|14RO-uq4KJwDtaVzU9XfFlJnPzOAx3Bgc|2|
|MIDORI2|COI|1m6SZWsSYDR5zYh3JO6OK-ErBW8ifFBfH|3|


## Setting up databases in the command line

It is possible to setup the databases through the command line using `gdown` following the steps below.

Starting in the MetaDiv_Builder directory:

1. Activate the MetaDiv conda environment:

```bash
conda activate metadiv
```

2. Install the `gdown` tool if  not available already:

```bash
conda install -c conda-forge gdown -y -q
```

3. Download a reference file using the Google Drive code:

```bash
gdown --id 1TU7WiKNN6cHG8_dn0GpDT6bb-YbBisF5 \
-O databases/taxonomic_reference_db/SINTAX_EUKARYOME_ITS_v2.0.fasta
```
## References

1. Tedersoo, L., Hosseyni Moghaddam, M. S., Mikryukov, V., Hakimzadeh, A., Bahram, M., Nilsson, R. H., ... & Anslan, S. (2024). EUKARYOME: the rRNA gene reference database for identification of all eukaryotes. *Database*, 2024, baae043.
2. Hsieh, Y. P., Hung, Y. M., Tsai, M. H., Lai, L. C., & Chuang, E. Y. (2022). 16S-ITGDB: an integrated database for improving species classification of prokaryotic 16S ribosomal RNA sequences. *Frontiers in Bioinformatics*, 2, 905489.
3. Leray, M., Knowlton, N., & Machida, R. J. (2022). MIDORI2: A collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. *Environmental Dna*, 4(4), 894-907.
