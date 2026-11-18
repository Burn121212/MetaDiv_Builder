# MetaDiv FunEco

**MetaDiv FunEco** is the functional ecology module of the MetaDiv framework.

It generates **site-level fungal functional structural metrics** from the functional annotations produced by **MetaDiv Builder**, allowing the exploration of taxonomic–functional patterns across metabarcoding datasets.

The current version is designed for **ITS fungal datasets annotated with FungalTraits**.

---

## Requirements

Before running **MetaDiv FunEco**, you must first complete the corresponding **MetaDiv Builder ITS workflow** with functional annotation enabled.

MetaDiv FunEco uses the FungalTraits-annotated database generated inside:

`output/ITS/Functional_Ecology/`

The expected input is a file such as:

`Final_Database_ITS_FungalTraits_full_annotated_species_only.csv`

The input database must contain:

- `SPPN`
- FungalTraits primary lifestyle annotation
- FungalTraits secondary lifestyle annotation
- Sample abundance columns

MetaDiv FunEco automatically detects the sample abundance columns.

---

## User configuration

Before running the script, change the dataset path in:

`DATASET_DIR`

For example:

`DATASET_DIR = Path(r"C:/path/to/your/MetaDiv/output/ITS")`

The script will automatically locate the `Functional_Ecology` directory inside the selected ITS output folder.

---

## Functional classification

MetaDiv FunEco groups FungalTraits lifestyles into generalized trophic guilds:

- `PATOTROPH`
- `SIMBIOTROPH`
- `SAPROTROPH`
- `UNCLASSIFIED`
- `SUPERPATOTROPH`
- `SUPERSIMBIOTROPH`
- `SUPERSAPROTROPH`

Primary and secondary lifestyle annotations are used to characterize the functional structure of each sample.

---

## Functional structural metrics

MetaDiv FunEco currently evaluates four complementary strategies.

### 1. Functional Redundancy Index

`Functional_Redundancy_Index = Taxonomic_Richness / Lifestyle_Richness`

This metric estimates the average number of taxa associated with each detected ecological lifestyle.

Higher values indicate greater taxonomic redundancy among ecological functions.

### 2. Functional Compression Difference

`Functional_Compression_Difference = Lifestyle_Richness_Normalized - Generalized_Guild_Richness_Normalized`

This metric evaluates whether detailed lifestyle diversity is compressed into a smaller number of generalized trophic guilds.

By default, lifestyle richness is normalized using the maximum lifestyle richness observed in the current dataset.

An optional known total lifestyle universe can also be supplied.

### 3. Functional Compartmentalization Index

`Functional_Compartmentalization_Index = Dominant_Guild_Richness / Taxonomic_Richness`

This metric estimates how strongly the taxonomic richness of a site is concentrated within one dominant trophic guild.

The module also calculates a lifestyle-based alternative:

`Lifestyle_Compartmentalization_Index`

### 4. Functional Reorganization Score

The **Functional Reorganization Score** is an exploratory composite metric calculated from the scaled values of:

- Functional Redundancy Index
- Functional Compression Difference
- Functional Compartmentalization Index

The score is calculated as:

`Functional_Reorganization_Score = mean(Functional_Redundancy_Index_scaled, Functional_Compression_Difference_scaled, Functional_Compartmentalization_Index_scaled)`

This score summarizes multiple dimensions of functional community structure.

> **Important:** The Functional Reorganization Score is an exploratory metric and should not be interpreted as a direct conservation, disturbance, or ecosystem-quality index unless independently validated using environmental metadata.

---

## Additional metrics

MetaDiv FunEco also calculates:

- Taxonomic richness
- Lifestyle richness
- Generalized guild richness
- Guild richness
- Guild abundance
- Lifestyle counts by guild
- Shannon diversity of guild richness
- Guild richness evenness
- Simpson dominance
- Lifestyle distribution among guilds
- Pathotroph richness compartment
- Saprotroph richness compartment
- Simbiotroph richness compartment
- Relative abundance of major trophic guilds
- Functional Shift Indices

---

## Functional Shift Indices

Directional functional shifts are estimated by combining richness and relative abundance.

For example:

`Pathotroph_Functional_Shift_Index = Pathotroph_Richness_Compartment × Pathotroph_Relative_Abundance`

Equivalent indices are calculated for:

- Pathotrophs
- Saprotrophs
- Simbiotrophs

These metrics describe whether both taxonomic richness and abundance are jointly shifted toward a particular trophic strategy.

---

## LOW, MEDIUM and HIGH classes

Several functional metrics are classified as:

- `LOW`
- `MEDIUM`
- `HIGH`

using the 33% and 66% quantiles of the current dataset.

These classifications are therefore **dataset-relative** and should not be interpreted as universal ecological thresholds.

---

## Output

Results are generated inside:

`output/ITS/Functional_Ecology/Functional_Structural_Metrics_V4/`

The main output table is:

`Functional_Ecology_by_Site_V4_species.csv`

A log file is also generated:

`LOG_Functional_Structural_Metrics_V4_species.txt`

The log records:

- Input file
- Output file
- Run date
- Metric definitions
- Mathematical formulas
- Interpretation of the four functional strategies
- Lifestyle normalization strategy

---

## Workflow

MetaDiv Builder  
→ ITS Final Database  
→ FungalTraits annotation  
→ `Final_Database_ITS_FungalTraits_full_annotated.csv`  
→ MetaDiv FunEco  
→ Functional richness metrics  
→ Functional redundancy  
→ Functional compression  
→ Functional compartmentalization  
→ Functional shift indices  
→ Functional reorganization score  
→ `Functional_Ecology_by_Site_V4.csv`

---

## Python dependencies

MetaDiv FunEco requires:

- `pandas`
- `numpy`
- `pathlib`
- `datetime`

`pathlib` and `datetime` are included in the standard Python library.

---

## Important

MetaDiv FunEco is intended as a **downstream module of MetaDiv Builder**.

The module does not assign fungal functional traits directly from raw sequences. Functional information must first be incorporated into the MetaDiv database using the corresponding **Functional Ecology / FungalTraits annotation workflow**.

The resulting metrics are intended for exploratory ecological analysis and for evaluating relationships between **taxonomic diversity, functional diversity, trophic structure, and community reorganization**.
