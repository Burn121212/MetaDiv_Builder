# MetaDiv ValStats

**MetaDiv ValStats** is the validation module of the MetaDiv framework. It is designed to evaluate the internal consistency and integrity of the databases generated with **MetaDiv Builder**.

## Requirements

To run **MetaDiv ValStats**, your MetaDiv Builder analysis must first be executed with **Developer Mode enabled** in the **MetaDiv Main Script**.

Developer Mode generates the additional intermediate and validation files required by MetaDiv ValStats. These files are not produced during a standard MetaDiv Builder run.

Therefore:

> **MetaDiv ValStats can only be used with MetaDiv Builder outputs generated with Developer Mode enabled.**

## Workflow

1. Open the **MetaDiv Main Script**.
2. Enable **Developer Mode**.
3. Run the complete MetaDiv Builder workflow.
4. Once the analysis is finished, use the generated Developer Mode outputs as input for **MetaDiv ValStats**.
5. Run the desired validation procedures.

## Validation modules

MetaDiv ValStats currently includes validation procedures for:

- **Bootstrap threshold validation**
- **Representative sequence (RepSeq) validation**
- **SPPN validation**
- **Abundance validation**

These analyses are intended to verify that the taxonomic grouping, sequence selection, SPPN assignment, and abundance aggregation performed by MetaDiv Builder remain internally consistent.

## Important

Running MetaDiv Builder without Developer Mode will still generate the standard final database and downstream analysis files, but the additional information required by **MetaDiv ValStats** will not be available.

For validation analyses, always enable **Developer Mode before starting the MetaDiv Builder run**.
