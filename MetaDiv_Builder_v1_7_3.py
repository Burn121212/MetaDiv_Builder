#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
MetaDiv Builder v1.7.3 - Heterogeneous Metabarcoding Integrator
from OTUs/ASVs/zOTUs/MOTUs to Diversity

Authors: Bernardo Águila, UNAM
         Miguel F. Romero, UNAM

A unified pipeline for integrating metabarcoding data from multiple sources.
General-purpose post-processing framework for large-scale metabarcoding meta-analyses.

================================================================================
INPUT FILE REQUIREMENTS
================================================================================

All raw input files must be placed inside a marker-specific input folder:

    input/ITS/
    input/16S/
    input/CO1/

The active folder is selected automatically from MODE.

Each dataset must contain three files with matching dataset prefixes:

    DatasetName_abundance.csv
    DatasetName_taxonomy.sintax
    DatasetName_sequences.fasta

File descriptions:

    *_abundance.csv
        Abundance table.
        The first column must contain ASV/OTU/zOTU/MOTU IDs.
        Remaining columns must correspond to samples.
        The first column name is flexible.

    *_sequences.fasta
        Representative sequence file.
        FASTA headers must match the IDs used in the abundance table
        and taxonomy file.

    *_taxonomy.sintax
        SINTAX taxonomy annotation file (generated from sequences.fasta).
        Must contain taxonomic ranks and confidence values when available.
        
================================================================================
COLLAPSE STRATEGIES
================================================================================

COLLAPSE_STRATEGY defines how ASVs/OTUs/zOTUs/MOTUs are grouped into
ecological diversity units.

Available options:

    "species_only"
        Species-level collapse only.
        Records are collapsed only when they share the same species-level
        annotation and pass the selected confidence threshold.

    "genus"
        Genus-level collapse only.
        Records are collapsed when they share the same genus-level annotation
        and pass the selected confidence threshold.

    "all"
        Lowest-rank recursive collapse.
        Records are recursively collapsed at the lowest available shared
        taxonomic rank that passes the selected confidence threshold.
        This is the most aggressive collapse strategy.
================================================================================
Key features:
- Always creates intermediate log files and concatenated tables for reliable collapse
- Removes temporary folder after successful completion (unless DEV_MODE=True)
- Supports ITS (fungi), 16S (prokaryotes), and CO1 (eukaryotes/metazoa) markers
- For ITS: integrates Fungal Traits Database (primary_lifestyle, Secondary_lifestyle)
- For 16S: directly integrates FAPROTAX.txt to generate multi-function ecological profiles
- For CO1: creates a Functional_Ecology placeholder for future ecological/trait annotation
- Exports Krona-compatible tables with scalable aggregation/row limits for Excel-template workflows
- First column of abundance table = ID (no column name restrictions)


"""

import pandas as pd
import shutil
from pathlib import Path
from collections import defaultdict
import time
import re
from datetime import datetime

# ============================================================================
# CONFIGURATION - USER MODIFIES THESE SETTINGS
# ============================================================================

MODE = "16S"  # "16S", "ITS", or "CO1"
SUBSET_MODE = "all_prokaryotes"  # ITS: "only_fungi"/"all_eukaryotes"; 16S: "only_bacteria"/"all_prokaryotes"; CO1: "only_metazoa"/"all_eukaryotes"
DEV_MODE = False  # If True, keeps intermediate concatenated tables; if False, deletes them
COLLAPSE_STRATEGY = "species_only" # species_only, genus, all
P_VALUE_THRESHOLD = 1.00 # Taxonomic collapse confidence threshold (higher = stricter, lower = more permissive). range from 1.00 to 0.00
SPPN_P_THRESHOLD = 0.80  # SPPN assignment confidence threshold (recommended = 0.80)
INPUT_ROOT = Path("input")
# Marker-specific input folder.
# Examples:
#   MODE = "ITS" -> input/ITS/
#   MODE = "16S" -> input/16S/
#   MODE = "CO1" -> input/CO1/
INPUT_DIR = INPUT_ROOT / MODE
DATABASE_DIR = Path("databases")
ECOLOGICAL_REFERENCE_DB_DIR = DATABASE_DIR / "ecological_reference_db"
TAXONOMIC_REFERENCE_DB_DIR = DATABASE_DIR / "taxonomic_reference_db"
OUTPUT_ROOT = Path("output")

# Optional ITS FungalTraits Functional Ecology module.
# Expected input file: databases/ecological_reference_db/Fungal_Traits_DB.txt
# When True, it exports a complete ITS annotation table in:
# output/ITS/Functional_Ecology/
# This keeps the regular FINAL_DB compact, but creates an additional version
# with all columns available in the FungalTraits database.
RUN_FUNGALTRAITS_ITS = False
FUNGAL_TRAITS_DB = ECOLOGICAL_REFERENCE_DB_DIR / "Fungal_Traits_DB.txt"

# Optional 16S FAPROTAX Functional Ecology module.
# Expected input file: databases/ecological_reference_db/FAPROTAX.txt
# This module reads the FAPROTAX database directly and does not require running
# collapse_table.py separately. It preserves the MetaDiv taxonomic table and
# generates additional functional summaries in output/16S/Functional_Ecology/.
RUN_FAPROTAX_16S = True
FAPROTAX_DB_16S = ECOLOGICAL_REFERENCE_DB_DIR / "FAPROTAX.txt"

# Optional Krona visualization module.
# This module writes Krona-compatible text tables but does not render HTML from Python.
# It is designed to work with either KronaTools or an external Excel/Krona template.
#
# KRONA_MODE options:
#   "FULL"    = use the full available hierarchy down to species.
#   "TOP"     = use the full hierarchy but keep only the most abundant rows.
#   "PHYLUM"  = aggregate paths to phylum.
#   "CLASS"   = aggregate paths to class.
#   "ORDER"   = aggregate paths to order.
#   "FAMILY"  = aggregate paths to family.
#   "GENUS"   = aggregate paths to genus. Recommended default for huge meta-analyses.
#   "SPECIES" = aggregate paths to species.
#
# If the resulting table exceeds KRONA_MAX_ROWS, MetaDiv Builder applies
# KRONA_LIMIT_ACTION. The safest option for Excel-template workflows is "TOP".
RUN_KRONA_EXPORT = True
RUN_KRONA_PER_SAMPLE_EXPORT = False
KRONA_MODE = "GENUS"
KRONA_MIN_ABUNDANCE = 1
KRONA_MAX_ROWS = 30000
KRONA_TOP_ROWS = 30000
KRONA_LIMIT_ACTION = "TOP"  # "TOP" or "WARNING_ONLY"
KRONA_FILL_UNCLASSIFIED = True
KRONA_EXPORT_FIXED_COLUMNS = True

# ============================================================================
# DERIVED CONFIGURATION
# ============================================================================

OUTPUT_DIR = OUTPUT_ROOT / MODE
FINAL_DB_DIR = OUTPUT_DIR / "FINAL_DB"
FOR_R_DIR = OUTPUT_DIR / "For_R"
FUNCTIONAL_ECOLOGY_DIR = OUTPUT_DIR / "Functional_Ecology"
KRONA_DIR = OUTPUT_DIR / "Krona"
CONCAT_DIR = OUTPUT_DIR / "concatenated_tables"  # Always created, deleted later if not DEV_MODE

# Create output folders for the active marker only.
# This prevents runs in one marker mode from modifying output folders from
# another marker mode.
for d in [OUTPUT_ROOT, OUTPUT_DIR, FINAL_DB_DIR, FOR_R_DIR, CONCAT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Define all possible subdirectories (initialize as None, set later)
FUNGI_DIR = ALL_EUK_DIR = NON_ANNOTATED_DIR = None
BACTERIA_DIR = ALL_PROK_DIR = None
CO1_DIR = METAZOA_DIR = None

# Subdirectories for filtered tables
if MODE == "ITS":
    FUNGI_DIR = CONCAT_DIR / "Fungi_concatenated"
    ALL_EUK_DIR = CONCAT_DIR / "AllEuk_concatenated"
    NON_ANNOTATED_DIR = CONCAT_DIR / "Non_annotated"
    if SUBSET_MODE == "only_fungi":
        TARGET_DIR = FUNGI_DIR
        PATTERN = "*_fungi.csv"
        SUBSET_NAME = "only_fungi"
        FILTER_DESCRIPTION = "only fungi"
    elif SUBSET_MODE == "all_eukaryotes":
        TARGET_DIR = ALL_EUK_DIR
        PATTERN = "*_all_eukaryotes.csv"
        SUBSET_NAME = "all_eukaryotes"
        FILTER_DESCRIPTION = "all eukaryotes (no filter)"
    else:
        raise ValueError("For MODE='ITS', SUBSET_MODE must be 'only_fungi' or 'all_eukaryotes'.")
    MARKER_NAME = "ITS"
    FUNGAL_TRAITS_DB = ECOLOGICAL_REFERENCE_DB_DIR / "Fungal_Traits_DB.txt"
    for d in [FUNGI_DIR, ALL_EUK_DIR, NON_ANNOTATED_DIR]:
        d.mkdir(parents=True, exist_ok=True)

elif MODE == "16S":
    BACTERIA_DIR = CONCAT_DIR / "OnlyBacteria_concatenated"
    ALL_PROK_DIR = CONCAT_DIR / "AllProkaryotes_concatenated"
    NON_ANNOTATED_DIR = CONCAT_DIR / "Non_annotated"
    if SUBSET_MODE == "only_bacteria":
        TARGET_DIR = BACTERIA_DIR
        PATTERN = "*_only_bacteria.csv"
        SUBSET_NAME = "only_bacteria"
        FILTER_DESCRIPTION = "only bacteria"
    elif SUBSET_MODE == "all_prokaryotes":
        TARGET_DIR = ALL_PROK_DIR
        PATTERN = "*_all_prokaryotes.csv"
        SUBSET_NAME = "all_prokaryotes"
        FILTER_DESCRIPTION = "all prokaryotes (bacteria + archaea)"
    else:
        raise ValueError("For MODE='16S', SUBSET_MODE must be 'only_bacteria' or 'all_prokaryotes'.")
    MARKER_NAME = "16S"
    FUNGAL_TRAITS_DB = None
    for d in [BACTERIA_DIR, ALL_PROK_DIR, NON_ANNOTATED_DIR]:
        d.mkdir(parents=True, exist_ok=True)

elif MODE == "CO1":
    CO1_DIR = CONCAT_DIR / "AllCO1_concatenated"
    METAZOA_DIR = CONCAT_DIR / "Metazoa_concatenated"
    NON_ANNOTATED_DIR = CONCAT_DIR / "Non_annotated"
    if SUBSET_MODE == "only_metazoa":
        TARGET_DIR = METAZOA_DIR
        PATTERN = "*_only_metazoa.csv"
        SUBSET_NAME = "only_metazoa"
        FILTER_DESCRIPTION = "only metazoa"
    elif SUBSET_MODE == "all_eukaryotes":
        TARGET_DIR = CO1_DIR
        PATTERN = "*_all_eukaryotes.csv"
        SUBSET_NAME = "all_eukaryotes"
        FILTER_DESCRIPTION = "all CO1 eukaryotes"
    else:
        raise ValueError("For MODE='CO1', SUBSET_MODE must be 'only_metazoa' or 'all_eukaryotes'.")
    MARKER_NAME = "CO1"
    FUNGAL_TRAITS_DB = None
    for d in [CO1_DIR, METAZOA_DIR, NON_ANNOTATED_DIR]:
        d.mkdir(parents=True, exist_ok=True)

else:
    raise ValueError("MODE must be 'ITS', '16S', or 'CO1'.")

TAX_RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PVALUE_COLS = [f"{r}_pvalue" for r in TAX_RANKS]

# Broad animal phyla used for CO1 "only_metazoa" filtering.
# MIDORI2/NCBI-style SINTAX annotations often place animal lineages directly
# at the phylum rank, for example p:Arthropoda_6656, rather than using
# p:Metazoa. Therefore, CO1 metazoan filtering cannot rely only on the word
# "Metazoa".
CO1_METAZOA_PHYLA = [
    "acanthocephala",
    "annelida",
    "arthropoda",
    "brachiopoda",
    "bryozoa",
    "chaetognatha",
    "chordata",
    "cnidaria",
    "ctenophora",
    "echinodermata",
    "entoprocta",
    "gastrotricha",
    "gnathostomulida",
    "hemichordata",
    "kinorhyncha",
    "loricifera",
    "micrognathozoa",
    "mollusca",
    "nematoda",
    "nematomorpha",
    "nemertea",
    "onychophora",
    "orthonectida",
    "phoronida",
    "placozoa",
    "platyhelminthes",
    "porifera",
    "priapulida",
    "rhombozoa",
    "rotifera",
    "sipuncula",
    "tardigrada",
    "xenacoelomorpha",
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def format_threshold_short(value: float) -> str:
    """
    Format confidence thresholds for compact and readable folder/file names.

    Examples
    --------
    1.0  -> "1"
    0.9  -> "09"
    0.8  -> "08"
    0.05 -> "005"
    """
    value = float(value)
    if value.is_integer():
        return str(int(value))
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    if text.startswith("0."):
        return "0" + text.split(".", 1)[1]
    return text.replace(".", "")

def format_threshold_for_name(value: float) -> str:
    """Return p-threshold tag, e.g. 1.0 -> p1 and 0.8 -> p08."""
    return f"p{format_threshold_short(value)}"

def build_final_db_filename(suffix: str) -> str:
    return f"Final_Database_{COLLAPSE_STRATEGY}_{suffix}_{format_threshold_for_name(P_VALUE_THRESHOLD)}.csv"

def build_collapse_log_filename(suffix: str) -> str:
    return f"Collapse_Log_{COLLAPSE_STRATEGY}_{suffix}_{format_threshold_for_name(P_VALUE_THRESHOLD)}.txt"

def build_for_r_tag(suffix: str) -> str:
    collapse_tag = format_threshold_for_name(P_VALUE_THRESHOLD)
    sppn_tag = f"sppn{format_threshold_short(SPPN_P_THRESHOLD)}"
    return f"{COLLAPSE_STRATEGY}_{collapse_tag}_{sppn_tag}"

def build_for_r_output_dir(suffix: str) -> Path:
    """
    Build a clean nested For_R folder structure.

    Example
    -------
    For_R/all_prokaryotes/species_only_p1_sppn08/
    """
    return FOR_R_DIR / suffix / build_for_r_tag(suffix)

def write_report(report_lines: list, report_path: Path):
    report_text = "\n".join(report_lines)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_text)
    print(report_text)


def validate_reference_databases():
    """
    Validate optional reference databases stored in the databases/ folder.

    input/<MODE>/ contains user datasets.
    databases/ecological_reference_db/ contains resources used directly by
    optional functional ecology modules.
    databases/taxonomic_reference_db/ contains preprocessing taxonomic resources
    such as Eukaryome, SILVA, ITGDB or MIDORI2. These are not used directly by
    this script because taxonomy.sintax is already expected as input.
    """
    print("\n📁 CHECKING MARKER-SPECIFIC INPUT FOLDER...")
    print(f"   Active MODE: {MODE}")
    print(f"   Expected input folder: {INPUT_DIR}")
    if INPUT_DIR.exists():
        print(f"   ✅ Input folder found: {INPUT_DIR}")
    else:
        print(f"   ⚠️ Input folder not found: {INPUT_DIR}")
        print("   ⚠️ Create this folder and place the active marker datasets inside it.")

    print("\n📚 CHECKING REFERENCE DATABASES...")
    if MODE == "ITS":
        expected = ECOLOGICAL_REFERENCE_DB_DIR / "Fungal_Traits_DB.txt"
        if expected.exists():
            print(f"   ✅ FungalTraits database found: {expected}")
        else:
            print(f"   ⚠️ FungalTraits database not found: {expected}")
            print("   ⚠️ ITS lifestyle annotation will be skipped or left empty.")
    if MODE == "16S" and RUN_FAPROTAX_16S:
        expected = ECOLOGICAL_REFERENCE_DB_DIR / "FAPROTAX.txt"
        if expected.exists():
            print(f"   ✅ FAPROTAX database found: {expected}")
        else:
            print(f"   ⚠️ FAPROTAX database not found: {expected}")
            print("   ⚠️ 16S functional annotation will be skipped.")


def build_co1_metazoa_mask(tax_lower: pd.Series) -> pd.Series:
    """
    Build a conservative CO1 Metazoa mask from SINTAX taxonomy strings.

    MIDORI2/NCBI-style CO1 annotations may not explicitly include "Metazoa" or
    "Animalia" as a rank. For example, animal records can appear as:
        k:Eukaryota, p:Arthropoda, c:Insecta, ...

    Therefore, this mask combines:
    - explicit Metazoa/Animalia labels when present
    - a curated list of major animal phyla at the phylum rank
    """
    explicit_metazoa = (
        tax_lower.str.contains("p:metazoa", na=False, case=False) |
        tax_lower.str.contains("p:animalia", na=False, case=False) |
        tax_lower.str.contains("k:metazoa", na=False, case=False) |
        tax_lower.str.contains("k:animalia", na=False, case=False)
    )

    phylum_pattern = "|".join([f"p:{re.escape(p)}" for p in CO1_METAZOA_PHYLA])
    animal_phyla = tax_lower.str.contains(
        phylum_pattern,
        na=False,
        case=False,
        regex=True
    )

    return explicit_metazoa | animal_phyla


# ============================================================================
# FUNGAL TRAITS INTEGRATION (ITS only)
# ============================================================================

def load_fungal_traits_db(traits_path: Path) -> pd.DataFrame:
    """
    Load the FungalTraits database from a tab-separated file.

    The regular FINAL_DB uses only:
    - primary_lifestyle
    - Secondary_lifestyle

    The ITS Functional_Ecology export uses all available FungalTraits columns.
    A normalized matching key named genus_clean is added internally.
    """
    print("\n📁 LOADING FUNGAL TRAITS DATABASE...")
    if not traits_path.exists():
        print(f"   ⚠️ Warning: {traits_path} not found. Lifestyle columns will be empty.")
        return pd.DataFrame()

    try:
        traits = pd.read_csv(traits_path, sep="\t", low_memory=False)
        traits.columns = [str(c).strip() for c in traits.columns]

        genus_col = None
        for col in traits.columns:
            if col.lower() == "genus":
                genus_col = col
                break
        if genus_col is None:
            for col in traits.columns:
                if "genus" in col.lower():
                    genus_col = col
                    break
        if genus_col is None:
            print("   ❌ ERROR: No genus column found in FungalTraits DB.")
            return pd.DataFrame()

        primary_col = None
        secondary_col = None
        for col in traits.columns:
            col_norm = col.lower().replace(" ", "_")
            if "primary_lifestyle" in col_norm:
                primary_col = col
            if "secondary_lifestyle" in col_norm:
                secondary_col = col

        if primary_col is None:
            print("   ⚠️ Warning: primary_lifestyle column not found. It will be created empty.")
            primary_col = "primary_lifestyle"
            traits[primary_col] = ""
        if secondary_col is None:
            print("   ⚠️ Warning: secondary_lifestyle column not found. It will be created empty.")
            secondary_col = "Secondary_lifestyle"
            traits[secondary_col] = ""

        rename_map = {}
        if genus_col != "genus":
            rename_map[genus_col] = "genus"
        if primary_col != "primary_lifestyle":
            rename_map[primary_col] = "primary_lifestyle"
        if secondary_col != "Secondary_lifestyle":
            rename_map[secondary_col] = "Secondary_lifestyle"
        traits = traits.rename(columns=rename_map)

        traits["genus"] = traits["genus"].astype(str).str.strip()
        traits["genus_clean"] = traits["genus"].str.lower()
        traits = traits[traits["genus_clean"].notna() & traits["genus_clean"].ne("")].copy()
        traits = traits.drop_duplicates(subset=["genus_clean"], keep="first").reset_index(drop=True)

        print(f"   ✅ Loaded {len(traits):,} unique genera with FungalTraits information.")
        print(f"   ✅ FungalTraits columns available for Functional_Ecology: {len(traits.columns) - 1:,}")
        return traits

    except Exception as e:
        print(f"   ❌ Error loading fungal traits database: {e}")
        return pd.DataFrame()


def prepare_functional_ecology_export_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare Functional_Ecology tables for user-facing ecological annotation.

    FINAL_DB keeps OTU_XX and Original_ID for traceability, but Functional_Ecology
    exports use SPPN as the stable feature identifier and hide intermediate IDs.

    Important:
    - Keep parsed taxonomy columns (domain, phylum, class, order, family, genus, species)
      because these are already filtered/parsed consistently by MetaDiv Builder.
    - Remove p-value columns from this ecological export to keep it clean.
      The p-values remain available in FINAL_DB for traceability.
    - Remove the raw sintax_taxonomy string from this export because it can contain
      low-confidence names below the selected threshold and should not be used as the
      final ecological taxonomy table.
    """
    out = df.copy()

    # If parsed taxonomy columns are missing for any reason, recover them from the
    # raw SINTAX string before dropping sintax_taxonomy.
    if "sintax_taxonomy" in out.columns:
        parsed_tax = out["sintax_taxonomy"].fillna("").apply(parse_tax_columns_from_sintax).apply(pd.Series)
        for col in parsed_tax.columns:
            if col not in out.columns:
                out[col] = parsed_tax[col]

    functional_drop_cols = ["OTU_XX", "Original_ID", "sintax_taxonomy"] + [f"{r}_pvalue" for r in TAX_RANKS]
    cols_to_drop = [c for c in functional_drop_cols if c in out.columns]
    if cols_to_drop:
        out = out.drop(columns=cols_to_drop)

    # Put stable ID and parsed taxonomy at the front; keep sequence, abundances,
    # FungalTraits, and other annotation columns after them.
    front_cols = [c for c in ["SPPN"] + TAX_RANKS if c in out.columns]
    remaining_cols = [c for c in out.columns if c not in front_cols]
    out = out[front_cols + remaining_cols]

    return out


def add_fungaltraits_functional_ecology_its(final_df: pd.DataFrame, output_dir: Path, traits_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create the ITS Functional_Ecology export using all FungalTraits columns.

    Output:
        output/ITS/Functional_Ecology/Final_Database_ITS_FungalTraits_full_annotated.csv
        output/ITS/Functional_Ecology/FungalTraits_Module_Log.txt

    The regular FINAL_DB remains compact and keeps only primary_lifestyle and
    Secondary_lifestyle. This function writes a separate complete annotation table.
    """
    if MODE != "ITS" or not RUN_FUNGALTRAITS_ITS:
        return final_df

    functional_dir = output_dir / "Functional_Ecology"
    functional_dir.mkdir(parents=True, exist_ok=True)

    if traits_df is None or traits_df.empty or "genus" not in final_df.columns:
        annotated = prepare_functional_ecology_export_table(final_df)
        annotated_path = functional_dir / "Final_Database_ITS_FungalTraits_full_annotated.csv"
        annotated.to_csv(annotated_path, index=False, encoding="utf-8")
        log_text = (
            "FungalTraits ITS Functional Ecology module\n"
            "Status: no complete annotation was added.\n"
            "Reason: FungalTraits database is missing/empty or final database has no genus column.\n"
        )
        (functional_dir / "FungalTraits_Module_Log.txt").write_text(log_text, encoding="utf-8")
        print(f"   ⚠️ ITS Functional_Ecology created without extra FungalTraits columns: {annotated_path}")
        return final_df

    annotated = final_df.copy()
    annotated["genus_clean"] = annotated["genus"].astype(str).str.strip().str.lower()

    ft = traits_df.copy()
    rename_ft = {}
    for col in ft.columns:
        if col == "genus_clean":
            continue
        if col in annotated.columns:
            rename_ft[col] = f"FungalTraits_{col}"
    ft = ft.rename(columns=rename_ft)

    annotated = annotated.merge(ft, on="genus_clean", how="left")
    annotated = annotated.drop(columns=["genus_clean"])

    ft_export_cols = [c for c in ft.columns if c != "genus_clean"]
    matched_col = "FungalTraits_genus" if "FungalTraits_genus" in annotated.columns else "genus"
    matched = int(annotated[matched_col].notna().sum()) if matched_col in annotated.columns else 0
    total = len(annotated)

    annotated = prepare_functional_ecology_export_table(annotated)

    annotated_path = functional_dir / "Final_Database_ITS_FungalTraits_full_annotated.csv"
    annotated.to_csv(annotated_path, index=False, encoding="utf-8")

    log_text = f"""FungalTraits ITS Functional Ecology module
Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MODE: {MODE}
SUBSET_MODE: {SUBSET_MODE}
COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}
P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}
SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}
Input FungalTraits DB: {FUNGAL_TRAITS_DB}
Output folder: {functional_dir.resolve()}
Output file: {annotated_path.name}
Final MetaDiv rows: {total:,}
Rows matched to FungalTraits genus: {matched:,}
Match percentage: {(matched / total * 100) if total else 0:.2f}%
FungalTraits columns exported: {len(ft_export_cols):,}

Notes:
- FINAL_DB keeps SPPN, OTU_XX and Original_ID for full feature traceability.
- Functional_Ecology exports SPPN as the user-facing feature identifier and removes OTU_XX/Original_ID.
- This Functional_Ecology file includes all available columns from FungalTraits.
- FungalTraits columns that had the same name as MetaDiv columns were prefixed with 'FungalTraits_'.
"""
    (functional_dir / "FungalTraits_Module_Log.txt").write_text(log_text, encoding="utf-8")

    print(f"   🍄 ITS Functional_Ecology full FungalTraits table saved: {annotated_path}")
    print(f"   • FungalTraits matches: {matched:,}/{total:,} rows ({(matched / total * 100) if total else 0:.1f}%)")
    print(f"   • FungalTraits columns exported: {len(ft_export_cols):,}")
    return final_df


# ============================================================================
# 16S FAPROTAX FUNCTIONAL ECOLOGY MODULE
# ============================================================================

FAPROTAX_RANK_COLUMNS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
FAPROTAX_OUTPUT_COLUMNS = [
    "faprotax_functions",
    "faprotax_function_count",
    "faprotax_matched_taxa",
    "faprotax_match_ranks",
    "faprotax_source",
]


def _normalize_faprotax_token(value) -> str:
    """
    Normalize taxon names for conservative FAPROTAX matching.

    The normalization is intentionally simple and transparent. It lowercases the
    string, removes confidence values and punctuation-like separators, and keeps
    only alphanumeric/underscore characters. This minimizes false positives in
    large meta-analysis datasets.
    """
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    text = re.sub(r"\([^)]*\)", "", text)
    text = text.replace(" ", "_")
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_").lower()
    if text in {"", "nan", "none", "na", "n_a", "unknown", "unclassified", "uncultured"}:
        return ""
    return text


def _parse_faprotax_group_header(line: str):
    """
    Parse a FAPROTAX group header.

    FAPROTAX group headers have the form:
        function_name    key:value; key:value

    Returns
    -------
    tuple[str, str]
        The functional group name and the raw metadata string.
    """
    parts = line.strip().split(None, 1)
    if not parts:
        return "", ""
    group_name = parts[0].strip()
    metadata = parts[1].strip() if len(parts) > 1 else ""
    return group_name, metadata


def parse_faprotax_database(faprotax_path: Path) -> tuple[list, dict]:
    """
    Read FAPROTAX.txt and convert it into matching rules.

    The official FAPROTAX database is organized as groups of functions followed
    by member taxa. Each member taxon is represented as a taxonomic path delimited
    by asterisks, for example:
        *Proteobacteria*Nitrosomonas*
        *Bacteria*Nitrosomonas*europaea*

    This parser also resolves FAPROTAX set operations:
        add_group:<function>
        subtract_group:<function>
        intersect_group:<function>

    Set operations are resolved in the order used by the database, following the
    FAPROTAX README. This makes the MetaDiv output closer to the official
    collapse_table.py behavior while avoiding the need to convert the table to
    BIOM or execute an external script.

    Returns
    -------
    rules : list[dict]
        Each rule contains function, tokens, terminal token and metadata.
    metadata_by_function : dict
        Raw metadata for each FAPROTAX function.
    """
    print("\n🌱 LOADING FAPROTAX DATABASE...")

    if not faprotax_path.exists():
        print(f"   ⚠️ FAPROTAX database not found: {faprotax_path}")
        print("   ⚠️ FAPROTAX module will be skipped.")
        return [], {}

    group_entries = defaultdict(list)
    group_order = []
    metadata_by_function = {}
    raw_rule_by_tokens = {}
    current_function = None
    current_metadata = ""

    with open(faprotax_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue

            lower = line.lower()
            if lower.startswith(("add_group:", "subtract_group:", "intersect_group:")):
                if current_function is None:
                    continue
                op, target = line.split(":", 1)
                group_entries[current_function].append((op.strip().lower(), target.strip(), None))
                continue

            if line.startswith("*"):
                if current_function is None:
                    continue
                tokens = tuple(_normalize_faprotax_token(t) for t in line.split("*") if t.strip())
                tokens = tuple(t for t in tokens if t)
                if not tokens:
                    continue
                group_entries[current_function].append(("rule", tokens, line))
                raw_rule_by_tokens.setdefault(tokens, line)
            else:
                current_function, current_metadata = _parse_faprotax_group_header(line)
                if current_function:
                    if current_function not in metadata_by_function:
                        group_order.append(current_function)
                    metadata_by_function[current_function] = current_metadata
                    group_entries.setdefault(current_function, [])

    resolved = {}
    for group in group_order:
        current_rules = set()
        for op, payload, raw_rule in group_entries.get(group, []):
            if op == "rule":
                current_rules.add(payload)
            elif op == "add_group":
                current_rules |= resolved.get(payload, set())
            elif op == "subtract_group":
                current_rules -= resolved.get(payload, set())
            elif op == "intersect_group":
                current_rules &= resolved.get(payload, set())
        resolved[group] = current_rules

    rules = []
    for function in group_order:
        for tokens in sorted(resolved.get(function, set())):
            if not tokens:
                continue
            rules.append({
                "function": function,
                "tokens": list(tokens),
                "terminal": tokens[-1],
                "metadata": metadata_by_function.get(function, ""),
                "raw_rule": raw_rule_by_tokens.get(tokens, "*" + "*".join(tokens) + "*"),
            })

    explicit_count = sum(1 for entries in group_entries.values() for e in entries if e[0] == "rule")
    print(f"   ✅ FAPROTAX functions detected: {len(metadata_by_function):,}")
    print(f"   ✅ Explicit taxonomic rules read: {explicit_count:,}")
    print(f"   ✅ Expanded function-rule associations: {len(rules):,}")
    return rules, metadata_by_function


def _row_taxonomy_tokens(row: pd.Series) -> tuple[list, dict]:
    """
    Build normalized taxonomic tokens from one MetaDiv row.

    Returns a list preserving rank order and a token->rank dictionary. Species are
    handled in two ways: as the full binomial produced by MetaDiv and, when
    possible, as the final epithet. This improves compatibility with FAPROTAX
    species-level rules such as *Genus*species*.
    """
    tokens = []
    token_to_rank = {}

    for rank in FAPROTAX_RANK_COLUMNS:
        value = row.get(rank, "")
        token = _normalize_faprotax_token(value)
        if token:
            tokens.append(token)
            token_to_rank.setdefault(token, rank)

        if rank == "species" and token and "_" in token:
            epithet = token.split("_")[-1]
            if epithet and epithet not in token_to_rank:
                tokens.append(epithet)
                token_to_rank[epithet] = "species"

    return tokens, token_to_rank


def _tokens_are_ordered_subset(rule_tokens: list, row_tokens: list) -> bool:
    """
    Check whether all rule tokens occur in the row taxonomy in the same order.

    FAPROTAX rules may omit intermediate ranks, so strict contiguous matching
    would be too restrictive. Ordered-subset matching is conservative enough to
    avoid most false positives while allowing paths such as phylum->genus.
    """
    if not rule_tokens or not row_tokens:
        return False
    pos = 0
    for rt in rule_tokens:
        found = False
        while pos < len(row_tokens):
            if row_tokens[pos] == rt:
                found = True
                pos += 1
                break
            pos += 1
        if not found:
            return False
    return True


def annotate_dataframe_with_faprotax(final_df: pd.DataFrame, rules: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Annotate the final MetaDiv 16S table with FAPROTAX functions.

    The function returns both the wide annotated database and a long-form table
    with one row per SPPN-function match. The long table is useful for auditing,
    filtering and building abundance profiles.
    """
    df = final_df.copy()

    if not rules:
        for col in FAPROTAX_OUTPUT_COLUMNS:
            df[col] = "unassigned" if col != "faprotax_function_count" else 0
        return df, pd.DataFrame()

    rules_by_terminal = defaultdict(list)
    for rule in rules:
        rules_by_terminal[rule["terminal"]].append(rule)

    function_strings = []
    function_counts = []
    matched_taxa_strings = []
    match_rank_strings = []
    long_records = []

    id_col = "SPPN" if "SPPN" in df.columns else "OTU_XX"

    print("\n🔎 MATCHING MetaDiv TAXA AGAINST FAPROTAX RULES...")

    for _, row in df.iterrows():
        row_tokens, token_to_rank = _row_taxonomy_tokens(row)
        candidate_rules = []

        # Candidate retrieval is efficient: only rules whose terminal taxon is
        # present in this row are evaluated.
        for token in set(row_tokens):
            candidate_rules.extend(rules_by_terminal.get(token, []))

        matched = {}
        matched_taxa = {}
        matched_ranks = {}
        matched_raw_rules = {}

        for rule in candidate_rules:
            if _tokens_are_ordered_subset(rule["tokens"], row_tokens):
                func = rule["function"]
                terminal = rule["terminal"]
                matched[func] = True
                matched_taxa[func] = terminal
                matched_ranks[func] = token_to_rank.get(terminal, "unknown")
                matched_raw_rules[func] = rule.get("raw_rule", "")

        if matched:
            funcs = sorted(matched.keys())
            taxa = [matched_taxa[f] for f in funcs]
            ranks = [matched_ranks[f] for f in funcs]
            function_strings.append(";".join(funcs))
            function_counts.append(len(funcs))
            matched_taxa_strings.append(";".join(taxa))
            match_rank_strings.append(";".join(ranks))

            for func in funcs:
                long_records.append({
                    id_col: row.get(id_col, ""),
                    "function": func,
                    "matched_taxon": matched_taxa[func],
                    "match_rank": matched_ranks[func],
                    "faprotax_rule": matched_raw_rules[func],
                })
        else:
            function_strings.append("unassigned")
            function_counts.append(0)
            matched_taxa_strings.append("unassigned")
            match_rank_strings.append("none")

    df["faprotax_functions"] = function_strings
    df["faprotax_function_count"] = function_counts
    df["faprotax_matched_taxa"] = matched_taxa_strings
    df["faprotax_match_ranks"] = match_rank_strings
    df["faprotax_source"] = "FAPROTAX.txt"

    long_df = pd.DataFrame(long_records)
    return df, long_df


def build_faprotax_functional_abundance(df: pd.DataFrame, long_df: pd.DataFrame, site_cols: list) -> pd.DataFrame:
    """
    Build a function-by-sample abundance table from FAPROTAX matches.

    If one taxon has multiple functions, its abundance contributes to each matched
    function. This follows the usual functional-screening interpretation of
    FAPROTAX-style assignments and keeps the output easy to analyze in R.
    """
    valid_sites = [c for c in site_cols if c in df.columns]
    if long_df.empty or not valid_sites:
        return pd.DataFrame(columns=["function"] + valid_sites)

    id_col = "SPPN" if "SPPN" in df.columns else "OTU_XX"
    abundance = df[[id_col] + valid_sites].copy()
    for c in valid_sites:
        abundance[c] = pd.to_numeric(abundance[c], errors="coerce").fillna(0)

    expanded = long_df[[id_col, "function"]].merge(abundance, on=id_col, how="left")
    profile = expanded.groupby("function", dropna=False)[valid_sites].sum().reset_index()
    return profile.sort_values("function").reset_index(drop=True)


def add_faprotax_functional_ecology_16s(final_df: pd.DataFrame, site_cols: list, output_dir: Path, faprotax_path: Path) -> pd.DataFrame:
    """
    Run the MetaDiv 16S FAPROTAX module.

    Outputs are saved in:
        output/16S/Functional_Ecology/

    Main outputs
    ------------
    - Final_Database_16S_FAPROTAX_annotated.csv
    - FAPROTAX_Functional_Abundance_Table.csv
    - FAPROTAX_Long_Format_Matches.csv
    - FAPROTAX_Unassigned_Taxa.csv
    - FAPROTAX_Module_Log.txt
    """
    if MODE != "16S" or not RUN_FAPROTAX_16S:
        return final_df

    functional_dir = output_dir / "Functional_Ecology"
    functional_dir.mkdir(parents=True, exist_ok=True)

    rules, metadata_by_function = parse_faprotax_database(faprotax_path)
    df, long_df = annotate_dataframe_with_faprotax(final_df, rules)

    assigned_mask = df["faprotax_function_count"].astype(int) > 0
    total = len(df)
    assigned = int(assigned_mask.sum())

    annotated_path = functional_dir / "Final_Database_16S_FAPROTAX_annotated.csv"
    df.to_csv(annotated_path, index=False, encoding="utf-8")

    long_path = functional_dir / "FAPROTAX_Long_Format_Matches.csv"
    long_df.to_csv(long_path, index=False, encoding="utf-8")

    profile = build_faprotax_functional_abundance(df, long_df, site_cols)
    profile_path = functional_dir / "FAPROTAX_Functional_Abundance_Table.csv"
    profile.to_csv(profile_path, index=False, encoding="utf-8")

    metadata_path = functional_dir / "FAPROTAX_Function_Metadata.csv"
    pd.DataFrame([
        {"function": k, "metadata": v} for k, v in sorted(metadata_by_function.items())
    ]).to_csv(metadata_path, index=False, encoding="utf-8")

    unassigned_cols = [
        c for c in [
            "SPPN", "OTU_XX", "Original_ID", "domain", "phylum", "class",
            "order", "family", "genus", "species", "Total_Abundance"
        ] if c in df.columns
    ]
    unassigned_path = functional_dir / "FAPROTAX_Unassigned_Taxa.csv"
    df.loc[~assigned_mask, unassigned_cols].to_csv(unassigned_path, index=False, encoding="utf-8")

    log_path = functional_dir / "FAPROTAX_Module_Log.txt"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("MetaDiv Builder - 16S FAPROTAX Functional Ecology Module\n")
        f.write("Author: Bernardo Águila, UNAM\n")
        f.write("=" * 70 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"FAPROTAX database used: {faprotax_path}\n")
        f.write(f"Parsed FAPROTAX functions: {len(metadata_by_function):,}\n")
        f.write(f"Parsed explicit matching rules: {len(rules):,}\n")
        f.write(f"Total MetaDiv features: {total:,}\n")
        f.write(f"Assigned features: {assigned:,}\n")
        f.write(f"Unassigned features: {total - assigned:,}\n")
        f.write(f"Assignment percentage: {(assigned / total * 100) if total else 0:.2f}%\n")
        f.write(f"Long-format function matches: {len(long_df):,}\n")
        f.write("\nMost frequent FAPROTAX functions among features:\n")
        if not long_df.empty:
            counts = long_df["function"].value_counts().head(50)
            for function, count in counts.items():
                f.write(f"  {function}: {int(count):,}\n")

    print("\n🌱 16S FAPROTAX MODULE COMPLETED")
    print(f"   • Assigned features: {assigned:,}/{total:,} ({(assigned / total * 100) if total else 0:.2f}%)")
    print(f"   • Annotated table: {annotated_path}")
    print(f"   • Functional abundance table: {profile_path}")
    print(f"   • Long-format matches: {long_path}")
    print(f"   • Module folder: {functional_dir.resolve()}")

    return df


# ============================================================================
# CO1 FUNCTIONAL ECOLOGY PLACEHOLDER
# ============================================================================

def add_co1_functional_ecology_placeholder(final_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """
    Create a placeholder Functional_Ecology export for CO1.

    CO1 is currently used in MetaDiv Builder as a taxonomic harmonization marker.
    A standardized CO1 ecological/functional reference database equivalent to
    FungalTraits or FAPROTAX is not implemented here. This function keeps the
    output architecture consistent and explicitly documents that functional
    ecological annotation was intentionally not inferred.
    """
    if MODE != "CO1":
        return final_df

    functional_dir = output_dir / "Functional_Ecology"
    functional_dir.mkdir(parents=True, exist_ok=True)

    placeholder = final_df.copy()

    placeholder["co1_functional_group"] = "not_available"
    placeholder["co1_trophic_group"] = "not_available"
    placeholder["co1_habitat_group"] = "not_available"
    placeholder["co1_functional_annotation_source"] = "not_available"
    placeholder["co1_functional_annotation_confidence"] = "not_available"

    out_file = functional_dir / "Final_Database_CO1_Functional_Ecology_placeholder.csv"

    placeholder.to_csv(
        out_file,
        index=False,
        encoding="utf-8"
    )

    log_file = functional_dir / "CO1_Functional_Ecology_README.txt"

    log_text = f"""MetaDiv Builder - CO1 Functional Ecology placeholder
Author: Bernardo Águila, UNAM
Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Status:
CO1 functional ecological annotation was not inferred.

Reason:
No standardized CO1 functional/ecological reference database equivalent to
FungalTraits or FAPROTAX is currently implemented in MetaDiv Builder.

Input:
MODE: {MODE}
SUBSET_MODE: {SUBSET_MODE}
COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}
P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}
SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}

Output:
{out_file.name}

Notes:
- CO1 taxonomic harmonization was performed normally.
- Functional columns were exported as placeholders.
- Future versions may support curated CO1 ecological annotations such as
  metazoan traits, trophic guilds, habitat categories, or arthropod-specific
  functional groups.
"""

    log_file.write_text(log_text, encoding="utf-8")

    print("\n🧬 CO1 FUNCTIONAL_ECOLOGY PLACEHOLDER CREATED")
    print(f"   • Placeholder table: {out_file}")
    print(f"   • README: {log_file}")

    return final_df

# ============================================================================
# PART 1: DATA CONCATENATION AND FILTERING (ALWAYS SAVES TO DISK)
# ============================================================================

def discover_datasets(input_dir: Path) -> dict:
    print("=== STEP 1: Automatically discovering datasets ===")
    datasets = defaultdict(dict)

    if not input_dir.exists():
        raise FileNotFoundError(
            f"Input folder not found for MODE={MODE}: {input_dir}. "
            f"Expected structure: input/ITS, input/16S, input/CO1."
        )

    for file_path in input_dir.glob("*"):
        filename = file_path.name
        if filename.endswith("_abundance.csv"):
            prefix = filename.replace("_abundance.csv", "")
            datasets[prefix]["abundance"] = file_path
        elif filename.endswith("_sequences.fasta"):
            prefix = filename.replace("_sequences.fasta", "")
            datasets[prefix]["sequences"] = file_path
        elif filename.endswith("_taxonomy.sintax.txt"):
            prefix = filename.replace("_taxonomy.sintax.txt", "")
            datasets[prefix]["taxonomy"] = file_path
        elif filename.endswith("_taxonomy.sintax"):
            prefix = filename.replace("_taxonomy.sintax", "")
            datasets[prefix]["taxonomy"] = file_path

    complete = {}
    required = {"abundance", "sequences", "taxonomy"}
    for prefix, files in datasets.items():
        missing = required - set(files.keys())
        if not missing:
            complete[prefix] = files
            print(f"✅ Complete dataset detected: {prefix}")
        else:
            print(f"⚠️ Incomplete dataset: {prefix} - missing: {missing}")

    print(f"\n📊 Total complete datasets found: {len(complete)}")
    return complete

def read_fasta_sequences(fasta_path: Path) -> dict:
    print(f"  Reading FASTA file: {fasta_path.name}")
    sequences = {}
    current_id = None
    current_sequence = []

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)
                current_id = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(current_sequence)

    print(f"    Sequences loaded: {len(sequences)}")
    return sequences

def read_taxonomy_sintax(taxonomy_path: Path) -> pd.DataFrame:
    print(f"  Reading SINTAX taxonomy: {taxonomy_path.name}")
    recs = []
    with open(taxonomy_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                recs.append({"ID": parts[0].strip(), "sintax_taxonomy": parts[1]})
    df = pd.DataFrame(recs)
    if not df.empty:
        df = df.set_index("ID")
    print(f"    Taxonomy records loaded: {len(df)}")
    return df

def process_single_dataset(prefix: str, files: dict, report_lines: list) -> pd.DataFrame:
    print(f"\n=== PROCESSING DATASET: {prefix} ===")

    print("  Reading abundance table (MASTER table)...")
    df_raw = pd.read_csv(files["abundance"], dtype=str, low_memory=False)
    id_col = df_raw.columns[0]
    print(f"    ID column detected: '{id_col}'")
    print(f"    Total rows in abundance table: {len(df_raw)}")
    
    abundance_df = df_raw.set_index(id_col)
    master_ids = set(abundance_df.index)
    print(f"    Master IDs count: {len(master_ids)}")

    print("  Reading SINTAX taxonomy...")
    taxonomy_df = read_taxonomy_sintax(files["taxonomy"])
    
    if not taxonomy_df.empty:
        merged_df = abundance_df.join(taxonomy_df, how="left")
        taxonomy_ids = set(taxonomy_df.index)
        missing_tax = master_ids - taxonomy_ids
        if missing_tax:
            warning = f"      ⚠️ WARNING: {len(missing_tax)} IDs in abundance missing from taxonomy file"
            print(warning)
            report_lines.append(f"  {prefix}: {warning}")
    else:
        merged_df = abundance_df.copy()
        merged_df["sintax_taxonomy"] = ""

    print("  Reading FASTA sequences...")
    sequences_dict = read_fasta_sequences(files["sequences"])
    merged_df["sequence"] = merged_df.index.map(sequences_dict)
    
    seq_ids = set(sequences_dict.keys())
    missing_seq = master_ids - seq_ids
    if missing_seq:
        warning = f"      ⚠️ WARNING: {len(missing_seq)} IDs in abundance missing from FASTA file"
        print(warning)
        report_lines.append(f"  {prefix}: {warning}")

    merged_df = merged_df.reset_index()
    original_id_name = merged_df.columns[0]
    merged_df = merged_df.rename(columns={original_id_name: "Original_ID"})
    merged_df["OTU_XX"] = [f"{prefix}_ID{str(i+1).zfill(5)}" for i in range(len(merged_df))]
    
    print("  Calculating and sorting by total abundance...")
    special = ["Original_ID", "sintax_taxonomy", "sequence", "OTU_XX"]
    site_columns = [c for c in merged_df.columns if c not in special]
    
    site_data = merged_df[site_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    merged_df["Total_Abundance"] = site_data.sum(axis=1)
    merged_df = merged_df.sort_values("Total_Abundance", ascending=False)
    
    final_cols = ["OTU_XX", "Original_ID", "sintax_taxonomy", "sequence"] + site_columns
    final_df = merged_df[final_cols].copy()
    
    # Save concatenated table
    concat_path = CONCAT_DIR / f"{prefix}_concatenated.csv"
    final_df.to_csv(concat_path, index=False)
    print(f"✅ Concatenated saved: {concat_path}")
    
    return final_df

def filter_and_save_datasets(concatenated_df: pd.DataFrame, prefix: str, report_lines: list) -> tuple:
    """
    Filter and save one concatenated dataset according to MODE and SUBSET_MODE.

    Supported marker modes
    ----------------------
    ITS:
        - only_fungi
        - all_eukaryotes

    16S:
        - only_bacteria
        - all_prokaryotes

    CO1:
        - only_metazoa
        - all_eukaryotes

    Notes
    -----
    CO1 is treated as a taxonomic harmonization marker. Functional ecological
    annotation is not inferred here because no standardized CO1 ecological
    reference database equivalent to FungalTraits or FAPROTAX is currently
    implemented.
    """
    print(f"=== FILTERING & SAVING DATASET: {prefix} ===")

    tax_lower = concatenated_df["sintax_taxonomy"].astype(str).str.lower()

    if MODE == "ITS":
        fungi_mask = (
            tax_lower.str.contains("d:fungi", na=False, case=False) |
            tax_lower.str.contains("k:fungi", na=False, case=False)
        )

        if SUBSET_MODE == "only_fungi":
            target_mask = fungi_mask
            target_name = "fungi"
            target_dir = FUNGI_DIR
        elif SUBSET_MODE == "all_eukaryotes":
            target_mask = pd.Series(True, index=concatenated_df.index)
            target_name = "all_eukaryotes"
            target_dir = ALL_EUK_DIR
        else:
            raise ValueError("For MODE='ITS', SUBSET_MODE must be 'only_fungi' or 'all_eukaryotes'.")

        all_name = "all_eukaryotes"
        all_dir = ALL_EUK_DIR
        other_dir = NON_ANNOTATED_DIR

    elif MODE == "16S":
        bacteria_mask = (
            tax_lower.str.contains("d:bacteria", na=False, case=False) |
            tax_lower.str.contains("k:bacteria", na=False, case=False)
        )

        archaea_mask = (
            tax_lower.str.contains("d:archaea", na=False, case=False) |
            tax_lower.str.contains("k:archaea", na=False, case=False)
        )

        if SUBSET_MODE == "only_bacteria":
            target_mask = bacteria_mask
            target_name = "only_bacteria"
            target_dir = BACTERIA_DIR
        elif SUBSET_MODE == "all_prokaryotes":
            target_mask = bacteria_mask | archaea_mask
            target_name = "all_prokaryotes"
            target_dir = ALL_PROK_DIR
        else:
            raise ValueError("For MODE='16S', SUBSET_MODE must be 'only_bacteria' or 'all_prokaryotes'.")

        all_name = "all_prokaryotes"
        all_dir = ALL_PROK_DIR
        other_dir = NON_ANNOTATED_DIR

    elif MODE == "CO1":
        eukaryota_mask = (
            tax_lower.str.contains("d:eukaryota", na=False, case=False) |
            tax_lower.str.contains("k:eukaryota", na=False, case=False)
        )

        metazoa_mask = build_co1_metazoa_mask(tax_lower)

        if SUBSET_MODE == "only_metazoa":
            target_mask = metazoa_mask
            target_name = "only_metazoa"
            target_dir = METAZOA_DIR
        elif SUBSET_MODE == "all_eukaryotes":
            target_mask = eukaryota_mask
            target_name = "all_eukaryotes"
            target_dir = CO1_DIR
        else:
            raise ValueError("For MODE='CO1', SUBSET_MODE must be 'only_metazoa' or 'all_eukaryotes'.")

        all_name = "all_eukaryotes"
        all_dir = CO1_DIR
        other_dir = NON_ANNOTATED_DIR

    else:
        raise ValueError("MODE must be 'ITS', '16S', or 'CO1'.")

    target_df = concatenated_df[target_mask].copy()
    other_df = concatenated_df[~target_mask].copy()

    target_path = target_dir / f"{prefix}_{target_name}.csv"
    other_path = other_dir / f"{prefix}_non_annotated.csv"
    all_path = all_dir / f"{prefix}_{all_name}.csv"

    target_df.to_csv(target_path, index=False)
    other_df.to_csv(other_path, index=False)
    concatenated_df.to_csv(all_path, index=False)

    stats = {
        "prefix": prefix,
        "total_ids": len(concatenated_df),
        "target_ids": len(target_df),
        "other_ids": len(other_df),
        "target_percentage": (len(target_df) / len(concatenated_df) * 100) if len(concatenated_df) else 0.0,
    }

    if MODE == "ITS":
        stats["fungi_ids"] = int(fungi_mask.sum())
        filter_desc = f"{SUBSET_MODE} ({stats['fungi_ids']} Fungi detected)"

    elif MODE == "16S":
        stats["bacteria_ids"] = int(bacteria_mask.sum())
        stats["archaea_ids"] = int(archaea_mask.sum())
        filter_desc = (
            f"{SUBSET_MODE} ({stats['bacteria_ids']} Bacteria, "
            f"{stats['archaea_ids']} Archaea detected)"
        )

    elif MODE == "CO1":
        stats["eukaryota_ids"] = int(eukaryota_mask.sum())
        stats["metazoa_ids"] = int(metazoa_mask.sum())
        filter_desc = (
            f"{SUBSET_MODE} ({stats['eukaryota_ids']} Eukaryota, "
            f"{stats['metazoa_ids']} Metazoa detected)"
        )

    else:
        filter_desc = SUBSET_MODE

    print(f"✅ {prefix}: {stats['target_ids']} {filter_desc} ({stats['target_percentage']:.1f}%)")
    report_lines.append(f"  {prefix}: {stats['target_ids']} {filter_desc} ({stats['target_percentage']:.1f}%)")

    return stats, target_df

# ============================================================================
# PART 2: TAXONOMIC COLLAPSE (READS FROM DISK)
# ============================================================================

def clean_taxonomy_string(tax_str):
    if not isinstance(tax_str, str):
        return tax_str
    return re.sub(r"([dkpcofgs]):([^(),]+)\([^)]+\)(?=\(\d+\.\d+\))", r"\1:\2", tax_str)

def extract_taxonomy_data(tax_str):
    if not isinstance(tax_str, str):
        return {}
    clean_str = clean_taxonomy_string(tax_str)
    tax_data = {}
    for prefix, taxon, p_value in re.findall(r"([dkpcofgs]):([^(]+)\(([\d.]+)\)", clean_str):
        prefix = prefix.lower()
        tax_data[prefix] = {"taxon": taxon.strip(), "p_value": float(p_value)}
    if "d" not in tax_data and "k" in tax_data:
        tax_data["d"] = tax_data["k"]
    return tax_data

def build_collapse_key_species_only(tax_data, p_thresh):
    if not tax_data or "s" not in tax_data:
        return None
    if tax_data["s"]["p_value"] < p_thresh:
        return None
    hierarchy = ["d", "p", "c", "o", "f", "g", "s"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == "s":
            break
    return "_".join(parts) if parts else None

def build_collapse_key_genus(tax_data, p_thresh):
    if not tax_data or "g" not in tax_data:
        return None
    if tax_data["g"]["p_value"] < p_thresh:
        return None
    hierarchy = ["d", "p", "c", "o", "f", "g"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == "g":
            break
    return "_".join(parts) if parts else None

def build_collapse_key_lowest_rank(tax_data, p_thresh):
    if not tax_data:
        return None
    ranked_levels = ["s", "g", "f", "o", "c", "p", "d"]
    chosen = None
    for level in ranked_levels:
        if level in tax_data and tax_data[level]["p_value"] >= p_thresh:
            chosen = level
            break
    if chosen is None:
        return None
    hierarchy = ["d", "p", "c", "o", "f", "g", "s"]
    parts = []
    for level in hierarchy:
        if level in tax_data:
            parts.append(f"{level}:{tax_data[level]['taxon']}")
        if level == chosen:
            break
    return "_".join(parts) if parts else None

def build_collapse_key(tax_data, strategy, p_thresh):
    if strategy == "species_only":
        return build_collapse_key_species_only(tax_data, p_thresh)
    elif strategy == "genus":
        return build_collapse_key_genus(tax_data, p_thresh)
    elif strategy == "all":
        return build_collapse_key_lowest_rank(tax_data, p_thresh)
    else:
        raise ValueError(f"Unsupported COLLAPSE_STRATEGY: {strategy}")

def parse_tax_columns_from_sintax(sintax_str):
    tax = extract_taxonomy_data(sintax_str)
    prefs = ["d", "p", "c", "o", "f", "g", "s"]
    out = {}
    for rank, pref in zip(TAX_RANKS, prefs):
        out[rank] = tax.get(pref, {}).get("taxon", "")
        out[f"{rank}_pvalue"] = tax.get(pref, {}).get("p_value", 0.0)
    return out

def has_complete_taxonomy(sintax_str):
    tax = extract_taxonomy_data(sintax_str)
    required = ["d", "p", "c", "o", "f", "g", "s"]
    for k in required:
        if k not in tax or not str(tax[k].get("taxon", "")).strip():
            return False
    return True

def get_deepest_confident_rank_label(tax_data, p_thresh):
    if not tax_data:
        return None
    if "s" in tax_data and tax_data["s"]["p_value"] >= p_thresh:
        return "species"
    if "g" in tax_data and tax_data["g"]["p_value"] >= p_thresh:
        return "genus"
    if "f" in tax_data and tax_data["f"]["p_value"] >= p_thresh:
        return "family"
    if "o" in tax_data and tax_data["o"]["p_value"] >= p_thresh:
        return "order"
    if "c" in tax_data and tax_data["c"]["p_value"] >= p_thresh:
        return "class"
    if "p" in tax_data and tax_data["p"]["p_value"] >= p_thresh:
        return "phylum"
    if "d" in tax_data and tax_data["d"]["p_value"] >= p_thresh:
        return "domain"
    return None

def create_binomial_species_inplace(df):
    created = 0
    no_genus = 0
    no_species = 0
    for idx, row in df.iterrows():
        genus = str(row.get("genus", "")).strip()
        species = str(row.get("species", "")).strip()
        if genus and species and genus != "nan" and species != "nan":
            df.at[idx, "species"] = f"{genus}_{species}"
            created += 1
        elif species and (not genus or genus == "nan"):
            no_genus += 1
        elif genus and (not species or species == "nan"):
            no_species += 1
    print(f"   • Binomials created: {created}")
    if no_genus or no_species:
        print(f"   • Warning: no genus for species: {no_genus} | no species for genus: {no_species}")

def sanitize_taxonomy_columns_by_threshold(df: pd.DataFrame, threshold: float, cascade_blank: bool = True) -> pd.DataFrame:
    print(f"\n🧽 SANITIZING TAXONOMY COLUMNS BY p≥{threshold:.2f} (cascade_blank={cascade_blank})...")
    ranks = TAX_RANKS
    for r in ranks:
        pcol = f"{r}_pvalue"
        if pcol in df.columns:
            df[pcol] = pd.to_numeric(df[pcol], errors="coerce").fillna(0.0)
        else:
            df[pcol] = 0.0
        if r not in df.columns:
            df[r] = ""
    blank_counts = {r: 0 for r in ranks}
    if not cascade_blank:
        for r in ranks:
            mask = df[f"{r}_pvalue"] < threshold
            blank_counts[r] = int((mask & (df[r].astype(str).str.strip() != "")).sum())
            df.loc[mask, r] = ""
    else:
        invalid_up_to = pd.Series(False, index=df.index)
        for r in ranks:
            invalid_up_to = invalid_up_to | (df[f"{r}_pvalue"] < threshold)
            mask = invalid_up_to
            blank_counts[r] = int((mask & (df[r].astype(str).str.strip() != "")).sum())
            df.loc[mask, r] = ""
    total_blank = sum(blank_counts.values())
    if total_blank:
        print("   • Blanked names per rank:")
        for r in ranks:
            if blank_counts[r]:
                print(f"     - {r}: {blank_counts[r]}")
    else:
        print("   • Nothing blanked.")
    return df

def assign_sppn_unique(final_df):
    print(f"\n🏷️ ASSIGNING FINAL SPPN CODES (PER TAXON, p≥{SPPN_P_THRESHOLD:.2f})...")
    if "SPPN" in final_df.columns:
        final_df = final_df.drop(columns=["SPPN"])
    sppn = []
    counter = {}
    rank_order = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    for _, row in final_df.iterrows():
        base_taxon = None
        for rank in rank_order:
            name = str(row.get(rank, "") or "").strip()
            pval = row.get(f"{rank}_pvalue", 0.0)
            try:
                pval = float(pval)
            except Exception:
                pval = 0.0
            if name and pval >= SPPN_P_THRESHOLD and name not in ["", "nan", "None"]:
                base_taxon = name
                break
        if not base_taxon:
            base_taxon = "Unknown"
        base = re.sub(r"[^\w]+", "_", base_taxon)
        counter[base] = counter.get(base, 0) + 1
        sppn.append(f"{base}_{counter[base]:04d}")
    final_df["SPPN"] = sppn
    return final_df

def prune_zero_abundance(final_df, site_cols):
    print("\n🧹 PRUNING ZERO-ABUNDANCE IDs...")
    site_data = final_df[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    final_df["Total_Abundance"] = site_data.sum(axis=1)
    zero_mask = final_df["Total_Abundance"] == 0
    n_zero = int(zero_mask.sum())
    if n_zero:
        print(f"   • Removing {n_zero} IDs with Total_Abundance == 0")
        final_df = final_df.loc[~zero_mask].reset_index(drop=True)
    else:
        print("   • No zeros to prune.")
    print(f"   • Remaining IDs: {len(final_df)}")
    return final_df

def load_datasets_with_stats_part2(report_lines):
    print("\n📁 PART 2 — LOADING FILES...")
    files = list(TARGET_DIR.glob(PATTERN))

    if not files:
        print(f"❌ No files found in: {TARGET_DIR} with pattern: {PATTERN}")
        return pd.DataFrame(), [], TARGET_DIR, PATTERN

    all_dfs = []
    site_union = set()
    fixed = {"OTU_XX", "Original_ID", "sintax_taxonomy", "sequence"}

    for csv_file in files:
        df = pd.read_csv(csv_file)
        subset = csv_file.stem.replace("_only_bacteria", "").replace("_all_prokaryotes", "").replace("_only_metazoa", "").replace("_fungi", "").replace("_all_eukaryotes", "")
        site_cols = [c for c in df.columns if c not in fixed]

        print(f"   🧹 Cleaning taxonomy strings for {subset}...")
        df["sintax_taxonomy"] = df["sintax_taxonomy"].fillna("").apply(clean_taxonomy_string)

        tx = df["sintax_taxonomy"].apply(extract_taxonomy_data)

        complete_tax = int(df["sintax_taxonomy"].apply(has_complete_taxonomy).sum())
        with_species = int(tx.apply(lambda x: ("s" in x) if x else False).sum())
        with_genus = int(tx.apply(lambda x: ("g" in x) if x else False).sum())
        confident_species = int(tx.apply(lambda x: ("s" in x and x["s"]["p_value"] >= P_VALUE_THRESHOLD) if x else False).sum())
        confident_genus = int(tx.apply(lambda x: ("g" in x and x["g"]["p_value"] >= P_VALUE_THRESHOLD) if x else False).sum())

        deepest_counts = {"species": 0, "genus": 0, "family": 0, "order": 0, "class": 0, "phylum": 0, "domain": 0, "none": 0}
        for tax_data in tx:
            label = get_deepest_confident_rank_label(tax_data, P_VALUE_THRESHOLD)
            if label is None:
                deepest_counts["none"] += 1
            else:
                deepest_counts[label] += 1

        report_lines.append(f"  {subset}: {len(df)} features, {len(site_cols)} sites")

        site_union.update(site_cols)
        all_dfs.append((subset, df))

    site_cols_all = sorted(site_union)
    unified = []

    for subset, df in all_dfs:
        current_site_cols = [c for c in df.columns if c not in fixed]
        missing = set(site_cols_all) - set(current_site_cols)
        if missing:
            df_missing = pd.DataFrame(0, index=df.index, columns=sorted(missing))
            df = pd.concat([df, df_missing], axis=1)
        df = df[["OTU_XX", "Original_ID", "sintax_taxonomy", "sequence"] + site_cols_all]
        unified.append(df)

    combined = pd.concat(unified, ignore_index=True, copy=False)

    print(f"\n📈 PART 2 — COMBINED SUMMARY")
    print(f"   • Total features: {len(combined):,}")
    print(f"   • Total sites: {len(site_cols_all):,}")
    return combined, site_cols_all, TARGET_DIR, PATTERN

def collapse_otus_simple(combined_df, site_cols, report_lines):
    print(f"\n🔬 APPLYING COLLAPSE STRATEGY: {COLLAPSE_STRATEGY} (p≥{P_VALUE_THRESHOLD})")

    tax_list = combined_df["sintax_taxonomy"].apply(extract_taxonomy_data)
    collapse_keys = tax_list.apply(lambda t: build_collapse_key(t, COLLAPSE_STRATEGY, P_VALUE_THRESHOLD))

    total = len(combined_df)
    collapsable = int(collapse_keys.notna().sum())
    non_collapsable = total - collapsable

    complete_in_collapsable = 0
    for idx in collapse_keys[collapse_keys.notna()].index:
        if has_complete_taxonomy(combined_df.loc[idx, "sintax_taxonomy"]):
            complete_in_collapsable += 1

    print("📊 FEATURE DISTRIBUTION:")
    print(f"   • Collapsable: {collapsable} ({(collapsable/total*100 if total else 0):.1f}%)")
    print(f"   • Complete among collapsable: {complete_in_collapsable} ({(complete_in_collapsable/collapsable*100 if collapsable else 0):.1f}%)")
    print(f"   • Individual: {non_collapsable} ({(non_collapsable/total*100 if total else 0):.1f}%)")

    report_lines.append(f"\n  Collapse Statistics:")
    report_lines.append(f"    Features before collapse: {total:,}")
    report_lines.append(f"    Collapsable groups: {collapsable:,}")
    report_lines.append(f"    Complete taxonomy in collapsable: {complete_in_collapsable:,}")

    result_rows = []
    logs = []

    group_index = {}
    for idx, key in collapse_keys.items():
        group_index.setdefault(key, []).append(idx)

    print(f"\n🔄 Processing {len([k for k in group_index if k is not None])} collapsable groups...")

    for key, idxs in group_index.items():
        grp = combined_df.loc[idxs]

        if key is None:
            for _, row in grp.iterrows():
                parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
                out = row.to_dict()
                out.update(parsed)
                result_rows.append(out)
            continue

        if len(grp) == 1:
            row = grp.iloc[0]
            parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
            out = row.to_dict()
            out.update(parsed)
            result_rows.append(out)
            continue

        grp = grp.copy()
        complete_mask = grp["sintax_taxonomy"].apply(has_complete_taxonomy)
        candidates = grp[complete_mask]

        if candidates.empty:
            for _, row in grp.iterrows():
                parsed = parse_tax_columns_from_sintax(row["sintax_taxonomy"])
                out = row.to_dict()
                out.update(parsed)
                result_rows.append(out)

            logs.append({
                "collapse_key": key, "kept_otu": None, "removed_count": 0,
                "removed_otus": [], "group_size": len(grp),
                "representative_seq_len": None, "representative_tax": None,
                "note": f"skip_collapse_incomplete_tax_{COLLAPSE_STRATEGY}",
            })
            continue

        seq_lengths = candidates["sequence"].astype(str).str.len().fillna(0)
        idx_keep = seq_lengths.idxmax()
        kept = candidates.loc[idx_keep].copy()

        summed = grp[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum()

        parsed_tax = parse_tax_columns_from_sintax(kept["sintax_taxonomy"])
        out = kept.to_dict()
        out.update(parsed_tax)

        for c, val in zip(site_cols, summed.values):
            out[c] = val
        result_rows.append(out)

        removed = [r for r in grp["OTU_XX"].tolist() if r != kept["OTU_XX"]]
        logs.append({
            "collapse_key": key, "kept_otu": kept["OTU_XX"], "removed_count": len(removed),
            "removed_otus": removed[:10], "group_size": len(grp),
            "representative_seq_len": int(seq_lengths.loc[idx_keep]),
            "representative_tax": kept["sintax_taxonomy"],
            "note": f"collapsed_with_complete_tax_{COLLAPSE_STRATEGY}",
        })

    final_df = pd.DataFrame(result_rows)
    print("\n✅ COLLAPSE COMPLETED")
    print(f"   • Initial features: {len(combined_df):,}")
    print(f"   • Final features:   {len(final_df):,}")
    print(f"   • Reduction: {len(combined_df) - len(final_df):,} features")

    final_df["sequence_lenght"] = final_df["sequence"].astype(str).str.len()
    
    report_lines.append(f"    Features after collapse: {len(final_df):,}")
    report_lines.append(f"    Reduction: {len(combined_df) - len(final_df):,} features")
    
    return final_df, logs

def compute_total_abundance_and_sort(final_df, site_cols):
    print("\n📊 COMPUTING TOTAL ABUNDANCE & SORTING...")
    site_data = final_df[site_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    final_df["Total_Abundance"] = site_data.sum(axis=1)
    final_df = final_df.sort_values("Total_Abundance", ascending=False).reset_index(drop=True)
    create_binomial_species_inplace(final_df)
    return final_df


# ============================================================================
# KRONA EXPORT MODULE
# ============================================================================

def _clean_krona_taxon(value, rank: str) -> str:
    """
    Clean a taxonomic label for Krona-compatible output.

    Krona text files are simple tab-separated files where the first column is a
    numeric magnitude and the remaining columns represent the hierarchical path.
    Empty labels are replaced by explicit unclassified labels to avoid broken
    taxonomic paths.
    """
    value = "" if pd.isna(value) else str(value).strip()
    bad_values = {"", "nan", "none", "null", "na", "unclassified", "uncultured"}
    if value.lower() in bad_values:
        return f"Unclassified_{rank}" if KRONA_FILL_UNCLASSIFIED else ""
    value = re.sub(r"[\t\r\n]+", " ", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value if value else (f"Unclassified_{rank}" if KRONA_FILL_UNCLASSIFIED else "")


def get_krona_terminal_rank() -> str:
    """
    Return the terminal taxonomic rank used for Krona aggregation.
    "FULL" and "TOP" both use the complete hierarchy available in TAX_RANKS.
    "GENUS" is the recommended default for large Excel-template workflows.
    """
    mode = str(KRONA_MODE).strip().lower()
    if mode in {"full", "top"}:
        return TAX_RANKS[-1]
    valid = set(TAX_RANKS)
    if mode not in valid:
        raise ValueError(
            f"Unsupported KRONA_MODE={KRONA_MODE}. Use FULL, TOP, "
            "PHYLUM, CLASS, ORDER, FAMILY, GENUS or SPECIES."
        )
    return mode


def get_krona_ranks() -> list:
    """
    Return the taxonomic ranks included in the Krona path.
    """
    terminal_rank = get_krona_terminal_rank()
    idx = TAX_RANKS.index(terminal_rank)
    return TAX_RANKS[:idx + 1]


def _prepare_krona_base_df(final_df: pd.DataFrame, site_cols: list) -> pd.DataFrame:
    """
    Create a compact working DataFrame with abundance and cleaned taxonomy columns.
    """
    df = final_df.copy()
    valid_sites = [c for c in site_cols if c in df.columns]

    if "Total_Abundance" in df.columns:
        abundance = pd.to_numeric(df["Total_Abundance"], errors="coerce").fillna(0)
    elif valid_sites:
        abundance = df[valid_sites].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    else:
        abundance = pd.Series(1, index=df.index)

    df = df.loc[:, [c for c in TAX_RANKS if c in df.columns]].copy()
    for rank in TAX_RANKS:
        if rank not in df.columns:
            df[rank] = ""
        df[rank] = df[rank].apply(lambda x: _clean_krona_taxon(x, rank))

    df["_krona_abundance"] = abundance.values
    if KRONA_MIN_ABUNDANCE is not None:
        df = df[df["_krona_abundance"] >= KRONA_MIN_ABUNDANCE].copy()

    return df


def _make_fixed_columns_krona_table(krona_table: pd.DataFrame, used_ranks: list) -> pd.DataFrame:
    """
    Return a table with the complete taxonomic column set.

    This is useful when an Excel template expects the same columns every time
    even if the main Krona export was aggregated to genus, family or another rank.
    Ranks below the chosen terminal rank are left blank rather than filled with
    artificial labels, so the user can clearly see the aggregation depth.
    """
    fixed = krona_table.copy()
    for rank in TAX_RANKS:
        if rank not in fixed.columns:
            fixed[rank] = ""
    return fixed[["abundance"] + TAX_RANKS]


def _apply_krona_row_limit(krona_table: pd.DataFrame) -> tuple:
    """
    Apply row limits for Excel-template safety.

    Returns
    -------
    limited_table : pandas.DataFrame
        The possibly truncated Krona table.
    limit_note : str
        Human-readable explanation of what happened.
    """
    original_rows = len(krona_table)
    mode = str(KRONA_MODE).strip().upper()

    if mode == "TOP":
        keep_n = min(int(KRONA_TOP_ROWS), original_rows)
        return krona_table.head(keep_n).copy(), (
            f"KRONA_MODE=TOP: retained top {keep_n:,} rows by abundance "
            f"from {original_rows:,} grouped rows."
        )

    if KRONA_MAX_ROWS is not None and original_rows > int(KRONA_MAX_ROWS):
        if str(KRONA_LIMIT_ACTION).strip().upper() == "TOP":
            keep_n = int(KRONA_MAX_ROWS)
            return krona_table.head(keep_n).copy(), (
                f"Row limit applied: retained top {keep_n:,} rows by abundance "
                f"from {original_rows:,} grouped rows because this exceeded "
                f"KRONA_MAX_ROWS={KRONA_MAX_ROWS:,}."
            )
        return krona_table, (
            f"WARNING_ONLY: Krona table has {original_rows:,} rows, exceeding "
            f"KRONA_MAX_ROWS={KRONA_MAX_ROWS:,}. No rows were removed."
        )

    return krona_table, f"No row limit applied: {original_rows:,} rows."


def build_krona_total_table(final_df: pd.DataFrame, site_cols: list, krona_dir: Path) -> pd.DataFrame:
    """
    Build a Krona-compatible total-abundance table aggregated across all samples.

    The export depth is controlled by KRONA_MODE. For large meta-analyses that
    are later rendered through an Excel template, KRONA_MODE="GENUS" and
    KRONA_MAX_ROWS=30000 are recommended.

    Main output format:
        abundance<TAB>domain<TAB>phylum<TAB>...<TAB>terminal_rank

    This file can also be rendered with KronaTools, for example:
        ktImportText Krona_Total_Abundance.tsv -o Krona_Total_Abundance.html
    """
    krona_dir.mkdir(parents=True, exist_ok=True)

    df = _prepare_krona_base_df(final_df, site_cols)
    original_features = len(df)
    used_ranks = get_krona_ranks()

    grouped = (
        df.groupby(used_ranks, dropna=False)["_krona_abundance"]
          .sum()
          .reset_index()
          .rename(columns={"_krona_abundance": "abundance"})
          .sort_values("abundance", ascending=False)
          .reset_index(drop=True)
    )

    grouped_rows_before_limit = len(grouped)
    krona_table, limit_note = _apply_krona_row_limit(grouped)

    krona_path = krona_dir / "Krona_Total_Abundance.tsv"
    krona_table[["abundance"] + used_ranks].to_csv(
        krona_path,
        sep="\t",
        index=False,
        header=False,
        encoding="utf-8"
    )

    readable_path = krona_dir / "Krona_Total_Abundance_with_header.tsv"
    krona_table[["abundance"] + used_ranks].to_csv(
        readable_path,
        sep="\t",
        index=False,
        header=True,
        encoding="utf-8"
    )

    fixed_path = None
    if KRONA_EXPORT_FIXED_COLUMNS:
        fixed = _make_fixed_columns_krona_table(krona_table, used_ranks)
        fixed_path = krona_dir / "Krona_Total_Abundance_fixed_columns_with_header.tsv"
        fixed.to_csv(
            fixed_path,
            sep="\t",
            index=False,
            header=True,
            encoding="utf-8"
        )

    print(f"   ✅ Krona total table: {krona_path}")
    print(f"   ✅ Krona readable table: {readable_path}")
    if fixed_path is not None:
        print(f"   ✅ Krona fixed-column readable table: {fixed_path}")
    print(f"   • KRONA_MODE: {KRONA_MODE}")
    print(f"   • Terminal rank: {used_ranks[-1]}")
    print(f"   • Features entering Krona module: {original_features:,}")
    print(f"   • Grouped rows before row-limit: {grouped_rows_before_limit:,}")
    print(f"   • Final Krona rows: {len(krona_table):,}")
    print(f"   • {limit_note}")

    # Store metadata for the info file and report.
    krona_table.attrs["used_ranks"] = used_ranks
    krona_table.attrs["original_features"] = original_features
    krona_table.attrs["grouped_rows_before_limit"] = grouped_rows_before_limit
    krona_table.attrs["limit_note"] = limit_note

    return krona_table


def build_krona_per_sample_tables(final_df: pd.DataFrame, site_cols: list, krona_dir: Path) -> int:
    """
    Optionally export one Krona-compatible text file per sample.

    This is disabled by default because it can create thousands of files or very
    large outputs in meta-analyses. Enable with RUN_KRONA_PER_SAMPLE_EXPORT=True.
    Per-sample exports use the same KRONA_MODE and row-limit logic as the total
    abundance table.
    """
    sample_dir = krona_dir / "per_sample"
    sample_dir.mkdir(parents=True, exist_ok=True)

    df = final_df.copy()
    valid_sites = [c for c in site_cols if c in df.columns]
    if not valid_sites:
        print("   ⚠️ No site columns found. Per-sample Krona export skipped.")
        return 0

    used_ranks = get_krona_ranks()
    for rank in TAX_RANKS:
        if rank not in df.columns:
            df[rank] = ""
        df[rank] = df[rank].apply(lambda x: _clean_krona_taxon(x, rank))

    exported = 0
    for sample in valid_sites:
        abund = pd.to_numeric(df[sample], errors="coerce").fillna(0)
        tmp = df.loc[abund >= KRONA_MIN_ABUNDANCE, used_ranks].copy()
        if tmp.empty:
            continue

        tmp.insert(0, "abundance", abund.loc[tmp.index].values)
        tmp = (
            tmp.groupby(used_ranks, dropna=False)["abundance"]
               .sum()
               .reset_index()
               .sort_values("abundance", ascending=False)
               .reset_index(drop=True)
        )
        tmp, _ = _apply_krona_row_limit(tmp)

        safe_sample = re.sub(r"[^\w.-]+", "_", str(sample))
        out_path = sample_dir / f"{safe_sample}_Krona.tsv"
        tmp[["abundance"] + used_ranks].to_csv(
            out_path,
            sep="\t",
            index=False,
            header=False,
            encoding="utf-8"
        )
        exported += 1

    print(f"   ✅ Per-sample Krona files exported: {exported:,}")
    return exported


def build_krona_info_file(krona_dir: Path, krona_table: pd.DataFrame, per_sample_files: int):
    """
    Write a short usage guide for KronaTools and Excel-template workflows.
    """
    used_ranks = krona_table.attrs.get("used_ranks", get_krona_ranks())
    original_features = krona_table.attrs.get("original_features", "NA")
    grouped_rows_before_limit = krona_table.attrs.get("grouped_rows_before_limit", "NA")
    limit_note = krona_table.attrs.get("limit_note", "NA")

    info = f"""# MetaDiv Builder Krona Export

Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MODE: {MODE}
SUBSET_MODE: {SUBSET_MODE}
COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}
P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}
SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}

Krona configuration:
- RUN_KRONA_EXPORT: {RUN_KRONA_EXPORT}
- KRONA_MODE: {KRONA_MODE}
- KRONA_MIN_ABUNDANCE: {KRONA_MIN_ABUNDANCE}
- KRONA_MAX_ROWS: {KRONA_MAX_ROWS}
- KRONA_TOP_ROWS: {KRONA_TOP_ROWS}
- KRONA_LIMIT_ACTION: {KRONA_LIMIT_ACTION}
- KRONA_EXPORT_FIXED_COLUMNS: {KRONA_EXPORT_FIXED_COLUMNS}
- Terminal rank used: {used_ranks[-1]}
- Path columns used: {', '.join(used_ranks)}

Krona table statistics:
- Features entering Krona module: {original_features}
- Grouped rows before row-limit: {grouped_rows_before_limit}
- Final rows in total table: {len(krona_table)}
- Row-limit note: {limit_note}

Main Krona input without header:
- Krona_Total_Abundance.tsv

Readable inspection table:
- Krona_Total_Abundance_with_header.tsv

Fixed-column inspection table for Excel templates:
- Krona_Total_Abundance_fixed_columns_with_header.tsv

Per-sample Krona files:
- {per_sample_files}

Render with KronaTools, if desired:
-------------------------------
ktImportText Krona_Total_Abundance.tsv -o Krona_Total_Abundance.html

Recommended settings for very large meta-analyses:
--------------------------------------------------
KRONA_MODE = "GENUS"
KRONA_MAX_ROWS = 30000
KRONA_LIMIT_ACTION = "TOP"

Notes:
------
This module does not render HTML directly from Python. It writes standard
Krona-compatible text files that can be used by KronaTools or external templates.
For Excel-template workflows, keeping the final Krona table below ~30,000 rows is
recommended to avoid unstable behavior in large visualizer templates.
"""
    with open(krona_dir / "KRONA_INFO.txt", "w", encoding="utf-8") as f:
        f.write(info)


def export_krona_package(final_df: pd.DataFrame, site_cols: list, output_dir: Path, report_lines: list):
    """
    Export Krona-compatible visualization tables.

    This module does not call KronaTools directly. It writes standard text files
    that can be rendered into interactive HTML plots with ktImportText or used by
    external Excel-based Krona templates.
    """
    if not RUN_KRONA_EXPORT:
        print("\n🌐 KRONA EXPORT skipped (RUN_KRONA_EXPORT=False).")
        return

    print("\n🌐 KRONA EXPORT MODULE")
    print("=" * 70)

    krona_dir = output_dir / "Krona"
    krona_dir.mkdir(parents=True, exist_ok=True)

    krona_table = build_krona_total_table(final_df, site_cols, krona_dir)
    per_sample_files = 0
    if RUN_KRONA_PER_SAMPLE_EXPORT:
        per_sample_files = build_krona_per_sample_tables(final_df, site_cols, krona_dir)

    build_krona_info_file(krona_dir, krona_table=krona_table, per_sample_files=per_sample_files)

    report_lines.append("\n  Krona Export:")
    report_lines.append(f"    Krona folder: {krona_dir.resolve()}")
    report_lines.append(f"    KRONA_MODE: {KRONA_MODE}")
    report_lines.append(f"    Terminal rank: {krona_table.attrs.get('used_ranks', get_krona_ranks())[-1]}")
    report_lines.append(f"    Grouped rows before row-limit: {krona_table.attrs.get('grouped_rows_before_limit', 'NA')}")
    report_lines.append(f"    Total Krona rows exported: {len(krona_table):,}")
    report_lines.append(f"    Row-limit note: {krona_table.attrs.get('limit_note', 'NA')}")
    report_lines.append(f"    Per-sample Krona files: {per_sample_files:,}")

    print(f"   ✅ Krona package ready: {krona_dir.resolve()}")

# ============================================================================
# PART 3: FOR_R EXPORT MODULE
# ============================================================================

def make_safe_sample_names(site_cols: list) -> tuple[list, pd.DataFrame]:
    """
    Deprecated helper retained for backward compatibility.

    v1.7.3-final keeps original sample names in For_R and relies on
    check.names = FALSE in R instead of writing sample_name_map.csv.
    """
    safe_names = []
    used = {}

    for original in site_cols:
        name = str(original).strip()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
        safe = re.sub(r"_+", "_", safe).strip("_")

        if not safe:
            safe = "Sample"
        if safe[0].isdigit():
            safe = f"S_{safe}"

        base = safe
        if base in used:
            used[base] += 1
            safe = f"{base}_{used[base]}"
        else:
            used[base] = 1

        safe_names.append(safe)

    sample_map = pd.DataFrame({
        "Original_SampleID": [str(c) for c in site_cols],
        "Safe_SampleID": safe_names,
    })
    return safe_names, sample_map

def export_for_r_package(final_csv_path: Path, suffix: str, site_cols: list):
    """
    Export MetaDiv Builder outputs for downstream ecological analyses in R.

    The function uses the site columns detected by the main pipeline instead of
    trying to rediscover samples from the final table. This avoids errors when
    sample names are numeric (1, 2, 3) or when functional annotation columns are
    also numeric.

    Generated files
    ---------------
    - abundance_table.csv
    - taxonomy_table.csv
    - sample_metadata.csv
    - sequences.fasta
    - FOR_R_INFO.txt
    """
    final_csv_path = Path(final_csv_path)
    print("\n" + "=" * 70)
    print(f"📦 FOR_R EXPORT: {suffix}")
    print("=" * 70)

    if not final_csv_path.exists():
        print(f"   ⚠️ Final database not found: {final_csv_path}")
        return

    df = pd.read_csv(final_csv_path, low_memory=False)

    id_col = "SPPN" if "SPPN" in df.columns else "OTU_XX"
    if id_col not in df.columns:
        raise ValueError("For_R export requires either SPPN or OTU_XX column.")

    valid_site_cols = [c for c in site_cols if c in df.columns]
    if not valid_site_cols:
        print("   ⚠️ No valid site columns found in final database. Exporting taxonomy and sequences only.")

    # Keep original sample names in For_R outputs. R can preserve numeric or unusual
    # names when files are read with check.names = FALSE.
    sample_names = [str(c) for c in valid_site_cols]

    outdir = build_for_r_output_dir(suffix)
    outdir.mkdir(parents=True, exist_ok=True)

    # Abundance table: rows = taxa, columns = samples.
    if valid_site_cols:
        abundance = df[[id_col] + valid_site_cols].copy()
        for col in valid_site_cols:
            abundance[col] = pd.to_numeric(abundance[col], errors="coerce").fillna(0).astype("int64")
        abundance = abundance.drop_duplicates(subset=[id_col], keep="first")
        abundance = abundance.set_index(id_col)
    else:
        abundance = df[[id_col]].drop_duplicates(subset=[id_col], keep="first").set_index(id_col)

    abundance_path = outdir / "abundance_table.csv"
    abundance.to_csv(abundance_path, encoding="utf-8", na_rep="0")

    # Taxonomy table, including marker-specific ecological annotations when present.
    tax_cols = [c for c in TAX_RANKS if c in df.columns]
    ecological_cols = [
        c for c in [
            "primary_lifestyle", "Secondary_lifestyle",
            "faprotax_functions", "faprotax_function_count",
            "faprotax_matched_taxa", "faprotax_match_ranks", "faprotax_source",
            "general_lifestyle", "secondary_lifestyle", "functional_group",
            "functional_confidence", "functional_source", "functional_match_rank",
        ] if c in df.columns
    ]

    taxonomy_cols = [id_col] + tax_cols + ecological_cols
    taxonomy = df[taxonomy_cols].copy()
    taxonomy = taxonomy.drop_duplicates(subset=[id_col], keep="first").set_index(id_col)

    taxonomy_table_path = outdir / "taxonomy_table.csv"
    taxonomy.to_csv(taxonomy_table_path, encoding="utf-8")

    # Sample metadata template uses the same original sample IDs as abundance_table.csv.
    metadata = pd.DataFrame({
        "sample_id": sample_names,
        "Group": "",
        "Site": "",
        "Treatment": "",
        "Latitude": "",
        "Longitude": "",
    })
    metadata_path = outdir / "sample_metadata.csv"
    metadata.to_csv(metadata_path, index=False, encoding="utf-8")

    # Representative sequences.
    fasta_path = outdir / "sequences.fasta"
    nseq = 0
    with open(fasta_path, "w", encoding="utf-8") as handle:
        if "sequence" in df.columns:
            seq_df = df[[id_col, "sequence"]].drop_duplicates(subset=[id_col], keep="first")
            for _, row in seq_df.iterrows():
                seq_id = str(row.get(id_col, "")).strip()
                seq = str(row.get("sequence", "")).strip().replace(" ", "")
                if not seq_id or seq_id.lower() in {"nan", "none"}:
                    continue
                if not seq or seq.lower() in {"nan", "none"}:
                    continue
                handle.write(f">{seq_id}\n")
                for i in range(0, len(seq), 80):
                    handle.write(seq[i:i+80] + "\n")
                nseq += 1

    # R loading instructions.
    info_path = outdir / "FOR_R_INFO.txt"
    info = f"""# MetaDiv Builder - For_R Export Package

Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MODE: {MODE}
SUBSET_MODE: {SUBSET_MODE}
COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}
P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}
SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}
Source final database: {final_csv_path}
ID column: {id_col}
Features exported: {len(abundance):,}
Samples exported: {len(valid_site_cols):,}

Folder naming:
{subdir_hint()}

Files generated:
- abundance_table.csv
- taxonomy_table.csv
- sample_metadata.csv
- sequences.fasta

Important note about sample names:
For_R keeps the original sample names exactly as they appear in FINAL_DB.
When reading the tables in R, use check.names = FALSE to prevent automatic renaming
of numeric or unusual sample names.

Minimal R example:
------------------
library(phyloseq)
library(Biostrings)

otu <- read.csv("abundance_table.csv", row.names = 1, check.names = FALSE)
tax <- read.csv("taxonomy_table.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("sample_metadata.csv", row.names = 1, check.names = FALSE)
seqs <- readDNAStringSet("sequences.fasta")

ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE),
  tax_table(as.matrix(tax)),
  sample_data(meta),
  refseq(seqs)
)
"""
    with open(info_path, "w", encoding="utf-8") as f:
        f.write(info)

    print(f"   ✅ abundance_table.csv: {abundance.shape[0]:,} features x {abundance.shape[1]:,} samples")
    print(f"   ✅ taxonomy_table.csv: {taxonomy.shape[0]:,} features")
    print(f"   ✅ sample_metadata.csv: {len(metadata):,} samples")
    print(f"   ✅ sequences.fasta: {nseq:,} sequences")
    print(f"   📁 For_R folder: {outdir.resolve()}")

def subdir_hint() -> str:
    """Return a human-readable description of the For_R folder scheme."""
    return f"For_R/<subset>/{build_for_r_tag(SUBSET_NAME)}/"

# ============================================================================
# SAVE FINAL OUTPUTS
# ============================================================================

def save_outputs(final_df, logs, site_cols, part2_input_dir, part2_pattern, report_lines) -> dict:
    produced = {}

    suffix_main = SUBSET_NAME
    out_csv_main = FINAL_DB_DIR / build_final_db_filename(suffix_main)

    rank_cols = TAX_RANKS
    pval_cols = [f"{r}_pvalue" for r in rank_cols]

    # Include lifestyle/function columns if present.
    lifestyle_cols = []
    if MODE == "ITS" and "primary_lifestyle" in final_df.columns:
        lifestyle_cols = ["primary_lifestyle", "Secondary_lifestyle"]
    if MODE == "16S" and "general_lifestyle" in final_df.columns:
        lifestyle_cols = [
            "general_lifestyle", "secondary_lifestyle", "functional_group",
            "functional_confidence", "functional_source", "functional_match_rank"
        ]
    if MODE == "16S" and "faprotax_functions" in final_df.columns:
        lifestyle_cols += [
            "faprotax_functions", "faprotax_function_count",
            "faprotax_matched_taxa", "faprotax_match_ranks", "faprotax_source"
        ]
    
    first_cols = (
        ["SPPN", "OTU_XX", "Original_ID"] + rank_cols + pval_cols +
        lifestyle_cols + ["sequence", "sequence_lenght", "Total_Abundance"]
    )
    existing = [c for c in first_cols if c in final_df.columns]
    final_order = existing + [c for c in site_cols if c in final_df.columns]

    final_df[final_order].to_csv(out_csv_main, index=False, encoding="utf-8")
    print(f"💾 Final database saved: {out_csv_main}")
    produced[suffix_main] = out_csv_main

    # For 16S with all_prokaryotes mode, also export only_bacteria subset
    if MODE == "16S" and SUBSET_MODE == "all_prokaryotes":
        bacteria_out = FINAL_DB_DIR / build_final_db_filename("only_bacteria")
        domain_series = final_df["domain"].astype(str).str.strip().str.lower()
        bacteria_df = final_df[domain_series.eq("bacteria")].copy()
        bacteria_df[final_order].to_csv(bacteria_out, index=False, encoding="utf-8")
        print(f"🦠 Bacteria-only subset saved: {bacteria_out} | Rows: {len(bacteria_df)}")
        produced["only_bacteria"] = bacteria_out
        report_lines.append(f"\n  Bacteria-only subset: {len(bacteria_df):,} features")

    # For ITS with all_eukaryotes mode, also export only_fungi subset
    if MODE == "ITS" and SUBSET_MODE == "all_eukaryotes":
        fungi_out = FINAL_DB_DIR / build_final_db_filename("only_fungi")
        domain_series = final_df["domain"].astype(str).str.strip().str.lower()
        fungi_df = final_df[domain_series.eq("fungi")].copy()
        fungi_df[final_order].to_csv(fungi_out, index=False, encoding="utf-8")
        print(f"🍄 Fungi-only subset saved: {fungi_out} | Rows: {len(fungi_df)}")
        produced["only_fungi"] = fungi_out
        report_lines.append(f"\n  Fungi-only subset: {len(fungi_df):,} features")


    # For CO1 with all_eukaryotes mode, also export only_metazoa subset
    if MODE == "CO1" and SUBSET_MODE == "all_eukaryotes":
        metazoa_out = FINAL_DB_DIR / build_final_db_filename("only_metazoa")
        tax_string_for_mask = (
            "k:" + final_df.get("domain", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "p:" + final_df.get("phylum", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "c:" + final_df.get("class", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "o:" + final_df.get("order", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "f:" + final_df.get("family", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "g:" + final_df.get("genus", pd.Series("", index=final_df.index)).astype(str).str.lower() + "," +
            "s:" + final_df.get("species", pd.Series("", index=final_df.index)).astype(str).str.lower()
        )
        metazoa_mask = build_co1_metazoa_mask(tax_string_for_mask)
        metazoa_df = final_df[metazoa_mask].copy()
        metazoa_df[final_order].to_csv(metazoa_out, index=False, encoding="utf-8")
        print(f"🧬 CO1 Metazoa-only subset saved: {metazoa_out} | Rows: {len(metazoa_df)}")
        produced["only_metazoa"] = metazoa_out
        report_lines.append(f"\n  CO1 Metazoa-only subset: {len(metazoa_df):,} features")

    log_file = FINAL_DB_DIR / build_collapse_log_filename(suffix_main)
    with open(log_file, "w", encoding="utf-8") as f:
        f.write(f"LOG — {COLLAPSE_STRATEGY} collapse (p≥{P_VALUE_THRESHOLD})\n")
        f.write(f"MODE: {MODE}\n")
        f.write(f"SUBSET_MODE: {SUBSET_MODE}\n")
        f.write("=" * 70 + "\n\n")
        for i, L in enumerate(logs, 1):
            f.write(f"[Group {i}] key: {L.get('collapse_key')}\n")
            f.write(f"  Representative: {L.get('kept_otu')}\n")
            f.write(f"  Representative length: {L.get('representative_seq_len')}\n")
            f.write(f"  Members: {L.get('group_size')} | Removed: {L.get('removed_count')}\n")
            if L.get('removed_otus'):
                f.write(f"  Removed (up to 10): {', '.join(L['removed_otus'])}\n")
            f.write(f"  Note: {L.get('note')}\n\n")
        f.write("SUMMARY\n")
        f.write(f"  Final features after pruning zeros: {len(final_df)}\n")

    print(f"📋 Collapse log saved: {log_file}")
    return produced

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append(f"MetaDiv Builder - Processing Report")
    report_lines.append(f"Author: Bernardo Águila, UNAM")
    report_lines.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("=" * 70)
    report_lines.append(f"MODE: {MODE}")
    report_lines.append(f"SUBSET_MODE: {SUBSET_MODE}")
    report_lines.append(f"COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    report_lines.append(f"P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}")
    report_lines.append(f"SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}")
    report_lines.append(f"DEV_MODE: {DEV_MODE} (temporary folder will be kept only if True)")
    report_lines.append(f"RUN_FUNGALTRAITS_ITS: {RUN_FUNGALTRAITS_ITS}")
    report_lines.append(f"RUN_FAPROTAX_16S: {RUN_FAPROTAX_16S}")
    report_lines.append("RUN_CO1_FUNCTIONAL_PLACEHOLDER: True if MODE == 'CO1' else False")
    report_lines.append(f"INPUT_ROOT: {INPUT_ROOT}")
    report_lines.append(f"INPUT_DIR_ACTIVE_MARKER: {INPUT_DIR}")
    report_lines.append(f"DATABASE_DIR: {DATABASE_DIR}")
    report_lines.append(f"ECOLOGICAL_REFERENCE_DB_DIR: {ECOLOGICAL_REFERENCE_DB_DIR}")
    report_lines.append(f"TAXONOMIC_REFERENCE_DB_DIR: {TAXONOMIC_REFERENCE_DB_DIR}")
    report_lines.append(f"FUNGAL_TRAITS_DB: {FUNGAL_TRAITS_DB}")
    report_lines.append(f"FAPROTAX_DB_16S: {FAPROTAX_DB_16S}")
    report_lines.append("=" * 70)
    
    print("=" * 70)
    print(f"🚀 MetaDiv Builder — {MODE} Pipeline")
    print("=" * 70)
    print(f"   MODE: {MODE}")
    print(f"   SUBSET_MODE: {SUBSET_MODE}")
    print(f"   DEV_MODE: {DEV_MODE}")
    print(f"   COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   Output directory: {OUTPUT_DIR.resolve()}")
    print(f"   Input root: {INPUT_ROOT.resolve()}")
    print(f"   Active input datasets folder: {INPUT_DIR.resolve()}")
    print(f"   Reference databases folder: {DATABASE_DIR.resolve()}")
    print(f"   Ecological reference DB folder: {ECOLOGICAL_REFERENCE_DB_DIR.resolve()}")
    print(f"   Taxonomic reference DB folder: {TAXONOMIC_REFERENCE_DB_DIR.resolve()}")
    print("=" * 70)

    validate_reference_databases()

    # PART 1: Discover and process datasets (always creates concatenated tables)
    print("\n🧩 PART 1 — LOADING AND CONCATENATING DATASETS")
    start_time = time.time()

    datasets = discover_datasets(INPUT_DIR)
    if not datasets:
        print("❌ No complete datasets found in 'input' folder")
        return

    all_stats = []
    all_target_dfs = []

    for prefix, files in datasets.items():
        try:
            concatenated_df = process_single_dataset(prefix, files, report_lines)
            stats, target_df = filter_and_save_datasets(concatenated_df, prefix, report_lines)
            all_stats.append(stats)
            all_target_dfs.append(target_df)
        except Exception as e:
            print(f"❌ Error processing dataset {prefix}: {e}")
            report_lines.append(f"\n  ERROR processing {prefix}: {e}")
            continue

    processing_time = time.time() - start_time
    
    total_ids = sum(s["total_ids"] for s in all_stats)
    total_target = sum(s["target_ids"] for s in all_stats)
    
    report_lines.append("\n" + "=" * 70)
    report_lines.append("PART 1 SUMMARY")
    report_lines.append("=" * 70)
    report_lines.append(f"Processed datasets: {len(all_stats)}")
    report_lines.append(f"Total features processed: {total_ids:,}")
    report_lines.append(f"Total {SUBSET_MODE}: {total_target:,} ({total_target/total_ids*100:.1f}%)")
    report_lines.append(f"Processing time: {processing_time:.2f} seconds")
    
    print("\n" + "=" * 60)
    print("📊 FINAL PROCESSING REPORT (PART 1)")
    print("=" * 60)
    print(f"📁 Processed datasets: {len(all_stats)}")
    print(f"🧬 Total features processed: {total_ids:,}")
    print(f"🎯 Total {SUBSET_MODE}: {total_target:,} ({total_target/total_ids*100 if total_ids else 0:.1f}%)")
    print(f"⏱️ Total time: {processing_time:.2f} seconds")

    # PART 2: Load filtered datasets from disk and collapse
    print("\n🧬 PART 2 — STARTING TAXONOMIC COLLAPSE")

    # Load from the saved concatenated tables (critical!)
    combined, site_cols, part2_input_dir, part2_pattern = load_datasets_with_stats_part2(report_lines)
    
    if combined.empty:
        print("❌ No data to collapse. Check that Part 1 produced files.")
        return

    final_df, logs = collapse_otus_simple(combined, site_cols, report_lines)
    final_df = compute_total_abundance_and_sort(final_df, site_cols)
    final_df = sanitize_taxonomy_columns_by_threshold(final_df, threshold=SPPN_P_THRESHOLD, cascade_blank=True)
    
    # Integrate fungal lifestyles if MODE is ITS.
    # The compact FINAL_DB receives only primary_lifestyle and Secondary_lifestyle.
    # A complete FungalTraits export is created later in Functional_Ecology.
    traits_df = pd.DataFrame()
    if MODE == "ITS" and FUNGAL_TRAITS_DB is not None:
        traits_df = load_fungal_traits_db(FUNGAL_TRAITS_DB)
        if not traits_df.empty and "genus" in final_df.columns:
            # Create a clean genus column for matching
            final_df["genus_clean"] = final_df["genus"].astype(str).str.strip().str.lower()
            # Left join only the compact lifestyle columns for FINAL_DB
            final_df = final_df.merge(
                traits_df[["genus_clean", "primary_lifestyle", "Secondary_lifestyle"]],
                on="genus_clean",
                how="left"
            )
            # Drop the temporary key column
            final_df = final_df.drop(columns=["genus_clean"])
            # Fill missing lifestyles with empty string
            final_df["primary_lifestyle"] = final_df["primary_lifestyle"].fillna("")
            final_df["Secondary_lifestyle"] = final_df["Secondary_lifestyle"].fillna("")
            matched = final_df["primary_lifestyle"].ne("").sum()
            print(f"   • Lifestyle info added for {matched} rows ({matched/len(final_df)*100:.1f}%)")
        else:
            # Add empty columns if no traits available or no genus column
            final_df["primary_lifestyle"] = ""
            final_df["Secondary_lifestyle"] = ""
            print("   • No lifestyle data added (missing traits or genus column)")
    else:
        # 16S functional ecology is added after SPPN assignment and zero pruning
        # so module-specific outputs include final stable identifiers.
        pass

    final_df = assign_sppn_unique(final_df)
    final_df = prune_zero_abundance(final_df, site_cols)

    if MODE == "ITS":
        final_df = add_fungaltraits_functional_ecology_its(
            final_df=final_df,
            output_dir=OUTPUT_DIR,
            traits_df=traits_df
        )

    if MODE == "16S":
        final_df = add_faprotax_functional_ecology_16s(
            final_df=final_df,
            site_cols=site_cols,
            output_dir=OUTPUT_DIR,
            faprotax_path=FAPROTAX_DB_16S
        )

    if MODE == "CO1":
        final_df = add_co1_functional_ecology_placeholder(
            final_df=final_df,
            output_dir=OUTPUT_DIR
        )

    produced_csvs = save_outputs(final_df, logs, site_cols, part2_input_dir, part2_pattern, report_lines)

    # KRONA: Export Krona-compatible visualization tables
    export_krona_package(final_df, site_cols, OUTPUT_DIR, report_lines)

    # PART 3: Export For_R packages
    print("\n🧬 PART 3 — EXPORTING For_R PACKAGES")

    if SUBSET_NAME in produced_csvs:
        export_for_r_package(produced_csvs[SUBSET_NAME], suffix=SUBSET_NAME, site_cols=site_cols)

    if MODE == "16S" and "only_bacteria" in produced_csvs:
        export_for_r_package(produced_csvs["only_bacteria"], suffix="only_bacteria", site_cols=site_cols)

    if MODE == "ITS" and "only_fungi" in produced_csvs:
        export_for_r_package(produced_csvs["only_fungi"], suffix="only_fungi", site_cols=site_cols)

    if MODE == "CO1" and "only_metazoa" in produced_csvs:
        export_for_r_package(produced_csvs["only_metazoa"], suffix="only_metazoa", site_cols=site_cols)

    # CLEANUP: Remove concatenated_tables folder if not in developer mode
    if not DEV_MODE and CONCAT_DIR.exists():
        print(f"\n🧹 Cleaning up temporary folder for active marker only: {CONCAT_DIR}")

        # Safety guard: never delete broad output folders.
        # Only remove output/<MODE>/concatenated_tables.
        expected_concat_dir = OUTPUT_ROOT / MODE / "concatenated_tables"
        if CONCAT_DIR.resolve() == expected_concat_dir.resolve():
            shutil.rmtree(CONCAT_DIR)
            print("   ✅ Temporary concatenated tables removed for active marker only.")
        else:
            print("   ⚠️ Cleanup skipped because CONCAT_DIR did not match the expected marker-specific path.")
            print(f"   CONCAT_DIR: {CONCAT_DIR.resolve()}")
            print(f"   Expected : {expected_concat_dir.resolve()}")

    # Final report
    report_lines.append("\n" + "=" * 70)
    report_lines.append("FINAL SUMMARY")
    report_lines.append("=" * 70)
    report_lines.append(f"Final features after collapse and pruning: {len(final_df):,}")
    report_lines.append(f"FINAL_DB folder: {FINAL_DB_DIR.resolve()}")
    report_lines.append(f"For_R folder: {FOR_R_DIR.resolve()}")
    if (MODE == "ITS" and RUN_FUNGALTRAITS_ITS) or (MODE == "16S" and RUN_FAPROTAX_16S) or MODE == "CO1":
        report_lines.append(f"Functional_Ecology folder: {(OUTPUT_DIR / 'Functional_Ecology').resolve()}")
    report_lines.append(f"Total execution time: {time.time() - start_time:.2f} seconds")
    report_lines.append("=" * 70)
    report_lines.append("\n✅ MetaDiv Builder pipeline completed successfully")
    
    report_path = OUTPUT_DIR / "METADIV_REPORT.txt"
    write_report(report_lines, report_path)
    
    print("\n🎉 ALL DONE")
    print(f"   MODE: {MODE}")
    print(f"   SUBSET_MODE: {SUBSET_MODE}")
    print(f"   COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   Final features: {len(final_df):,}")
    print(f"   FINAL_DB folder: {FINAL_DB_DIR.resolve()}")
    print(f"   For_R folder: {FOR_R_DIR.resolve()}")
    print(f"   Report saved: {report_path.resolve()}")

if __name__ == "__main__":
    main()