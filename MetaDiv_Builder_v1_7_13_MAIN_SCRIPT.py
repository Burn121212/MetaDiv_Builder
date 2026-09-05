#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
MetaDiv Builder v1.7.13 - Heterogeneous Metabarcoding Integrator
from OTUs/ASVs/zOTUs/MOTUs to Diversity

Authors: Bernardo Águila, UNAM
         Miguel F. Romero-Guiterrez, UNAM

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

import argparse
import pandas as pd
import numpy as np
import os
import shutil
import gc
import html as html_lib
from pathlib import Path
from collections import defaultdict
import time
import re
from datetime import datetime

# ============================================================================
# CONFIGURATION - COMMAND-LINE ARGUMENTS
# ============================================================================

# ----------------------------------------------------------------------------
# Valid option sets
# ----------------------------------------------------------------------------

VALID_MODES = {"ITS", "16S", "CO1"}

VALID_SUBSET_MODES = {
    "ITS": {
        "only_fungi",
        "all_eukaryotes",
    },
    "16S": {
        "only_bacteria",
        "all_prokaryotes",
    },
    "CO1": {
        "only_metazoa",
        "all_eukaryotes",
    },
}

DEFAULT_SUBSET_MODE = {
    "ITS": "all_eukaryotes",
    "16S": "all_prokaryotes",
    "CO1": "all_eukaryotes",
}

VALID_COLLAPSE_STRATEGIES = {
    "species_only",
    "genus",
    "all",
}

VALID_KRONA_MODES = {
    "FULL",
    "TOP",
    "PHYLUM",
    "CLASS",
    "ORDER",
    "FAMILY",
    "GENUS",
    "SPECIES",
}

VALID_KRONA_LIMIT_ACTIONS = {
    "TOP",
    "WARNING_ONLY",
}


# ----------------------------------------------------------------------------
# Argument parsing helpers
# ----------------------------------------------------------------------------

class ArgumentFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Formatter that preserves examples and also prints default values."""
    pass


def probability(value):
    """
    Parse a probability/confidence threshold.

    Accepted range: 0.0 to 1.0
    """
    try:
        parsed_value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid numeric value: {value}"
        )

    if parsed_value < 0.0 or parsed_value > 1.0:
        raise argparse.ArgumentTypeError(
            f"Value must be between 0.0 and 1.0, got: {parsed_value}"
        )

    return parsed_value


def positive_int(value):
    """
    Parse a positive integer.
    """
    try:
        parsed_value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid integer value: {value}"
        )

    if parsed_value <= 0:
        raise argparse.ArgumentTypeError(
            f"Value must be greater than 0, got: {parsed_value}"
        )

    return parsed_value


def nonnegative_int(value):
    """
    Parse a non-negative integer.
    """
    try:
        parsed_value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid integer value: {value}"
        )

    if parsed_value < 0:
        raise argparse.ArgumentTypeError(
            f"Value must be >= 0, got: {parsed_value}"
        )

    return parsed_value


def nonnegative_float(value):
    """
    Parse a non-negative float.
    """
    try:
        parsed_value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid numeric value: {value}"
        )

    if parsed_value < 0:
        raise argparse.ArgumentTypeError(
            f"Value must be >= 0, got: {parsed_value}"
        )

    return parsed_value


NUMPY_UINT_DTYPES = {
    "uint16": np.uint16,
    "uint32": np.uint32,
    "uint64": np.uint64,
}


def normalize_numpy_uint_dtype(value):
    """
    Normalize unsigned integer dtype names used for abundance matrices.
    """
    normalized = str(value).lower().replace("np.", "")
    if normalized not in NUMPY_UINT_DTYPES:
        valid = ", ".join(sorted(NUMPY_UINT_DTYPES))
        raise argparse.ArgumentTypeError(
            f"Invalid dtype '{value}'. Valid options: {valid}"
        )
    return normalized


def normalize_mode(value):
    """
    Normalize marker mode.

    Accepts COI as an alias for CO1.
    """
    normalized = str(value).upper()

    if normalized == "COI":
        normalized = "CO1"

    if normalized not in VALID_MODES:
        valid = ", ".join(sorted(VALID_MODES))
        raise argparse.ArgumentTypeError(
            f"Invalid mode '{value}'. Valid options: {valid}"
        )

    return normalized


def normalize_subset_mode(value):
    """
    Normalize subset mode.

    Allows either hyphens or underscores.
    """
    return str(value).lower().replace("-", "_")


def normalize_collapse_strategy(value):
    """
    Normalize and validate collapse strategy.
    """
    normalized = str(value).lower().replace("-", "_")

    if normalized not in VALID_COLLAPSE_STRATEGIES:
        valid = ", ".join(sorted(VALID_COLLAPSE_STRATEGIES))
        raise argparse.ArgumentTypeError(
            f"Invalid collapse strategy '{value}'. Valid options: {valid}"
        )

    return normalized


def normalize_krona_mode(value):
    """
    Normalize and validate Krona mode.
    """
    normalized = str(value).upper()

    if normalized not in VALID_KRONA_MODES:
        valid = ", ".join(sorted(VALID_KRONA_MODES))
        raise argparse.ArgumentTypeError(
            f"Invalid Krona mode '{value}'. Valid options: {valid}"
        )

    return normalized


def normalize_krona_limit_action(value):
    """
    Normalize and validate Krona limit action.
    """
    normalized = str(value).upper()

    if normalized not in VALID_KRONA_LIMIT_ACTIONS:
        valid = ", ".join(sorted(VALID_KRONA_LIMIT_ACTIONS))
        raise argparse.ArgumentTypeError(
            f"Invalid Krona limit action '{value}'. Valid options: {valid}"
        )

    return normalized


def resolve_path_argument(value):
    """
    Expand environment variables, expand '~', and resolve a path.

    The path does not need to exist yet.
    """
    if value is None:
        return None

    return Path(os.path.expandvars(str(value))).expanduser().resolve(strict=False)


def add_bool_argument(
    parser,
    positive_flag,
    negative_flag,
    dest,
    default,
    positive_help,
    negative_help,
):
    """
    Add a pair of boolean flags, for example:

        --run-krona-export
        --no-run-krona-export
    """
    group = parser.add_mutually_exclusive_group()

    group.add_argument(
        positive_flag,
        dest=dest,
        action="store_true",
        help=positive_help,
    )

    group.add_argument(
        negative_flag,
        dest=dest,
        action="store_false",
        help=negative_help,
    )

    parser.set_defaults(**{dest: default})


# ----------------------------------------------------------------------------
# Project directory resolution
# ----------------------------------------------------------------------------

def resolve_project_directory(configured_directory=None) -> Path:
    """
    Resolve the MetaDiv Builder project directory.

    Priority order:
        1. --project-dir command-line argument
        2. METADIV_PROJECT_DIR environment variable
        3. Search current working directory and parents for input/ or databases/
        4. Current working directory

    Returns
    -------
    pathlib.Path
        Absolute path to the project directory.
    """
    if configured_directory:
        project_directory = resolve_path_argument(configured_directory)
        if not project_directory.exists():
            raise FileNotFoundError(
                "The directory provided with --project-dir does not exist: "
                f"{project_directory}"
            )
        return project_directory

    environment_directory = os.environ.get("METADIV_PROJECT_DIR")
    if environment_directory:
        project_directory = resolve_path_argument(environment_directory)
        if not project_directory.exists():
            raise FileNotFoundError(
                "METADIV_PROJECT_DIR does not exist: "
                f"{project_directory}"
            )
        return project_directory

    current_directory = Path.cwd().resolve()

    for candidate in (current_directory, *current_directory.parents):
        if (candidate / "input").is_dir() or (candidate / "databases").is_dir():
            return candidate

    return current_directory


# ----------------------------------------------------------------------------
# Parser definition
# ----------------------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="metadiv_builder.py",
        formatter_class=ArgumentFormatter,
        description="""
MetaDiv Builder v1.7.13

Heterogeneous metabarcoding integrator for building standardized diversity
matrices from OTUs, ASVs, zOTUs, MOTUs, representative sequences, and SINTAX
taxonomy annotations.
        """,
        epilog="""
Examples:

  ITS, all eukaryotes, strict species-level collapse:
    python metadiv_builder.py \\
      --mode ITS \\
      --subset-mode all_eukaryotes \\
      --collapse-strategy species_only \\
      --p-value-threshold 1.0 \\
      --sppn-p-threshold 0.8

  ITS, fungi only, recursive lowest-rank collapse:
    python metadiv_builder.py \\
      --mode ITS \\
      --subset-mode only_fungi \\
      --collapse-strategy all \\
      --p-value-threshold 0.8

  16S, bacteria only, with FAPROTAX:
    python metadiv_builder.py \\
      --mode 16S \\
      --subset-mode only_bacteria \\
      --run-faprotax-16s

  CO1, metazoa only, no intermediate development files:
    python metadiv_builder.py \\
      --mode CO1 \\
      --subset-mode only_metazoa \\
      --no-dev-mode

  Use a custom project directory:
    python metadiv_builder.py \\
      --project-dir /path/to/metadiv_project \\
      --mode ITS
        """,
    )

    # ------------------------------------------------------------------------
    # Core pipeline options
    # ------------------------------------------------------------------------

    parser.add_argument(
        "-m",
        "--mode",
        type=normalize_mode,
        default="ITS",
        help=(
            "Marker mode. Options: ITS, 16S, CO1. "
            "COI is accepted as an alias for CO1."
        ),
    )

    parser.add_argument(
        "-s",
        "--subset-mode",
        type=normalize_subset_mode,
        default=None,
        help=(
            "Subset mode. Valid options depend on --mode. "
            "ITS: only_fungi, all_eukaryotes. "
            "16S: only_bacteria, all_prokaryotes. "
            "CO1: only_metazoa, all_eukaryotes. "
            "If omitted, a mode-specific default is used."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--dev-mode",
        negative_flag="--no-dev-mode",
        dest="dev_mode",
        default=False,
        positive_help=(
            "Keep intermediate concatenated tables and temporary files."
        ),
        negative_help=(
            "Delete intermediate temporary files after successful completion."
        ),
    )

    parser.add_argument(
        "-c",
        "--collapse-strategy",
        type=normalize_collapse_strategy,
        default="species_only",
        help=(
            "Taxonomic collapse strategy. "
            "Options: species_only, genus, all."
        ),
    )

    parser.add_argument(
        "-p",
        "--p-value-threshold",
        "--p_value_threshold",
        dest="p_value_threshold",
        type=probability,
        default=1.00,
        help=(
            "Taxonomic collapse confidence threshold, from 0.0 to 1.0."
        ),
    )

    parser.add_argument(
        "--sppn-p-threshold",
        "--sppn_p_threshold",
        dest="sppn_p_threshold",
        type=probability,
        default=0.80,
        help=(
            "SPPN assignment confidence threshold, from 0.0 to 1.0."
        ),
    )

    # ------------------------------------------------------------------------
    # Directory and file path options
    # ------------------------------------------------------------------------

    parser.add_argument(
        "--project-dir",
        default=None,
        help=(
            "MetaDiv Builder project directory. "
            "If omitted, METADIV_PROJECT_DIR is used if available; otherwise "
            "the current directory and its parents are searched for input/ or "
            "databases/."
        ),
    )

    parser.add_argument(
        "--input-root",
        default=None,
        help=(
            "Root input directory. Default: PROJECT_DIR/input."
        ),
    )

    parser.add_argument(
        "--input-dir",
        default=None,
        help=(
            "Marker-specific input directory. "
            "If omitted, INPUT_ROOT/MODE is used."
        ),
    )

    parser.add_argument(
        "--database-dir",
        default=None,
        help=(
            "Database directory. Default: PROJECT_DIR/databases."
        ),
    )

    parser.add_argument(
        "--ecological-reference-db-dir",
        default=None,
        help=(
            "Ecological reference database directory. "
            "Default: DATABASE_DIR/ecological_reference_db."
        ),
    )

    parser.add_argument(
        "--taxonomic-reference-db-dir",
        default=None,
        help=(
            "Taxonomic reference database directory. "
            "Default: DATABASE_DIR/taxonomic_reference_db."
        ),
    )

    parser.add_argument(
        "--output-root",
        default=None,
        help=(
            "Output root directory. Default: PROJECT_DIR/output."
        ),
    )

    # ------------------------------------------------------------------------
    # Reference database filenames
    # ------------------------------------------------------------------------

    parser.add_argument(
        "--fungal-traits-filename",
        default="Fungal_Traits_DB.txt",
        help=(
            "Filename for the ITS FungalTraits database inside "
            "ECOLOGICAL_REFERENCE_DB_DIR. Ignored if --fungal-traits-db is supplied."
        ),
    )

    parser.add_argument(
        "--faprotax-filename",
        default="FAPROTAX.txt",
        help=(
            "Filename for the 16S FAPROTAX database inside "
            "ECOLOGICAL_REFERENCE_DB_DIR. Ignored if --faprotax-db/--faprotax-db-16s is supplied."
        ),
    )

    # ------------------------------------------------------------------------
    # Functional ecology modules
    # ------------------------------------------------------------------------

    add_bool_argument(
        parser=parser,
        positive_flag="--run-fungaltraits-its",
        negative_flag="--no-run-fungaltraits-its",
        dest="run_fungaltraits_its",
        default=True,
        positive_help=(
            "Run the optional ITS FungalTraits functional ecology module."
        ),
        negative_help=(
            "Disable the ITS FungalTraits functional ecology module."
        ),
    )

    parser.add_argument(
        "--fungal-traits-db",
        default=None,
        help=(
            "Path to Fungal_Traits_DB.txt. "
            "Default: ECOLOGICAL_REFERENCE_DB_DIR/Fungal_Traits_DB.txt."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--run-faprotax-16s",
        negative_flag="--no-run-faprotax-16s",
        dest="run_faprotax_16s",
        default=False,
        positive_help=(
            "Run the optional 16S FAPROTAX functional ecology module."
        ),
        negative_help=(
            "Disable the 16S FAPROTAX functional ecology module."
        ),
    )

    parser.add_argument(
        "--faprotax-db-16s",
        "--faprotax-db",
        dest="faprotax_db_16s",
        default=None,
        help=(
            "Path to FAPROTAX.txt. "
            "Default: ECOLOGICAL_REFERENCE_DB_DIR/FAPROTAX.txt."
        ),
    )

    # ------------------------------------------------------------------------
    # Krona export module
    # ------------------------------------------------------------------------

    add_bool_argument(
        parser=parser,
        positive_flag="--run-krona-export",
        negative_flag="--no-run-krona-export",
        dest="run_krona_export",
        default=True,
        positive_help=(
            "Write Krona-compatible output tables."
        ),
        negative_help=(
            "Disable Krona-compatible output tables."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--run-krona-per-sample-export",
        negative_flag="--no-run-krona-per-sample-export",
        dest="run_krona_per_sample_export",
        default=False,
        positive_help=(
            "Write per-sample Krona-compatible output tables."
        ),
        negative_help=(
            "Disable per-sample Krona-compatible output tables."
        ),
    )

    parser.add_argument(
        "--krona-mode",
        type=normalize_krona_mode,
        default="GENUS",
        help=(
            "Krona aggregation mode. Options: FULL, TOP, PHYLUM, CLASS, "
            "ORDER, FAMILY, GENUS, SPECIES."
        ),
    )

    parser.add_argument(
        "--krona-min-abundance",
        type=nonnegative_float,
        default=1,
        help=(
            "Minimum abundance retained in Krona-compatible output tables."
        ),
    )

    parser.add_argument(
        "--krona-max-rows",
        type=positive_int,
        default=30000,
        help=(
            "Maximum number of rows allowed before applying the Krona "
            "row-limit action."
        ),
    )

    parser.add_argument(
        "--krona-top-rows",
        type=positive_int,
        default=30000,
        help=(
            "Number of top-abundance rows retained when KRONA_LIMIT_ACTION "
            "is TOP."
        ),
    )

    parser.add_argument(
        "--krona-limit-action",
        type=normalize_krona_limit_action,
        default="TOP",
        help=(
            "Action if Krona table exceeds KRONA_MAX_ROWS. "
            "Options: TOP, WARNING_ONLY."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--krona-fill-unclassified",
        negative_flag="--no-krona-fill-unclassified",
        dest="krona_fill_unclassified",
        default=True,
        positive_help=(
            "Fill empty taxonomic ranks with unclassified labels in Krona output."
        ),
        negative_help=(
            "Do not fill empty taxonomic ranks in Krona output."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--krona-export-fixed-columns",
        negative_flag="--no-krona-export-fixed-columns",
        dest="krona_export_fixed_columns",
        default=True,
        positive_help=(
            "Export Krona tables using a fixed-column layout."
        ),
        negative_help=(
            "Export Krona tables without forcing the fixed-column layout."
        ),
    )

    # ------------------------------------------------------------------------
    # For_R and secondary subset exports
    # ------------------------------------------------------------------------

    add_bool_argument(
        parser=parser,
        positive_flag="--run-for-r-export",
        negative_flag="--no-run-for-r-export",
        dest="run_for_r_export",
        default=True,
        positive_help=(
            "Write For_R export packages for downstream MetaDiv R scripts."
        ),
        negative_help=(
            "Disable For_R export packages."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--run-secondary-subset-exports",
        negative_flag="--no-run-secondary-subset-exports",
        dest="run_secondary_subset_exports",
        default=True,
        positive_help=(
            "When a broad subset is selected, also export the nested subset "
            "for that marker, for example all_prokaryotes -> only_bacteria."
        ),
        negative_help=(
            "Do not write secondary nested subset CSVs. Useful for very large wide tables."
        ),
    )

    # ------------------------------------------------------------------------
    # ReInput module
    # ------------------------------------------------------------------------

    add_bool_argument(
        parser=parser,
        positive_flag="--run-reinput-export",
        negative_flag="--no-run-reinput-export",
        dest="run_reinput_export",
        default=False,
        positive_help=(
            "Write a MetaDiv-compatible ReInput package from FINAL_DB."
        ),
        negative_help=(
            "Disable ReInput package export."
        ),
    )

    # ------------------------------------------------------------------------
    # Memory / performance options
    # ------------------------------------------------------------------------

    add_bool_argument(
        parser=parser,
        positive_flag="--ram-safe-mode",
        negative_flag="--no-ram-safe-mode",
        dest="ram_safe_mode",
        default=False,
        positive_help=(
            "Enable RAM-safe mode flag in reports. The v1.7.13 CLI uses "
            "the RAM-safe implementations by default."
        ),
        negative_help=(
            "Disable the RAM-safe mode flag in reports."
        ),
    )

    parser.add_argument(
        "--ram-safe-chunk-rows",
        type=positive_int,
        default=20000,
        help=(
            "Row chunk size used by RAM-safe CSV reading and row compaction."
        ),
    )

    parser.add_argument(
        "--ram-safe-site-chunk",
        type=positive_int,
        default=32,
        help=(
            "Number of sample/abundance columns processed at once in RAM-safe operations."
        ),
    )

    parser.add_argument(
        "--ram-safe-group-rows",
        type=positive_int,
        default=20000,
        help=(
            "Maximum number of rows from one collapse group materialized at once."
        ),
    )

    parser.add_argument(
        "--ram-safe-abundance-dtype",
        type=normalize_numpy_uint_dtype,
        default="uint32",
        help=(
            "Unsigned integer dtype used for abundance matrices. Options: uint16, uint32, uint64."
        ),
    )

    parser.add_argument(
        "--ram-safe-sum-dtype",
        type=normalize_numpy_uint_dtype,
        default="uint64",
        help=(
            "Unsigned integer dtype used for abundance sums. Options: uint32, uint64."
        ),
    )

    add_bool_argument(
        parser=parser,
        positive_flag="--run-stage-timing",
        negative_flag="--no-run-stage-timing",
        dest="run_stage_timing",
        default=True,
        positive_help=(
            "Write detailed stage timing entries to the report."
        ),
        negative_help=(
            "Disable detailed stage timing report entries."
        ),
    )

    # ------------------------------------------------------------------------
    # Reporting option
    # ------------------------------------------------------------------------

    parser.add_argument(
        "--quiet-config",
        action="store_true",
        help=(
            "Do not print the resolved configuration at startup."
        ),
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------------
    # Mode-specific validation
    # ------------------------------------------------------------------------

    if args.subset_mode is None:
        args.subset_mode = DEFAULT_SUBSET_MODE[args.mode]

    if args.subset_mode not in VALID_SUBSET_MODES[args.mode]:
        valid = ", ".join(sorted(VALID_SUBSET_MODES[args.mode]))
        parser.error(
            f"Invalid --subset-mode '{args.subset_mode}' for --mode {args.mode}. "
            f"Valid options are: {valid}"
        )

    return args


# ----------------------------------------------------------------------------
# Parse command-line arguments
# ----------------------------------------------------------------------------

ARGS = parse_arguments()


# ----------------------------------------------------------------------------
# Core configuration values used by the rest of the pipeline
# ----------------------------------------------------------------------------

MODE = ARGS.mode
SUBSET_MODE = ARGS.subset_mode
DEV_MODE = ARGS.dev_mode
COLLAPSE_STRATEGY = ARGS.collapse_strategy
P_VALUE_THRESHOLD = ARGS.p_value_threshold
SPPN_P_THRESHOLD = ARGS.sppn_p_threshold
FUNGAL_TRAITS_FILENAME = ARGS.fungal_traits_filename
FAPROTAX_FILENAME = ARGS.faprotax_filename
RUN_FOR_R_EXPORT = ARGS.run_for_r_export
RUN_SECONDARY_SUBSET_EXPORTS = ARGS.run_secondary_subset_exports
RAM_SAFE_MODE = ARGS.ram_safe_mode
RAM_SAFE_CHUNK_ROWS = ARGS.ram_safe_chunk_rows
RAM_SAFE_SITE_CHUNK = ARGS.ram_safe_site_chunk
RAM_SAFE_GROUP_ROWS = ARGS.ram_safe_group_rows
RAM_SAFE_ABUNDANCE_DTYPE = NUMPY_UINT_DTYPES[ARGS.ram_safe_abundance_dtype]
RAM_SAFE_SUM_DTYPE = NUMPY_UINT_DTYPES[ARGS.ram_safe_sum_dtype]
RUN_STAGE_TIMING = ARGS.run_stage_timing


# ----------------------------------------------------------------------------
# Directory configuration
# ----------------------------------------------------------------------------

PROJECT_DIR = resolve_project_directory(ARGS.project_dir)

INPUT_ROOT = (
    resolve_path_argument(ARGS.input_root)
    if ARGS.input_root
    else PROJECT_DIR / "input"
)

INPUT_DIR = (
    resolve_path_argument(ARGS.input_dir)
    if ARGS.input_dir
    else INPUT_ROOT / MODE
)

DATABASE_DIR = (
    resolve_path_argument(ARGS.database_dir)
    if ARGS.database_dir
    else PROJECT_DIR / "databases"
)

ECOLOGICAL_REFERENCE_DB_DIR = (
    resolve_path_argument(ARGS.ecological_reference_db_dir)
    if ARGS.ecological_reference_db_dir
    else DATABASE_DIR / "ecological_reference_db"
)

TAXONOMIC_REFERENCE_DB_DIR = (
    resolve_path_argument(ARGS.taxonomic_reference_db_dir)
    if ARGS.taxonomic_reference_db_dir
    else DATABASE_DIR / "taxonomic_reference_db"
)

OUTPUT_ROOT = (
    resolve_path_argument(ARGS.output_root)
    if ARGS.output_root
    else PROJECT_DIR / "output"
)


# ----------------------------------------------------------------------------
# ITS FungalTraits Functional Ecology module
# ----------------------------------------------------------------------------

RUN_FUNGALTRAITS_ITS = ARGS.run_fungaltraits_its

FUNGAL_TRAITS_DB = (
    resolve_path_argument(ARGS.fungal_traits_db)
    if ARGS.fungal_traits_db
    else ECOLOGICAL_REFERENCE_DB_DIR / FUNGAL_TRAITS_FILENAME
)


# ----------------------------------------------------------------------------
# 16S FAPROTAX Functional Ecology module
# ----------------------------------------------------------------------------

RUN_FAPROTAX_16S = ARGS.run_faprotax_16s

FAPROTAX_DB_16S = (
    resolve_path_argument(ARGS.faprotax_db_16s)
    if ARGS.faprotax_db_16s
    else ECOLOGICAL_REFERENCE_DB_DIR / FAPROTAX_FILENAME
)


# ----------------------------------------------------------------------------
# Krona export module
# ----------------------------------------------------------------------------

RUN_KRONA_EXPORT = ARGS.run_krona_export
RUN_KRONA_PER_SAMPLE_EXPORT = ARGS.run_krona_per_sample_export
KRONA_MODE = ARGS.krona_mode
KRONA_MIN_ABUNDANCE = ARGS.krona_min_abundance
KRONA_MAX_ROWS = ARGS.krona_max_rows
KRONA_TOP_ROWS = ARGS.krona_top_rows
KRONA_LIMIT_ACTION = ARGS.krona_limit_action
KRONA_FILL_UNCLASSIFIED = ARGS.krona_fill_unclassified
KRONA_EXPORT_FIXED_COLUMNS = ARGS.krona_export_fixed_columns


# ----------------------------------------------------------------------------
# ReInput export module
# ----------------------------------------------------------------------------

RUN_REINPUT_EXPORT = ARGS.run_reinput_export


# ----------------------------------------------------------------------------
# Optional configuration checks and startup report
# ----------------------------------------------------------------------------

def print_configuration():
    """
    Print resolved MetaDiv Builder configuration.
    """
    print("=" * 80)
    print("MetaDiv Builder configuration")
    print("=" * 80)

    print("\nCore settings")
    print(f"  MODE:                         {MODE}")
    print(f"  SUBSET_MODE:                  {SUBSET_MODE}")
    print(f"  DEV_MODE:                     {DEV_MODE}")
    print(f"  RUN_STAGE_TIMING:             {RUN_STAGE_TIMING}")
    print(f"  COLLAPSE_STRATEGY:            {COLLAPSE_STRATEGY}")
    print(f"  P_VALUE_THRESHOLD:            {P_VALUE_THRESHOLD}")
    print(f"  SPPN_P_THRESHOLD:             {SPPN_P_THRESHOLD}")

    print("\nDirectories")
    print(f"  PROJECT_DIR:                  {PROJECT_DIR}")
    print(f"  INPUT_ROOT:                   {INPUT_ROOT}")
    print(f"  INPUT_DIR:                    {INPUT_DIR}")
    print(f"  DATABASE_DIR:                 {DATABASE_DIR}")
    print(f"  ECOLOGICAL_REFERENCE_DB_DIR:  {ECOLOGICAL_REFERENCE_DB_DIR}")
    print(f"  TAXONOMIC_REFERENCE_DB_DIR:   {TAXONOMIC_REFERENCE_DB_DIR}")
    print(f"  OUTPUT_ROOT:                  {OUTPUT_ROOT}")

    print("\nFunctional ecology")
    print(f"  RUN_FUNGALTRAITS_ITS:         {RUN_FUNGALTRAITS_ITS}")
    print(f"  FUNGAL_TRAITS_DB:             {FUNGAL_TRAITS_DB}")
    print(f"  RUN_FAPROTAX_16S:             {RUN_FAPROTAX_16S}")
    print(f"  RUN_FOR_R_EXPORT:             {RUN_FOR_R_EXPORT}")
    print(f"  RUN_SECONDARY_SUBSET_EXPORTS: {RUN_SECONDARY_SUBSET_EXPORTS}")
    print(f"  FAPROTAX_DB_16S:              {FAPROTAX_DB_16S}")

    print("\nKrona export")
    print(f"  RUN_KRONA_EXPORT:             {RUN_KRONA_EXPORT}")
    print(f"  RUN_KRONA_PER_SAMPLE_EXPORT:  {RUN_KRONA_PER_SAMPLE_EXPORT}")
    print(f"  KRONA_MODE:                   {KRONA_MODE}")
    print(f"  KRONA_MIN_ABUNDANCE:          {KRONA_MIN_ABUNDANCE}")
    print(f"  KRONA_MAX_ROWS:               {KRONA_MAX_ROWS}")
    print(f"  KRONA_TOP_ROWS:               {KRONA_TOP_ROWS}")
    print(f"  KRONA_LIMIT_ACTION:           {KRONA_LIMIT_ACTION}")
    print(f"  KRONA_FILL_UNCLASSIFIED:      {KRONA_FILL_UNCLASSIFIED}")
    print(f"  KRONA_EXPORT_FIXED_COLUMNS:   {KRONA_EXPORT_FIXED_COLUMNS}")

    print("\nReInput export")
    print(f"  RUN_REINPUT_EXPORT:           {RUN_REINPUT_EXPORT}")

    print("\nMemory / performance")
    print(f"  RAM_SAFE_MODE:                {RAM_SAFE_MODE}")
    print(f"  RAM_SAFE_CHUNK_ROWS:          {RAM_SAFE_CHUNK_ROWS}")
    print(f"  RAM_SAFE_SITE_CHUNK:          {RAM_SAFE_SITE_CHUNK}")
    print(f"  RAM_SAFE_GROUP_ROWS:          {RAM_SAFE_GROUP_ROWS}")
    print(f"  RAM_SAFE_ABUNDANCE_DTYPE:     {RAM_SAFE_ABUNDANCE_DTYPE}")
    print(f"  RAM_SAFE_SUM_DTYPE:           {RAM_SAFE_SUM_DTYPE}")

    print("=" * 80)
    print()


def warn_about_configuration():
    """
    Print non-fatal configuration warnings.
    """
    warnings = []

    if not INPUT_DIR.exists():
        warnings.append(
            f"INPUT_DIR does not exist yet: {INPUT_DIR}"
        )

    if RUN_FUNGALTRAITS_ITS and MODE != "ITS":
        warnings.append(
            "RUN_FUNGALTRAITS_ITS is enabled, but MODE is not ITS. "
            "This is usually harmless if the downstream code only runs this "
            "module in ITS mode."
        )

    if RUN_FUNGALTRAITS_ITS and MODE == "ITS" and not FUNGAL_TRAITS_DB.exists():
        warnings.append(
            f"FungalTraits database was not found: {FUNGAL_TRAITS_DB}"
        )

    if RUN_FAPROTAX_16S and MODE != "16S":
        warnings.append(
            "RUN_FAPROTAX_16S is enabled, but MODE is not 16S. "
            "This is usually harmless if the downstream code only runs this "
            "module in 16S mode."
        )

    if RUN_FAPROTAX_16S and MODE == "16S" and not FAPROTAX_DB_16S.exists():
        warnings.append(
            f"FAPROTAX database was not found: {FAPROTAX_DB_16S}"
        )

    if KRONA_LIMIT_ACTION == "TOP" and KRONA_TOP_ROWS > KRONA_MAX_ROWS:
        warnings.append(
            "KRONA_TOP_ROWS is larger than KRONA_MAX_ROWS. "
            "This is allowed, but usually KRONA_TOP_ROWS should be <= "
            "KRONA_MAX_ROWS."
        )

    if warnings:
        print("Configuration warnings:")
        for warning in warnings:
            print(f"  - {warning}")
        print()


if not ARGS.quiet_config:
    print_configuration()
    warn_about_configuration()

# ============================================================================
# RAM-SAFE HELPERS PORTED FROM METADIV v1.7.13 NOTEBOOK
# ============================================================================

RAM_SAFE_MEMMAP_DIRNAME = "_ram_safe_memmap"
_RAM_SAFE_TEMP_FILES = []

def _ram_safe_memmap_dir() -> Path:
    path = OUTPUT_ROOT / MODE / RAM_SAFE_MEMMAP_DIRNAME
    path.mkdir(parents=True, exist_ok=True)
    return path

def _register_ram_safe_temp_file(path: Path):
    path = Path(path)
    if path not in _RAM_SAFE_TEMP_FILES:
        _RAM_SAFE_TEMP_FILES.append(path)

def _cleanup_ram_safe_temp_files():
    gc.collect()
    for path in list(_RAM_SAFE_TEMP_FILES):
        try:
            if Path(path).exists():
                Path(path).unlink()
        except Exception as e:
            print(f"   ⚠️ Could not remove RAM-safe temp file yet: {path} ({e})")
    try:
        temp_dir = OUTPUT_ROOT / MODE / RAM_SAFE_MEMMAP_DIRNAME
        if temp_dir.exists() and not any(temp_dir.iterdir()):
            temp_dir.rmdir()
    except Exception:
        pass

def _ram_safe_filter_wide_dataframe(source_df, keep_mask, site_cols, tag):
    keep_mask = np.asarray(keep_mask, dtype=bool)
    if keep_mask.size != len(source_df):
        raise ValueError("RAM-safe filter received a keep mask with the wrong length.")

    valid_sites = [c for c in site_cols if c in source_df.columns]
    metadata_cols = [c for c in source_df.columns if c not in valid_sites]
    keep_positions = np.flatnonzero(keep_mask).astype(np.int64, copy=False)
    retained_n = int(keep_positions.size)

    print(f"\n💾 RAM-SAFE ROW COMPACTION ({tag})...")
    print(f"   • Source rows: {len(source_df):,}")
    print(f"   • Retained rows: {retained_n:,}")
    print(f"   • Abundance columns: {len(valid_sites):,}")

    if metadata_cols:
        metadata_positions = source_df.columns.get_indexer(metadata_cols)
        metadata_df = source_df.iloc[keep_positions, metadata_positions].copy()
        metadata_df.reset_index(drop=True, inplace=True)
    else:
        metadata_df = pd.DataFrame(index=np.arange(retained_n))

    if not valid_sites:
        return metadata_df

    temp_dir = _ram_safe_memmap_dir()
    memmap_path = temp_dir / f"{tag}_{os.getpid()}_{int(time.time() * 1000)}.dat"
    _register_ram_safe_temp_file(memmap_path)

    abundance_map = np.memmap(
        memmap_path,
        dtype=RAM_SAFE_ABUNDANCE_DTYPE,
        mode="w+",
        shape=(retained_n, len(valid_sites))
    )

    site_positions = source_df.columns.get_indexer(valid_sites)

    for out_start in range(0, retained_n, RAM_SAFE_CHUNK_ROWS):
        out_stop = min(out_start + RAM_SAFE_CHUNK_ROWS, retained_n)
        src_positions = keep_positions[out_start:out_stop]
        block = source_df.iloc[src_positions, site_positions].to_numpy(
            dtype=RAM_SAFE_ABUNDANCE_DTYPE,
            copy=True
        )
        abundance_map[out_start:out_stop, :] = block
        del block
        if out_start == 0 or out_stop == retained_n or out_stop % 50000 == 0:
            print(f"   • Compacted rows: {out_stop:,}/{retained_n:,}")
        gc.collect()

    abundance_map.flush()

    # Remove the complete abundance block from the source before constructing
    # the returned wide DataFrame. This prevents two full dense matrices from
    # coexisting in physical RAM.
    source_df.drop(columns=valid_sites, inplace=True)
    gc.collect()

    abundance_df = pd.DataFrame(abundance_map, columns=valid_sites, copy=False)
    final_df = abundance_df
    for col in reversed(metadata_cols):
        final_df.insert(0, col, metadata_df[col].to_numpy(copy=False))

    del metadata_df
    del abundance_map
    gc.collect()
    print(f"   ✅ RAM-safe compaction completed: {retained_n:,} rows")
    return final_df

def _ram_safe_filter_wide_dataframe_from_matrix(
    source_df,
    abundance_matrix,
    keep_mask,
    site_cols,
    tag
):
    """
    RAM-safe row compaction using the original numeric abundance matrix directly.

    This v1.7.13 path avoids thousands of repeated pandas selections during
    taxonomic collapse. Metadata are retained from ``source_df`` while abundance
    rows are copied from the NumPy matrix into a disk-backed memmap.
    """
    keep_mask = np.asarray(keep_mask, dtype=bool)
    if keep_mask.size != len(source_df):
        raise ValueError("RAM-safe matrix filter received a keep mask with the wrong length.")

    valid_sites = [c for c in site_cols if c in source_df.columns]
    metadata_cols = [c for c in source_df.columns if c not in valid_sites]
    keep_positions = np.flatnonzero(keep_mask).astype(np.int64, copy=False)
    retained_n = int(keep_positions.size)

    if abundance_matrix is None:
        raise ValueError("The v1.7.13 matrix compaction path requires an abundance matrix.")
    if abundance_matrix.shape != (len(source_df), len(valid_sites)):
        raise ValueError(
            "Abundance matrix shape does not match the source table: "
            f"{abundance_matrix.shape} vs {(len(source_df), len(valid_sites))}."
        )

    print(f"\n💾 RAM-SAFE ROW COMPACTION ({tag}; direct NumPy matrix)...")
    print(f"   • Source rows: {len(source_df):,}")
    print(f"   • Retained rows: {retained_n:,}")
    print(f"   • Abundance columns: {len(valid_sites):,}")

    if metadata_cols:
        metadata_positions = source_df.columns.get_indexer(metadata_cols)
        metadata_df = source_df.iloc[keep_positions, metadata_positions].copy()
        metadata_df.reset_index(drop=True, inplace=True)
    else:
        metadata_df = pd.DataFrame(index=np.arange(retained_n))

    if not valid_sites:
        return metadata_df

    temp_dir = _ram_safe_memmap_dir()
    memmap_path = temp_dir / f"{tag}_{os.getpid()}_{int(time.time() * 1000)}.dat"
    _register_ram_safe_temp_file(memmap_path)

    abundance_map = np.memmap(
        memmap_path,
        dtype=RAM_SAFE_ABUNDANCE_DTYPE,
        mode="w+",
        shape=(retained_n, len(valid_sites))
    )

    for out_start in range(0, retained_n, RAM_SAFE_CHUNK_ROWS):
        out_stop = min(out_start + RAM_SAFE_CHUNK_ROWS, retained_n)
        src_positions = keep_positions[out_start:out_stop]
        abundance_map[out_start:out_stop, :] = abundance_matrix[src_positions, :]
        if out_start == 0 or out_stop == retained_n or out_stop % 50000 == 0:
            print(f"   • Compacted rows: {out_stop:,}/{retained_n:,}")

    abundance_map.flush()

    # Release the wide numeric block held by the source DataFrame before the
    # returned memmap-backed DataFrame is assembled.
    source_df.drop(columns=valid_sites, inplace=True)
    gc.collect()

    abundance_df = pd.DataFrame(abundance_map, columns=valid_sites, copy=False)
    final_df = abundance_df
    for col in reversed(metadata_cols):
        final_df.insert(0, col, metadata_df[col].to_numpy(copy=False))

    del metadata_df
    del abundance_map
    gc.collect()
    print(f"   ✅ RAM-safe compaction completed: {retained_n:,} rows")
    return final_df


def _abundance_order_positions(df):
    if "Total_Abundance" not in df.columns:
        return np.arange(len(df), dtype=np.int64)
    return (
        df["Total_Abundance"]
        .sort_values(ascending=False, kind="stable")
        .index
        .to_numpy(dtype=np.int64, copy=False)
    )

def _write_dataframe_in_order_chunks(df, row_positions, output_path, columns, chunk_rows=None):
    chunk_rows = int(chunk_rows or RAM_SAFE_CHUNK_ROWS)
    row_positions = np.asarray(row_positions, dtype=np.int64)
    columns = [c for c in columns if c in df.columns]
    col_positions = df.columns.get_indexer(columns)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    first = True
    total_written = 0
    for start in range(0, len(row_positions), chunk_rows):
        stop = min(start + chunk_rows, len(row_positions))
        positions = row_positions[start:stop]
        chunk = df.iloc[positions, col_positions]
        chunk.to_csv(
            output_path,
            mode="w" if first else "a",
            header=first,
            index=False,
            encoding="utf-8"
        )
        total_written += len(chunk)
        first = False
        del chunk
        gc.collect()
    return total_written

# ============================================================================
# DERIVED CONFIGURATION
# ============================================================================

OUTPUT_DIR = OUTPUT_ROOT / MODE
FINAL_DB_DIR = OUTPUT_DIR / "FINAL_DB"
FOR_R_DIR = OUTPUT_DIR / "For_R"
FUNCTIONAL_ECOLOGY_DIR = OUTPUT_DIR / "Functional_Ecology"
KRONA_DIR = OUTPUT_DIR / "Krona"
REINPUT_DIR = FINAL_DB_DIR / "ReInput"
CONCAT_DIR = OUTPUT_DIR / "concatenated_tables"  # Always created, deleted later if not DEV_MODE

# Create output folders for the active marker only.
# This prevents runs in one marker mode from modifying output folders from
# another marker mode.
for d in [OUTPUT_ROOT, OUTPUT_DIR, FINAL_DB_DIR, FOR_R_DIR, REINPUT_DIR, CONCAT_DIR]:
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
    # FUNGAL_TRAITS_DB is resolved from CLI settings; preserve custom path/filename.
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
    """
    Build a run-specific FINAL_DB filename containing both confidence thresholds.

    Example
    -------
    Final_Database_genus_all_eukaryotes_p1_sppn08.csv
    """
    sintax_tag = format_threshold_for_name(P_VALUE_THRESHOLD)
    sppn_tag = f"sppn{format_threshold_short(SPPN_P_THRESHOLD)}"
    return f"Final_Database_{COLLAPSE_STRATEGY}_{suffix}_{sintax_tag}_{sppn_tag}.csv"

def build_collapse_log_filename(suffix: str) -> str:
    """
    Build a run-specific collapse-log filename containing both confidence thresholds.

    Example
    -------
    Collapse_Log_species_only_all_eukaryotes_p1_sppn08.txt
    """
    sintax_tag = format_threshold_for_name(P_VALUE_THRESHOLD)
    sppn_tag = f"sppn{format_threshold_short(SPPN_P_THRESHOLD)}"
    return f"Collapse_Log_{COLLAPSE_STRATEGY}_{suffix}_{sintax_tag}_{sppn_tag}.txt"

def build_report_filename(suffix: str) -> str:
    """
    Build a run-specific report filename containing both confidence thresholds.

    This prevents reports generated with the same SINTAX threshold but different
    SPPN parsing thresholds from overwriting one another.

    Examples
    --------
    METADIV_REPORT_genus_all_eukaryotes_p1_sppn08.txt
    METADIV_REPORT_genus_all_eukaryotes_p08_sppn05.txt
    """
    sintax_tag = format_threshold_for_name(P_VALUE_THRESHOLD)
    sppn_tag = f"sppn{format_threshold_short(SPPN_P_THRESHOLD)}"
    return (
        f"METADIV_REPORT_{COLLAPSE_STRATEGY}_{suffix}_"
        f"{sintax_tag}_{sppn_tag}.txt"
    )

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

    # Output names include collapse strategy
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()
    functional_output_name = (
        f"Final_Database_ITS_FungalTraits_full_annotated_{collapse_label}.csv"
    )
    functional_log_name = (
        f"FungalTraits_Module_Log_{collapse_label}.txt"
    )

    if traits_df is None or traits_df.empty or "genus" not in final_df.columns:
        annotated = prepare_functional_ecology_export_table(final_df)
        annotated_path = functional_dir / functional_output_name
        annotated.to_csv(annotated_path, index=False, encoding="utf-8")
        log_text = (
            "FungalTraits ITS Functional Ecology module\n"
            "Status: no complete annotation was added.\n"
            "Reason: FungalTraits database is missing/empty or final database has no genus column.\n"
        )
        (functional_dir / functional_log_name).write_text(log_text, encoding="utf-8")
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

    annotated_path = functional_dir / functional_output_name
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
    (functional_dir / functional_log_name).write_text(log_text, encoding="utf-8")

    print(f"   🍄 ITS Functional_Ecology full FungalTraits table saved: {annotated_path}")
    print(f"   • FungalTraits matches: {matched:,}/{total:,} rows ({(matched / total * 100) if total else 0:.1f}%)")
    print(f"   • FungalTraits columns exported: {len(ft_export_cols):,}")
    return final_df


# ============================================================================
# 16S FAPROTAX FUNCTIONAL ECOLOGY MODULE
# ============================================================================

FAPROTAX_RANK_COLUMNS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

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


def match_faprotax_compact(
    taxonomy_df: pd.DataFrame,
    rules: list
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Match MetaDiv 16S taxa against FAPROTAX without modifying FINAL_DB.

    v1.7.13 keeps only the information required to build the functional
    abundance table: feature ID -> FAPROTAX function. Feature-level annotation
    strings, matched-taxon columns, match-rank columns and source columns are
    intentionally not constructed.

    Matching semantics are unchanged: candidate rules are selected by terminal
    taxon and accepted only when their tokens form an ordered subset of the
    feature taxonomy.
    """
    total = len(taxonomy_df)
    assigned_mask = np.zeros(total, dtype=bool)

    id_col = "SPPN" if "SPPN" in taxonomy_df.columns else "OTU_XX"
    compact_cols = [id_col] + [
        rank for rank in FAPROTAX_RANK_COLUMNS
        if rank in taxonomy_df.columns
    ]
    compact = taxonomy_df.loc[:, compact_cols]

    if not rules or compact.empty:
        return pd.DataFrame(columns=[id_col, "function"]), assigned_mask

    rules_by_terminal = defaultdict(list)
    for rule in rules:
        rules_by_terminal[rule["terminal"]].append(rule)

    rank_names = compact_cols[1:]
    annotation_cache = {}
    long_records = []

    print("\n🔎 MATCHING MetaDiv TAXA AGAINST FAPROTAX RULES (v1.7.13 COMPACT MODE)...")
    print(f"   • Taxonomic columns scanned: {len(rank_names)}")
    print("   • Abundance columns are NOT scanned during rule matching.")
    print("   • Feature-level FAPROTAX annotation columns are NOT constructed.")

    for row_position, values in enumerate(
        compact.itertuples(index=False, name=None)
    ):
        row_id = values[0]
        raw_tax_values = values[1:]

        normalized_rank_tokens = tuple(
            _normalize_faprotax_token(value)
            for value in raw_tax_values
        )

        matched_functions = annotation_cache.get(normalized_rank_tokens)

        if matched_functions is None:
            row_tokens = []
            for rank, token in zip(rank_names, normalized_rank_tokens):
                if token:
                    row_tokens.append(token)

                # FAPROTAX rules can encode the species epithet separately
                # from the full MetaDiv binomial.
                if rank == "species" and token and "_" in token:
                    epithet = token.split("_")[-1]
                    if epithet:
                        row_tokens.append(epithet)

            candidate_rules = []
            for token in set(row_tokens):
                candidate_rules.extend(rules_by_terminal.get(token, []))

            functions = set()
            for rule in candidate_rules:
                if _tokens_are_ordered_subset(rule["tokens"], row_tokens):
                    functions.add(rule["function"])

            matched_functions = tuple(sorted(functions))
            annotation_cache[normalized_rank_tokens] = matched_functions

        if matched_functions:
            assigned_mask[row_position] = True
            for function in matched_functions:
                long_records.append({
                    id_col: row_id,
                    "function": function,
                })

        row_number = row_position + 1
        if row_number % 100000 == 0:
            print(
                f"   • FAPROTAX matched {row_number:,}/{total:,} features "
                f"| cached taxonomies: {len(annotation_cache):,}"
            )

    long_df = pd.DataFrame(long_records, columns=[id_col, "function"])

    del compact, annotation_cache, long_records
    gc.collect()
    return long_df, assigned_mask


def build_faprotax_functional_abundance(
    df: pd.DataFrame,
    long_df: pd.DataFrame,
    site_cols: list
) -> pd.DataFrame:
    """
    Build a FAPROTAX function-by-sample abundance table without the large
    long_df x abundance merge used by v1.7.12.

    Only small blocks of sample columns are materialized at a time.
    SciPy is not used.
    """
    valid_sites = [c for c in site_cols if c in df.columns]
    if long_df.empty or not valid_sites:
        return pd.DataFrame(columns=["function"] + valid_sites)

    print("\n🧠 Building FAPROTAX functional abundance table (RAM-SAFE)...")

    id_col = "SPPN" if "SPPN" in df.columns else "OTU_XX"

    # Convert feature IDs to row positions.
    id_to_position = pd.Series(
        np.arange(len(df), dtype=np.int64),
        index=df[id_col].fillna("").astype(str)
    )

    work = long_df.loc[:, [id_col, "function"]].copy()
    work["_row_position"] = (
        work[id_col].fillna("").astype(str).map(id_to_position)
    )
    work = work.dropna(subset=["_row_position"])
    work["_row_position"] = work["_row_position"].astype(np.int64)

    grouped_positions = {
        function: group["_row_position"].to_numpy(dtype=np.int64, copy=True)
        for function, group in work.groupby("function", sort=True)
    }
    functions = list(grouped_positions.keys())

    result_matrix = np.zeros(
        (len(functions), len(valid_sites)),
        dtype=RAM_SAFE_SUM_DTYPE
    )

    for start in range(0, len(valid_sites), RAM_SAFE_SITE_CHUNK):
        stop = min(start + RAM_SAFE_SITE_CHUNK, len(valid_sites))
        chunk_sites = valid_sites[start:stop]

        block = df.loc[:, chunk_sites].to_numpy(
            dtype=RAM_SAFE_SUM_DTYPE,
            copy=True
        )

        for function_index, function in enumerate(functions):
            positions = grouped_positions[function]
            if positions.size:
                result_matrix[function_index, start:stop] = (
                    block[positions, :]
                    .sum(axis=0, dtype=RAM_SAFE_SUM_DTYPE)
                )

        del block
        gc.collect()

    profile = pd.DataFrame(result_matrix, columns=valid_sites)
    profile.insert(0, "function", functions)

    del result_matrix, grouped_positions, work, id_to_position
    gc.collect()
    print("   ✅ FAPROTAX functional abundance table completed.")
    return profile

def add_faprotax_functional_ecology_16s(
    final_df: pd.DataFrame,
    site_cols: list,
    output_dir: Path,
    faprotax_path: Path
) -> pd.DataFrame:
    """
    Run the compact MetaDiv Builder v1.7.13 16S FAPROTAX module.

    The functional module uses a narrow taxonomy-only table for FAPROTAX
    matching and writes only the function-by-sample abundance profile plus a
    concise audit log. Two compact feature-level annotations are retained in
    FINAL_DB because they are directly useful for downstream ecology:

        - faprotax_functions
        - faprotax_function_count

    Diagnostic fields such as matched taxa, match ranks and source are not
    created or exported.

    Outputs
    -------
    output/16S/Functional_Ecology/
        - FAPROTAX_Functional_Abundance_Table.csv
        - FAPROTAX_Module_Log.txt
    """
    if MODE != "16S" or not RUN_FAPROTAX_16S:
        return final_df

    functional_dir = output_dir / "Functional_Ecology"
    functional_dir.mkdir(parents=True, exist_ok=True)

    # Remove legacy v1.7.12/v4 files from the same output directory so a rerun
    # cannot be mistaken for having generated them under v1.7.13.
    legacy_names = [
        "Final_Database_16S_FAPROTAX_annotated.csv",
        "FAPROTAX_Long_Format_Matches.csv",
        "FAPROTAX_Unassigned_Taxa.csv",
        "FAPROTAX_Function_Metadata.csv",
    ]
    for legacy_name in legacy_names:
        legacy_path = functional_dir / legacy_name
        if legacy_path.exists():
            try:
                legacy_path.unlink()
            except Exception as exc:
                print(f"   ⚠️ Could not remove legacy FAPROTAX output {legacy_path.name}: {exc}")

    rules, metadata_by_function = parse_faprotax_database(faprotax_path)

    id_col = "SPPN" if "SPPN" in final_df.columns else "OTU_XX"
    annotation_cols = [id_col] + [
        rank for rank in FAPROTAX_RANK_COLUMNS
        if rank in final_df.columns
    ]
    annotation_df = final_df.loc[:, annotation_cols].copy()

    long_df, assigned_mask = match_faprotax_compact(
        annotation_df,
        rules
    )

    total = len(annotation_df)
    assigned = int(assigned_mask.sum())

    # ------------------------------------------------------------------
    # Compact feature-level annotation for FINAL_DB.
    # Reuse long_df already created for functional abundance aggregation;
    # no second FAPROTAX matching pass and no abundance recalculation.
    # ------------------------------------------------------------------
    if long_df.empty:
        final_df["faprotax_functions"] = ""
        final_df["faprotax_function_count"] = np.zeros(
            len(final_df), dtype=np.uint16
        )
    else:
        summary_work = long_df.loc[:, [id_col, "function"]].copy()
        summary_work[id_col] = summary_work[id_col].fillna("").astype(str)
        function_summary = (
            summary_work.groupby(id_col, sort=False)["function"]
            .agg(lambda values: ";".join(sorted(set(map(str, values)))))
        )
        function_counts = (
            summary_work.groupby(id_col, sort=False)["function"]
            .nunique()
        )

        feature_ids = final_df[id_col].fillna("").astype(str)
        final_df["faprotax_functions"] = (
            feature_ids.map(function_summary).fillna("")
        )
        final_df["faprotax_function_count"] = (
            feature_ids.map(function_counts).fillna(0).astype(np.uint16)
        )
        del summary_work, function_summary, function_counts, feature_ids
        gc.collect()

    # The only CSV generated by the 16S functional module.
    profile = build_faprotax_functional_abundance(
        final_df,
        long_df,
        site_cols
    )
    profile_path = functional_dir / "FAPROTAX_Functional_Abundance_Table.csv"
    profile.to_csv(profile_path, index=False, encoding="utf-8")

    # One concise audit log documents assignment and output statistics.
    log_path = functional_dir / "FAPROTAX_Module_Log.txt"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("MetaDiv Builder v1.7.13 - 16S FAPROTAX Functional Ecology Module\n")
        f.write("Author: Bernardo Águila, UNAM\n")
        f.write("=" * 70 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"FAPROTAX database used: {faprotax_path}\n")
        f.write(f"Parsed FAPROTAX functions: {len(metadata_by_function):,}\n")
        f.write(f"Parsed explicit matching rules: {len(rules):,}\n")
        f.write(f"Total MetaDiv features evaluated: {total:,}\n")
        f.write(f"Assigned features: {assigned:,}\n")
        f.write(f"Unassigned features: {total - assigned:,}\n")
        f.write(
            f"Assignment percentage: "
            f"{(assigned / total * 100) if total else 0:.2f}%\n"
        )
        f.write(f"Feature-function matches used for aggregation: {len(long_df):,}\n")
        f.write(f"Functional groups exported: {len(profile):,}\n")
        f.write("\nOutput files:\n")
        f.write(f"  {profile_path.name}\n")
        f.write(f"  {log_path.name}\n")
        f.write("\nFINAL_DB policy:\n")
        f.write("  Included: faprotax_functions, faprotax_function_count\n")
        f.write(
            "  Excluded: faprotax_matched_taxa, faprotax_match_ranks, "
            "faprotax_source and other diagnostic feature-level fields.\n"
        )
        f.write("\nMost frequent FAPROTAX functions among matched features:\n")
        if not long_df.empty and "function" in long_df.columns:
            counts = long_df["function"].value_counts().head(50)
            for function, count in counts.items():
                f.write(f"  {function}: {int(count):,}\n")

    print("\n🌱 16S FAPROTAX MODULE COMPLETED — COMPACT v1.7.13")
    print(
        f"   • Assigned features: {assigned:,}/{total:,} "
        f"({(assigned / total * 100) if total else 0:.2f}%)"
    )
    print(f"   • Functional abundance table: {profile_path}")
    print(f"   • Module log: {log_path}")
    print("   • FINAL_DB annotations: faprotax_functions + faprotax_function_count")
    print("   • Diagnostic FAPROTAX columns: excluded from FINAL_DB")
    print("   • Redundant FAPROTAX CSV outputs: disabled by design")

    del annotation_df, long_df, profile, assigned_mask
    gc.collect()

    return final_df


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

    # Output names include collapse strategy
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()
    functional_output_name = (
        f"Final_Database_ITS_FungalTraits_full_annotated_{collapse_label}.csv"
    )
    functional_log_name = (
        f"FungalTraits_Module_Log_{collapse_label}.txt"
    )

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


def process_single_dataset(prefix: str, files: dict, report_lines: list) -> dict:
    """
    RAM-safe Part 1 processor.

    The original v1.7.12 loaded the full feature x sample abundance table as
    Python strings, which can require several gigabytes for wide datasets.
    This implementation reads abundance tables in row chunks, joins taxonomy
    and sequence information chunk-by-chunk, converts abundances to uint32,
    and writes intermediate files incrementally.

    Biological behavior is unchanged:
    - OTU_XX identifiers follow original abundance-table row order.
    - taxonomy and representative sequences are matched by original feature ID.
    - no abundance filtering or normalization is performed.
    - target-subset filtering follows the same MODE / SUBSET_MODE logic.
    """
    print(f"\n=== PROCESSING DATASET: {prefix} (RAM-SAFE) ===")

    # ------------------------------------------------------------------
    # Read only the header first.
    # ------------------------------------------------------------------
    header = pd.read_csv(files["abundance"], nrows=0)
    if len(header.columns) < 2:
        raise ValueError(f"Abundance table for {prefix} must contain an ID column and at least one sample column.")

    id_col = header.columns[0]
    site_columns = list(header.columns[1:])
    print(f"  ID column detected: '{id_col}'")
    print(f"  Sample columns detected: {len(site_columns):,}")

    # ------------------------------------------------------------------
    # Taxonomy and FASTA are much smaller than the wide abundance matrix.
    # Keep only the mappings needed for chunk-wise lookup.
    # ------------------------------------------------------------------
    print("  Reading SINTAX taxonomy...")
    taxonomy_df = read_taxonomy_sintax(files["taxonomy"])
    if taxonomy_df.empty:
        taxonomy_map = pd.Series(dtype="object")
        taxonomy_ids = set()
    else:
        taxonomy_map = taxonomy_df["sintax_taxonomy"].fillna("")
        taxonomy_ids = set(taxonomy_df.index.astype(str))

    print("  Reading FASTA sequences...")
    sequences_dict = read_fasta_sequences(files["sequences"])
    seq_ids = set(str(k) for k in sequences_dict.keys())

    # ------------------------------------------------------------------
    # Resolve output paths using the same names as v1.7.12.
    # ------------------------------------------------------------------
    concat_path = CONCAT_DIR / f"{prefix}_concatenated.csv"

    if MODE == "ITS":
        target_name = "fungi" if SUBSET_MODE == "only_fungi" else "all_eukaryotes"
        target_dir = FUNGI_DIR if SUBSET_MODE == "only_fungi" else ALL_EUK_DIR
        all_name = "all_eukaryotes"
        all_dir = ALL_EUK_DIR
    elif MODE == "16S":
        target_name = "only_bacteria" if SUBSET_MODE == "only_bacteria" else "all_prokaryotes"
        target_dir = BACTERIA_DIR if SUBSET_MODE == "only_bacteria" else ALL_PROK_DIR
        all_name = "all_prokaryotes"
        all_dir = ALL_PROK_DIR
    elif MODE == "CO1":
        target_name = "only_metazoa" if SUBSET_MODE == "only_metazoa" else "all_eukaryotes"
        target_dir = METAZOA_DIR if SUBSET_MODE == "only_metazoa" else CO1_DIR
        all_name = "all_eukaryotes"
        all_dir = CO1_DIR
    else:
        raise ValueError("MODE must be 'ITS', '16S', or 'CO1'.")

    target_path = target_dir / f"{prefix}_{target_name}.csv"
    all_path = all_dir / f"{prefix}_{all_name}.csv"
    other_path = NON_ANNOTATED_DIR / f"{prefix}_non_annotated.csv"

    # Remove stale intermediate outputs BEFORE processing. If Part 1 fails,
    # Part 2 can no longer silently reuse a file from an older run.
    paths_to_reset = {concat_path.resolve(), target_path.resolve(), all_path.resolve(), other_path.resolve()}
    for p in paths_to_reset:
        p = Path(p)
        if p.exists():
            p.unlink()

    total_ids = 0
    target_ids = 0
    other_ids = 0
    bacteria_ids = 0
    archaea_ids = 0
    fungi_ids = 0
    eukaryota_ids = 0
    metazoa_ids = 0
    missing_tax_count = 0
    missing_seq_count = 0

    first_concat = True
    first_target = True
    first_all = True
    first_other = True
    row_offset = 0

    print(f"  Reading abundance in chunks of {RAM_SAFE_CHUNK_ROWS:,} rows...")

    reader = pd.read_csv(
        files["abundance"],
        chunksize=RAM_SAFE_CHUNK_ROWS,
        dtype={id_col: str},
        low_memory=False
    )

    for chunk_number, raw_chunk in enumerate(reader, start=1):
        raw_chunk[id_col] = raw_chunk[id_col].astype(str)
        ids = raw_chunk[id_col]
        n_chunk = len(raw_chunk)
        total_ids += n_chunk

        # Map taxonomy and sequences without joining the full wide table.
        taxonomy_values = ids.map(taxonomy_map).fillna("")
        sequence_values = ids.map(sequences_dict).fillna("")

        if taxonomy_ids:
            missing_tax_count += int((~ids.isin(taxonomy_ids)).sum())
        else:
            missing_tax_count += n_chunk
        missing_seq_count += int((~ids.isin(seq_ids)).sum())

        # Convert only the current abundance block.
        abundance_numeric = (
            raw_chunk.loc[:, site_columns]
            .apply(pd.to_numeric, errors="coerce")
            .fillna(0)
        )

        values = abundance_numeric.to_numpy(copy=False)
        if values.size:
            min_value = np.nanmin(values)
            max_value = np.nanmax(values)
            if min_value < 0:
                raise ValueError(
                    f"Negative abundance detected in {prefix}. MetaDiv expects non-negative read counts."
                )
            if max_value > np.iinfo(RAM_SAFE_ABUNDANCE_DTYPE).max:
                raise OverflowError(
                    f"Abundance {max_value} exceeds {RAM_SAFE_ABUNDANCE_DTYPE}. "
                    "Change RAM_SAFE_ABUNDANCE_DTYPE to np.uint64."
                )

        abundance_numeric = abundance_numeric.astype(RAM_SAFE_ABUNDANCE_DTYPE, copy=False)

        # Build only the current output chunk. pd.concat is used once here
        # instead of inserting hundreds/thousands of sample columns one-by-one.
        metadata_chunk = pd.DataFrame({
            "OTU_XX": [
                f"{prefix}_ID{str(i + 1).zfill(5)}"
                for i in range(row_offset, row_offset + n_chunk)
            ],
            "Original_ID": ids.to_numpy(copy=False),
            "sintax_taxonomy": taxonomy_values.to_numpy(copy=False),
            "sequence": sequence_values.to_numpy(copy=False),
        })
        abundance_numeric.reset_index(drop=True, inplace=True)
        out = pd.concat(
            [metadata_chunk.reset_index(drop=True), abundance_numeric],
            axis=1,
            copy=False
        )

        row_offset += n_chunk

        tax_lower = out["sintax_taxonomy"].astype(str).str.lower()

        if MODE == "ITS":
            fungi_mask = (
                tax_lower.str.contains("d:fungi", na=False, case=False) |
                tax_lower.str.contains("k:fungi", na=False, case=False)
            )
            fungi_ids += int(fungi_mask.sum())
            if SUBSET_MODE == "only_fungi":
                target_mask = fungi_mask
            else:
                # Preserve v1.7.12 behavior: all_eukaryotes mode is the broad
                # integrated table and is not additionally filtered here.
                target_mask = pd.Series(True, index=out.index)

        elif MODE == "16S":
            bacteria_mask = (
                tax_lower.str.contains("d:bacteria", na=False, case=False) |
                tax_lower.str.contains("k:bacteria", na=False, case=False)
            )
            archaea_mask = (
                tax_lower.str.contains("d:archaea", na=False, case=False) |
                tax_lower.str.contains("k:archaea", na=False, case=False)
            )
            bacteria_ids += int(bacteria_mask.sum())
            archaea_ids += int(archaea_mask.sum())
            target_mask = bacteria_mask if SUBSET_MODE == "only_bacteria" else (bacteria_mask | archaea_mask)

        elif MODE == "CO1":
            eukaryota_mask = (
                tax_lower.str.contains("d:eukaryota", na=False, case=False) |
                tax_lower.str.contains("k:eukaryota", na=False, case=False)
            )
            metazoa_mask = build_co1_metazoa_mask(tax_lower)
            eukaryota_ids += int(eukaryota_mask.sum())
            metazoa_ids += int(metazoa_mask.sum())
            target_mask = metazoa_mask if SUBSET_MODE == "only_metazoa" else eukaryota_mask

        target_chunk = out.loc[target_mask]
        other_chunk = out.loc[~target_mask]
        target_ids += len(target_chunk)
        other_ids += len(other_chunk)

        # Concatenated full table.
        out.to_csv(
            concat_path,
            mode="w" if first_concat else "a",
            header=first_concat,
            index=False
        )
        first_concat = False

        # Broad marker table. If this is the same path as target_path,
        # write only the target table to avoid duplicate writes.
        if all_path.resolve() != target_path.resolve():
            out.to_csv(
                all_path,
                mode="w" if first_all else "a",
                header=first_all,
                index=False
            )
            first_all = False

        # Preserve v1.7.12 broad-subset behavior. When target_path and
        # all_path are the same file (e.g. all_prokaryotes), the original code
        # ultimately wrote the complete concatenated table to that path.
        target_write_chunk = (
            out
            if target_path.resolve() == all_path.resolve()
            else target_chunk
        )
        target_write_chunk.to_csv(
            target_path,
            mode="w" if first_target else "a",
            header=first_target,
            index=False
        )
        first_target = False

        if len(other_chunk):
            other_chunk.to_csv(
                other_path,
                mode="w" if first_other else "a",
                header=first_other,
                index=False
            )
            first_other = False

        if chunk_number == 1 or chunk_number % 10 == 0:
            print(f"    • chunk {chunk_number:,}: {total_ids:,} rows processed")

        del raw_chunk, abundance_numeric, metadata_chunk, out, target_chunk, target_write_chunk, other_chunk
        gc.collect()

    if missing_tax_count:
        warning = f"      ⚠️ WARNING: {missing_tax_count} IDs in abundance missing from taxonomy file"
        print(warning)
        report_lines.append(f"  {prefix}: {warning}")
    if missing_seq_count:
        warning = f"      ⚠️ WARNING: {missing_seq_count} IDs in abundance missing from FASTA file"
        print(warning)
        report_lines.append(f"  {prefix}: {warning}")

    stats = {
        "prefix": prefix,
        "total_ids": total_ids,
        "target_ids": target_ids,
        "other_ids": other_ids,
        "target_percentage": (target_ids / total_ids * 100) if total_ids else 0.0,
    }

    if MODE == "ITS":
        stats["fungi_ids"] = fungi_ids
        detected = fungi_ids
        filter_desc = f"{SUBSET_MODE} ({fungi_ids} Fungi detected)"
    elif MODE == "16S":
        stats["bacteria_ids"] = bacteria_ids
        stats["archaea_ids"] = archaea_ids
        detected = bacteria_ids if SUBSET_MODE == "only_bacteria" else bacteria_ids + archaea_ids
        filter_desc = f"{SUBSET_MODE} ({bacteria_ids} Bacteria, {archaea_ids} Archaea detected)"
    else:
        stats["eukaryota_ids"] = eukaryota_ids
        stats["metazoa_ids"] = metazoa_ids
        detected = metazoa_ids if SUBSET_MODE == "only_metazoa" else eukaryota_ids
        filter_desc = f"{SUBSET_MODE} ({eukaryota_ids} Eukaryota, {metazoa_ids} Metazoa detected)"

    stats["detected_percentage"] = (detected / total_ids * 100) if total_ids else 0.0

    print(f"✅ Concatenated saved incrementally: {concat_path}")
    print(
        f"✅ {prefix}: {target_ids} {filter_desc} "
        f"({stats['detected_percentage']:.1f}%)"
    )
    report_lines.append(
        f"  {prefix}: {target_ids} {filter_desc} "
        f"({stats['detected_percentage']:.1f}%)"
    )

    del taxonomy_df, taxonomy_map, sequences_dict
    gc.collect()
    return stats

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
        stats["detected_percentage"] = (
            stats["fungi_ids"] / stats["total_ids"] * 100
            if stats["total_ids"]
            else 0.0
        )
        filter_desc = (
            f"{SUBSET_MODE} "
            f"({stats['fungi_ids']} Fungi detected)"
        )

    elif MODE == "16S":
        stats["bacteria_ids"] = int(bacteria_mask.sum())
        stats["archaea_ids"] = int(archaea_mask.sum())
        detected_ids = (
            stats["bacteria_ids"]
            if SUBSET_MODE == "only_bacteria"
            else stats["bacteria_ids"] + stats["archaea_ids"]
        )
        stats["detected_percentage"] = (
            detected_ids / stats["total_ids"] * 100
            if stats["total_ids"]
            else 0.0
        )
        filter_desc = (
            f"{SUBSET_MODE} ({stats['bacteria_ids']} Bacteria, "
            f"{stats['archaea_ids']} Archaea detected)"
        )

    elif MODE == "CO1":
        stats["eukaryota_ids"] = int(eukaryota_mask.sum())
        stats["metazoa_ids"] = int(metazoa_mask.sum())
        detected_ids = (
            stats["metazoa_ids"]
            if SUBSET_MODE == "only_metazoa"
            else stats["eukaryota_ids"]
        )
        stats["detected_percentage"] = (
            detected_ids / stats["total_ids"] * 100
            if stats["total_ids"]
            else 0.0
        )
        filter_desc = (
            f"{SUBSET_MODE} ({stats['eukaryota_ids']} Eukaryota, "
            f"{stats['metazoa_ids']} Metazoa detected)"
        )

    else:
        stats["detected_percentage"] = stats["target_percentage"]
        filter_desc = SUBSET_MODE

    print(
        f"✅ {prefix}: {stats['target_ids']} {filter_desc} "
        f"({stats['detected_percentage']:.1f}%)"
    )
    report_lines.append(
        f"  {prefix}: {stats['target_ids']} {filter_desc} "
        f"({stats['detected_percentage']:.1f}%)"
    )

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

def create_binomial_species_inplace(df, mode=None):
    """
    Normalize species labels before SPPN assignment.

    ITS and 16S may provide only the specific epithet in the parsed ``species``
    column (for example ``coli``). For these markers, MetaDiv Builder creates a
    binomial label using the parsed genus whenever both ranks are available:

        Escherichia + coli -> Escherichia_coli

    Species labels that already contain the complete genus are preserved. CO1
    labels are left unchanged because CO1 reference databases commonly provide
    complete species names and identifier-rich taxonomic strings that should not
    be reconstructed.

    The elementwise prefix test intentionally avoids pandas ``str.startswith``
    with a Series pattern, which is not supported by pandas.
    """
    active_mode = str(mode or "").upper()

    if active_mode == "CO1":
        print("   • CO1 mode: existing species labels preserved; binomial reconstruction skipped.")
        return

    if active_mode not in {"ITS", "16S"}:
        print(f"   • {active_mode or 'UNKNOWN'} mode: binomial reconstruction skipped.")
        return

    if "genus" not in df.columns or "species" not in df.columns:
        print("   • Binomial reconstruction skipped: genus/species column missing.")
        return

    genus = df["genus"].fillna("").astype(str).str.strip()
    species = df["species"].fillna("").astype(str).str.strip()

    invalid_tokens = {"", "nan", "None", "none", "NA", "na", "N/A", "n/a"}
    genus_valid = ~genus.isin(invalid_tokens)
    species_valid = ~species.isin(invalid_tokens)

    genus_values = genus.str.lower().to_numpy(dtype=object, copy=False)
    species_values = species.str.lower().to_numpy(dtype=object, copy=False)

    # Per-row comparison is required because each species has a different genus.
    # np.fromiter is memory-light and avoids DataFrame.iterrows().
    already_binomial_array = np.fromiter(
        (
            bool(g) and bool(s) and (
                s == g or s.startswith(g + "_") or s.startswith(g + " ")
            )
            for g, s in zip(genus_values, species_values)
        ),
        dtype=bool,
        count=len(df)
    )
    already_binomial = pd.Series(already_binomial_array, index=df.index)
    already_binomial = already_binomial & genus_valid & species_valid

    build_mask = genus_valid & species_valid & ~already_binomial

    if build_mask.any():
        df.loc[build_mask, "species"] = (
            genus.loc[build_mask] + "_" + species.loc[build_mask]
        )

    created = int(build_mask.sum())
    preserved = int(already_binomial.sum())
    no_genus = int((species_valid & ~genus_valid).sum())
    no_species = int((genus_valid & ~species_valid).sum())

    print(f"   • {active_mode} binomials created: {created:,}")
    print(f"   • Existing binomials preserved: {preserved:,}")
    if no_genus or no_species:
        print(
            f"   • Warning: no genus for species: {no_genus:,} | "
            f"no species for genus: {no_species:,}"
        )

    del genus, species, genus_valid, species_valid
    del genus_values, species_values, already_binomial_array
    del already_binomial, build_mask
    gc.collect()

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
    """
    Assign SPPNs in descending Total_Abundance order without physically sorting
    the complete feature x sample DataFrame.
    """
    print(f"\n🏷️ ASSIGNING FINAL SPPN CODES (PER TAXON, p≥{SPPN_P_THRESHOLD:.2f})...")

    if "SPPN" in final_df.columns:
        final_df.drop(columns=["SPPN"], inplace=True)

    rank_order = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    abundance_order = _abundance_order_positions(final_df)
    sppn = np.empty(len(final_df), dtype=object)
    counter = {}
    col_pos = {c: final_df.columns.get_loc(c) for c in final_df.columns}

    for n, pos in enumerate(abundance_order, start=1):
        base_taxon = None
        for rank in rank_order:
            if rank not in col_pos:
                continue

            raw_name = final_df.iat[pos, col_pos[rank]]
            name = "" if pd.isna(raw_name) else str(raw_name).strip()

            pcol = f"{rank}_pvalue"
            pval = final_df.iat[pos, col_pos[pcol]] if pcol in col_pos else 0.0
            try:
                pval = float(pval)
            except Exception:
                pval = 0.0

            if name and pval >= SPPN_P_THRESHOLD and name not in {"nan", "None"}:
                base_taxon = name
                break

        if not base_taxon:
            base_taxon = "Unknown"

        base = re.sub(r"[^\w]+", "_", base_taxon)
        counter[base] = counter.get(base, 0) + 1
        sppn[pos] = f"{base}_{counter[base]:04d}"

        if n % 100000 == 0:
            print(f"   • SPPNs assigned: {n:,}/{len(final_df):,}")

    final_df["SPPN"] = sppn
    del sppn, abundance_order
    gc.collect()
    return final_df

def prune_zero_abundance(final_df, site_cols):
    print("\n🧹 PRUNING ZERO-ABUNDANCE IDs...")

    if "Total_Abundance" not in final_df.columns:
        totals = np.zeros(len(final_df), dtype=RAM_SAFE_SUM_DTYPE)
        for start in range(0, len(site_cols), RAM_SAFE_SITE_CHUNK):
            cols = [c for c in site_cols[start:start + RAM_SAFE_SITE_CHUNK] if c in final_df.columns]
            if not cols:
                continue
            block = final_df.loc[:, cols].to_numpy(dtype=RAM_SAFE_SUM_DTYPE, copy=True)
            totals += block.sum(axis=1, dtype=RAM_SAFE_SUM_DTYPE)
            del block
        final_df["Total_Abundance"] = totals
        del totals

    zero_mask = final_df["Total_Abundance"].to_numpy(copy=False) == 0
    n_zero = int(zero_mask.sum())

    if n_zero:
        print(f"   • Removing {n_zero} IDs with Total_Abundance == 0")
        final_df = _ram_safe_filter_wide_dataframe(
            source_df=final_df,
            keep_mask=~zero_mask,
            site_cols=site_cols,
            tag="zero_pruned"
        )
    else:
        print("   • No zeros to prune.")

    print(f"   • Remaining IDs: {len(final_df)}")
    gc.collect()
    return final_df

def load_datasets_with_stats_part2(report_lines):
    """
    RAM-safe Part 2 loader.

    Reads the filtered intermediate CSV files in small row chunks and stores
    abundance values in one preallocated uint32 matrix. This avoids the
    all_dfs -> unified -> concat chain of full DataFrames used by v1.7.12.
    """
    print("\n📁 PART 2 — LOADING FILES (RAM-SAFE)...")
    files = sorted(TARGET_DIR.glob(PATTERN))

    if not files:
        print(f"❌ No files found in: {TARGET_DIR} with pattern: {PATTERN}")
        return pd.DataFrame(), [], TARGET_DIR, PATTERN, None

    # Columns that are metadata, never sample-abundance columns.
    # NOTE: ``sequence`` is the canonical sequence field. Some input/intermediate
    # tables may also contain a legacy ``Sequence`` column; it must be ignored
    # rather than interpreted as a sample (it otherwise appears as a final
    # all-zero abundance column).
    fixed = {"OTU_XX", "Original_ID", "sintax_taxonomy", "sequence", "Sequence"}

    # First pass: discover all site columns and count rows cheaply.
    file_info = []
    site_union = set()
    total_rows = 0

    for csv_file in files:
        header = pd.read_csv(csv_file, nrows=0)
        current_sites = [c for c in header.columns if c not in fixed]
        site_union.update(current_sites)

        # Count rows using one lightweight column only.
        first_col = header.columns[0]
        nrows = 0
        for id_chunk in pd.read_csv(
            csv_file,
            usecols=[first_col],
            chunksize=max(RAM_SAFE_CHUNK_ROWS * 10, 50000),
            dtype=str
        ):
            nrows += len(id_chunk)

        subset = (
            csv_file.stem
            .replace("_only_bacteria", "")
            .replace("_all_prokaryotes", "")
            .replace("_only_metazoa", "")
            .replace("_fungi", "")
            .replace("_all_eukaryotes", "")
        )
        file_info.append((subset, csv_file, current_sites, nrows))
        total_rows += nrows

    site_cols_all = sorted(site_union)
    site_position = {c: i for i, c in enumerate(site_cols_all)}

    print(f"   • Preallocating abundance matrix: {total_rows:,} features x {len(site_cols_all):,} sites")
    abundance_matrix = np.zeros(
        (total_rows, len(site_cols_all)),
        dtype=RAM_SAFE_ABUNDANCE_DTYPE
    )

    otu_values = np.empty(total_rows, dtype=object)
    original_values = np.empty(total_rows, dtype=object)
    taxonomy_values = np.empty(total_rows, dtype=object)
    sequence_values = np.empty(total_rows, dtype=object)

    write_pos = 0

    for subset, csv_file, current_sites, expected_rows in file_info:
        print(f"   🧹 Loading {subset} in chunks...")
        subset_start = write_pos

        usecols = ["OTU_XX", "Original_ID", "sintax_taxonomy", "sequence"] + current_sites
        for chunk in pd.read_csv(
            csv_file,
            usecols=usecols,
            chunksize=RAM_SAFE_CHUNK_ROWS,
            low_memory=False
        ):
            n = len(chunk)
            sl = slice(write_pos, write_pos + n)

            otu_values[sl] = chunk["OTU_XX"].fillna("").astype(str).to_numpy(copy=False)
            original_values[sl] = chunk["Original_ID"].fillna("").astype(str).to_numpy(copy=False)
            taxonomy_clean = (
                chunk["sintax_taxonomy"]
                .fillna("")
                .astype(str)
                .map(clean_taxonomy_string)
            )
            taxonomy_values[sl] = taxonomy_clean.to_numpy(copy=False)
            sequence_values[sl] = chunk["sequence"].fillna("").astype(str).to_numpy(copy=False)

            if current_sites:
                numeric = (
                    chunk.loc[:, current_sites]
                    .apply(pd.to_numeric, errors="coerce")
                    .fillna(0)
                )
                arr = numeric.to_numpy(copy=False)
                if arr.size:
                    if np.nanmin(arr) < 0:
                        raise ValueError(f"Negative abundance detected in {csv_file.name}")
                    if np.nanmax(arr) > np.iinfo(RAM_SAFE_ABUNDANCE_DTYPE).max:
                        raise OverflowError(
                            f"Abundance exceeds {RAM_SAFE_ABUNDANCE_DTYPE} in {csv_file.name}. "
                            "Use np.uint64 for RAM_SAFE_ABUNDANCE_DTYPE."
                        )
                arr = arr.astype(RAM_SAFE_ABUNDANCE_DTYPE, copy=False)
                cols_idx = [site_position[c] for c in current_sites]
                abundance_matrix[sl, cols_idx] = arr
                del numeric, arr

            write_pos += n
            del chunk
            gc.collect()

        loaded_rows = write_pos - subset_start
        if loaded_rows != expected_rows:
            raise RuntimeError(
                f"Row-count mismatch while loading {csv_file.name}: "
                f"expected {expected_rows}, loaded {loaded_rows}."
            )
        report_lines.append(f"  {subset}: {loaded_rows} features, {len(current_sites)} sites")

    # Build a DataFrame directly over the compact uint32 matrix.
    combined = pd.DataFrame(
        abundance_matrix,
        columns=site_cols_all,
        copy=False
    )
    combined.insert(0, "sequence", sequence_values)
    combined.insert(0, "sintax_taxonomy", taxonomy_values)
    combined.insert(0, "Original_ID", original_values)
    combined.insert(0, "OTU_XX", otu_values)

    # Keep abundance_matrix alive: it is returned explicitly and reused by
    # collapse_otus_simple() to avoid rebuilding the full wide numeric matrix.
    # Only the temporary metadata arrays can be released here.
    del otu_values, original_values, taxonomy_values, sequence_values
    gc.collect()

    print("\n📈 PART 2 — COMBINED SUMMARY")
    print(f"   • Total features: {len(combined):,}")
    print(f"   • Total sites: {len(site_cols_all):,}")
    print(f"   • Abundance dtype: {RAM_SAFE_ABUNDANCE_DTYPE}")

    expected_shape = (len(combined), len(site_cols_all))
    if abundance_matrix.shape != expected_shape:
        raise RuntimeError(
            "Internal Part 2 abundance-matrix shape mismatch: "
            f"{abundance_matrix.shape} vs {expected_shape}."
        )

    return combined, site_cols_all, TARGET_DIR, PATTERN, abundance_matrix


def collapse_otus_simple(combined_df, site_cols, report_lines, abundance_matrix=None):
    """
    MetaDiv Builder v1.7.13 fast RAM-safe taxonomic collapse.

    Biological semantics are unchanged:
        representative = longest sequence
                         -> highest ORIGINAL total abundance
                         -> lexicographically smallest OTU_XX
                         -> lowest original row position (final deterministic tie)

    Performance changes in v1.7.13
    -----------------------------
    1. Sequence lengths, OTU labels and original row abundances are prepared once.
    2. Representative selection uses compact NumPy arrays rather than repeated
       ``DataFrame.at`` calls and repeated regular-expression cleaning.
    3. Group abundance aggregation is performed directly on the numeric NumPy
       abundance matrix returned by Part 2, avoiding one pandas wide-table
       selection per collapse group.
    4. Retained rows are compacted directly from that numeric matrix into the
       existing disk-backed memmap architecture.

    The taxonomic collapse keys, confidence thresholds, completeness rule and
    abundance-summing semantics are unchanged from v1.7.12/v4.1.
    """
    print(f"\n🔬 APPLYING COLLAPSE STRATEGY: {COLLAPSE_STRATEGY} (p≥{P_VALUE_THRESHOLD})")

    total = len(combined_df)
    valid_sites = [c for c in site_cols if c in combined_df.columns]

    # Production runs receive the numeric matrix directly from Part 2.
    # A compact fallback is kept for programmatic use of this function alone.
    if abundance_matrix is None:
        print("   ⚠️ Direct abundance matrix not supplied; creating a temporary numeric matrix.")
        abundance_matrix = combined_df.loc[:, valid_sites].to_numpy(
            dtype=RAM_SAFE_ABUNDANCE_DTYPE,
            copy=True
        )

    abundance_matrix = np.asarray(abundance_matrix)
    expected_shape = (total, len(valid_sites))
    if abundance_matrix.shape != expected_shape:
        raise ValueError(
            "Abundance matrix shape does not match the combined table: "
            f"{abundance_matrix.shape} vs {expected_shape}."
        )

    # ------------------------------------------------------------------
    # Taxonomic keys: preserve the established MetaDiv collapse semantics.
    # ------------------------------------------------------------------
    collapse_keys = np.empty(total, dtype=object)
    complete_flags = np.zeros(total, dtype=bool)

    taxonomy_series = combined_df["sintax_taxonomy"].fillna("").astype(str)
    for i, tx_string in enumerate(taxonomy_series):
        tax_data = extract_taxonomy_data(tx_string)
        collapse_keys[i] = build_collapse_key(
            tax_data, COLLAPSE_STRATEGY, P_VALUE_THRESHOLD
        )
        complete_flags[i] = has_complete_taxonomy(tx_string)

    taxonomy_values = taxonomy_series.to_numpy(copy=False)
    del taxonomy_series

    collapsable_mask = np.fromiter(
        (key is not None for key in collapse_keys),
        dtype=bool,
        count=total
    )
    collapsable_features = int(collapsable_mask.sum())
    non_collapsable = total - collapsable_features
    complete_in_collapsable = int(
        np.logical_and(collapsable_mask, complete_flags).sum()
    )

    group_index = defaultdict(list)
    for idx, key in enumerate(collapse_keys):
        if key is not None:
            group_index[key].append(idx)

    total_groups = len(group_index)
    multi_member_groups = sum(1 for idxs in group_index.values() if len(idxs) > 1)

    print("📊 FEATURE DISTRIBUTION:")
    print(
        f"   • Collapsable features: {collapsable_features:,} "
        f"({(collapsable_features/total*100 if total else 0):.1f}%)"
    )
    print(
        f"   • Complete among collapsable: {complete_in_collapsable:,} "
        f"({(complete_in_collapsable/collapsable_features*100 if collapsable_features else 0):.1f}%)"
    )
    print(
        f"   • Individual/non-collapsable: {non_collapsable:,} "
        f"({(non_collapsable/total*100 if total else 0):.1f}%)"
    )
    print(f"   • Taxonomic collapse groups: {total_groups:,}")
    print(f"   • Multi-feature groups requiring collapse: {multi_member_groups:,}")

    report_lines.append("\n  Collapse Statistics:")
    report_lines.append(f"    Features before collapse: {total:,}")
    report_lines.append(f"    Collapsable features: {collapsable_features:,}")
    report_lines.append(f"    Complete taxonomy in collapsable features: {complete_in_collapsable:,}")
    report_lines.append(f"    Taxonomic collapse groups: {total_groups:,}")
    report_lines.append(f"    Multi-feature groups requiring collapse: {multi_member_groups:,}")

    # ------------------------------------------------------------------
    # v1.7.13 FAST PATH: prepare representative-ranking data only once.
    # ------------------------------------------------------------------
    print("\n⚡ PRECOMPUTING REPRESENTATIVE-SELECTION ARRAYS...")

    sequence_lengths = (
        combined_df["sequence"]
        .fillna("")
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.len()
        .fillna(0)
        .to_numpy(dtype=np.int32, copy=True)
    )
    otu_values = (
        combined_df["OTU_XX"]
        .fillna("")
        .astype(str)
        .to_numpy(copy=True)
    )

    original_row_totals = np.zeros(total, dtype=RAM_SAFE_SUM_DTYPE)
    if valid_sites:
        for row_start in range(0, total, RAM_SAFE_CHUNK_ROWS):
            row_stop = min(row_start + RAM_SAFE_CHUNK_ROWS, total)
            original_row_totals[row_start:row_stop] = abundance_matrix[
                row_start:row_stop, :
            ].sum(axis=1, dtype=RAM_SAFE_SUM_DTYPE)

    print("   ✅ Sequence lengths and original abundance totals ready.")

    keep_mask = np.ones(total, dtype=bool)
    logs = []
    dtype_max = np.iinfo(RAM_SAFE_ABUNDANCE_DTYPE).max

    print(f"\n🔄 Processing {multi_member_groups:,} multi-feature collapse groups...")

    processed_multi = 0
    for key, idxs in group_index.items():
        if len(idxs) <= 1:
            continue

        processed_multi += 1
        idx_array = np.asarray(idxs, dtype=np.int64)
        complete_idx_array = idx_array[complete_flags[idx_array]]

        if complete_idx_array.size == 0:
            logs.append({
                "collapse_key": key,
                "kept_otu": None,
                "removed_count": 0,
                "removed_otus": [],
                "group_size": len(idxs),
                "representative_seq_len": None,
                "representative_tax": None,
                "note": f"skip_collapse_incomplete_tax_{COLLAPSE_STRATEGY}",
            })
            continue

        # Rule 1: longest sequence.
        candidate_lengths = sequence_lengths[complete_idx_array]
        max_length = int(candidate_lengths.max())
        tied = complete_idx_array[candidate_lengths == max_length]

        # Rule 2: if sequence length is tied, highest ORIGINAL abundance.
        if tied.size > 1:
            tied_abundance = original_row_totals[tied]
            max_abundance = tied_abundance.max()
            tied = tied[tied_abundance == max_abundance]

        # Rule 3: deterministic lexical OTU_XX, then original row position.
        if tied.size > 1:
            idx_keep = min(
                (int(idx) for idx in tied),
                key=lambda idx: (str(otu_values[idx]), idx)
            )
        else:
            idx_keep = int(tied[0])

        representative_seq_len = int(sequence_lengths[idx_keep])
        representative_abundance = int(original_row_totals[idx_keep])

        # Sum every member of the collapse group across every sample.
        # NumPy performs the wide numeric work directly; pandas is not called
        # inside the group loop.
        if valid_sites:
            summed = np.zeros(len(valid_sites), dtype=RAM_SAFE_SUM_DTYPE)

            for group_start in range(0, idx_array.size, RAM_SAFE_GROUP_ROWS):
                group_stop = min(group_start + RAM_SAFE_GROUP_ROWS, idx_array.size)
                group_rows = idx_array[group_start:group_stop]
                summed += abundance_matrix[group_rows, :].sum(
                    axis=0,
                    dtype=RAM_SAFE_SUM_DTYPE
                )

            if summed.size and summed.max() > dtype_max:
                raise OverflowError(
                    "Collapsed abundance exceeds uint32. "
                    "Set RAM_SAFE_ABUNDANCE_DTYPE = np.uint64."
                )

            abundance_matrix[idx_keep, :] = summed.astype(
                RAM_SAFE_ABUNDANCE_DTYPE,
                copy=False
            )
            del summed

        removed_idx_array = idx_array[idx_array != idx_keep]
        keep_mask[removed_idx_array] = False

        removed_otus = [
            str(otu_values[int(idx)])
            for idx in removed_idx_array[:10]
        ]

        logs.append({
            "collapse_key": key,
            "kept_otu": str(otu_values[idx_keep]),
            "removed_count": int(removed_idx_array.size),
            "removed_otus": removed_otus,
            "group_size": int(idx_array.size),
            "representative_seq_len": representative_seq_len,
            "representative_original_abundance": float(representative_abundance),
            "representative_selection_rule": "longest_sequence_then_highest_abundance_then_OTU_XX",
            "representative_tax": str(taxonomy_values[idx_keep]),
            "note": f"collapsed_with_complete_tax_{COLLAPSE_STRATEGY}",
        })

        if processed_multi % 500 == 0 or processed_multi == multi_member_groups:
            print(
                f"   • Multi-feature groups processed: "
                f"{processed_multi:,}/{multi_member_groups:,}"
            )

    del original_row_totals, sequence_lengths, otu_values, collapsable_mask
    del taxonomy_values

    # Compact directly from the numeric matrix. This avoids reading the wide
    # abundance block back through pandas after collapse.
    final_df = _ram_safe_filter_wide_dataframe_from_matrix(
        source_df=combined_df,
        abundance_matrix=abundance_matrix,
        keep_mask=keep_mask,
        site_cols=site_cols,
        tag="postcollapse"
    )

    del collapse_keys, complete_flags, keep_mask, group_index
    gc.collect()

    # Parse taxonomy only for retained features.
    n_final = len(final_df)
    tax_name_arrays = {
        rank: np.empty(n_final, dtype=object)
        for rank in TAX_RANKS
    }
    tax_p_arrays = {
        rank: np.zeros(n_final, dtype=np.float32)
        for rank in TAX_RANKS
    }

    for i, tx_string in enumerate(final_df["sintax_taxonomy"].fillna("").astype(str)):
        parsed = parse_tax_columns_from_sintax(tx_string)
        for rank in TAX_RANKS:
            tax_name_arrays[rank][i] = parsed.get(rank, "")
            tax_p_arrays[rank][i] = parsed.get(f"{rank}_pvalue", 0.0)

    for rank in TAX_RANKS:
        final_df[rank] = tax_name_arrays[rank]
        final_df[f"{rank}_pvalue"] = tax_p_arrays[rank]

    del tax_name_arrays, tax_p_arrays
    gc.collect()

    print("\n✅ COLLAPSE COMPLETED")
    print(f"   • Initial features: {total:,}")
    print(f"   • Final features:   {len(final_df):,}")
    print(f"   • Reduction: {total - len(final_df):,} features")

    final_df["sequence_lenght"] = (
        final_df["sequence"]
        .fillna("")
        .astype(str)
        .str.replace(r"\s+", "", regex=True)
        .str.len()
    )

    report_lines.append(f"    Features after collapse: {len(final_df):,}")
    report_lines.append(f"    Reduction: {total - len(final_df):,} features")
    return final_df, logs

def compute_total_abundance_and_sort(final_df, site_cols):
    """
    Compute Total_Abundance without physically sorting the full wide DataFrame.

    RAM-SAFE v1.7.13 processes row blocks across all sample columns. This avoids
    repeatedly selecting the complete 644k-row table in many tiny site-column
    blocks and substantially reduces pandas indexing overhead.
    """
    print("\n📊 COMPUTING TOTAL ABUNDANCE (RAM-SAFE v1.7.13; LOGICAL SORT ORDER)...")

    valid_sites = [c for c in site_cols if c in final_df.columns]
    totals = np.zeros(len(final_df), dtype=RAM_SAFE_SUM_DTYPE)

    if valid_sites:
        site_positions = final_df.columns.get_indexer(valid_sites)

        for row_start in range(0, len(final_df), RAM_SAFE_CHUNK_ROWS):
            row_stop = min(row_start + RAM_SAFE_CHUNK_ROWS, len(final_df))

            block = final_df.iloc[
                row_start:row_stop, site_positions
            ].to_numpy(
                dtype=RAM_SAFE_ABUNDANCE_DTYPE,
                copy=True
            )

            totals[row_start:row_stop] = block.sum(
                axis=1,
                dtype=RAM_SAFE_SUM_DTYPE
            )

            del block

            if row_start == 0 or row_stop == len(final_df) or row_stop % 100000 == 0:
                print(f"   • Total abundance: {row_stop:,}/{len(final_df):,} features")

    final_df["Total_Abundance"] = totals
    del totals

    print("   • Physical row sorting skipped to avoid a full wide-table copy.")
    print("   • Descending abundance order retained for SPPN assignment and FINAL_DB export.")

    create_binomial_species_inplace(final_df, MODE)
    gc.collect()
    return final_df


# ============================================================================
# EMBEDDED KRONA HTML RESOURCES
# ============================================================================
# Embedded from the provided Krona Excel template during development. Runtime
# requires no Excel, KronaTools, openpyxl, or additional package.

_KRONA_LEGACY_JS = '{//-----------------------------------------------------------------------------\n// \n// PURPOSE\n// \n// Krona is a flexible tool for exploring the relative proportions of\n// hierarchical data, such as metagenomic classifications, using a\n// radial, space-filling display. It is implemented using HTML5 and\n// JavaScript, allowing charts to be explored locally or served over the\n// Internet, requiring only a current version of any major web\n// browser. Krona charts can be created using an Excel template or from\n// common bioinformatic formats using the provided conversion scripts.\n// \n// \n// COPYRIGHT LICENSE\n// \n// Copyright (c) 2011, Battelle National Biodefense Institute (BNBI);\n// all rights reserved. Authored by: Brian Ondov, Nicholas Bergman, and\n// Adam Phillippy\n// \n// This Software was prepared for the Department of Homeland Security\n// (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as\n// part of contract HSHQDC-07-C-00020 to manage and operate the National\n// Biodefense Analysis and Countermeasures Center (NBACC), a Federally\n// Funded Research and Development Center.\n// \n// Redistribution and use in source and binary forms, with or without\n// modification, are permitted provided that the following conditions are\n// met:\n// \n// * Redistributions of source code must retain the above copyright\n//   notice, this list of conditions and the following disclaimer.\n// \n// * Redistributions in binary form must reproduce the above copyright\n//   notice, this list of conditions and the following disclaimer in the\n//   documentation and/or other materials provided with the distribution.\n// \n// * Neither the name of the Battelle National Biodefense Institute nor\n//   the names of its contributors may be used to endorse or promote\n//   products derived from this software without specific prior written\n//   permission.\n// \n// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n// \n// \n// TRADEMARK LICENSE\n// \n// KRONA(TM) is a trademark of the Department of Homeland Security, and use\n// of the trademark is subject to the following conditions:\n// \n// * Distribution of the unchanged, official code/software using the\n//   KRONA(TM) mark is hereby permitted by the Department of Homeland\n//   Security, provided that the software is distributed without charge\n//   and modification.\n// \n// * Distribution of altered source code/software using the KRONA(TM) mark\n//   is not permitted unless written permission has been granted by the\n//   Department of Homeland Security.\n// \n// \n// FOR MORE INFORMATION VISIT\n// \n// https://github.com/marbl/Krona/wiki/\n// \n//-----------------------------------------------------------------------------\n}\nvar canvas;\nvar context;\nvar svg; // for snapshot mode\nvar collapse = true;\nvar collapseCheckBox;\nvar collapseLast;\nvar compress;\nvar compressCheckBox;\nvar maxAbsoluteDepthText;\nvar maxAbsoluteDepthButtonDecrease;\nvar maxAbsoluteDepthButtonIncrease;\nvar fontSize = 11;\nvar fontSizeText;\nvar fontSizeButtonDecrease;\nvar fontSizeButtonIncrease;\nvar fontSizeLast;\nvar radiusButtonDecrease;\nvar radiusButtonIncrease;\nvar shorten;\nvar shortenCheckBox;\nvar maxAbsoluteDepth;\nvar backButton;\nvar upButton;\nvar forwardButton;\nvar snapshotButton;\nvar snapshotMode = false;\nvar details;\nvar detailsName;\nvar search;\nvar searchResults;\nvar nSearchResults;\nvar useHueCheckBox;\nvar useHueDiv;\nvar datasetDropDown;\nvar datasetButtonLast;\nvar datasetButtonPrev;\nvar datasetButtonNext;\nvar keyControl;\nvar showKeys = true;\nvar linkButton;\nvar linkText;\nvar frame;\n// Node references. Note that the meanings of \'selected\' and \'focused\' are\n// swapped in the docs.\n//\nvar head; // the root of the entire tree\nvar selectedNode = 0; // the root of the current view\nvar focusNode = 0; // a node chosen for more info (single-click)\nvar highlightedNode = 0; // mouse hover node\nvar highlightingHidden = false;\nvar nodes = new Array();\nvar currentNodeID = 0; // to iterate while loading\nvar nodeHistory = new Array();\nvar nodeHistoryPosition = 0;\nvar dataEnabled = false; // true when supplemental files are present\n// store non-Krona GET variables so they can be passed on to links\n//\nvar getVariables = new Array();\n// selectedNodeLast is separate from the history, since we need to check\n// properties of the last node viewed when browsing through the history\n//\nvar selectedNodeLast = 0;\nvar zoomOut = false;\n// temporary zoom-in while holding the mouse button on a wedge\n//\nvar quickLook = false; // true when in quick look state\nvar mouseDown = false;\nvar mouseDownTime; // to detect mouse button hold\nvar quickLookHoldLength = 200;\nvar imageWidth;\nvar imageHeight;\nvar centerX;\nvar centerY;\nvar gRadius;\nvar updateViewNeeded = false;\n// Determines the angle that the pie chart starts at.  90 degrees makes the\n// center label consistent with the children.\n//\nvar rotationOffset = Math.PI / 2;\nvar buffer;\nvar bufferFactor = .1;\n// The maps are the small pie charts showing the current slice being viewed.\n//\nvar mapBuffer = 10;\nvar mapRadius = 0;\nvar maxMapRadius = 25;\nvar mapWidth = 150;\nvar maxLabelOverhang = Math.PI * 4.18;\n// Keys are the labeled boxes for slices in the highest level that are too thin\n// to label.\n//\nvar maxKeySizeFactor = 2; // will be multiplied by font size\nvar keySize;\nvar keys;\nvar keyBuffer = 10;\nvar currentKey;\nvar keyMinTextLeft;\nvar keyMinAngle;\nvar minRingWidthFactor = 5; // will be multiplied by font size\nvar maxPossibleDepth; // the theoretical max that can be displayed\nvar maxDisplayDepth; // the actual depth that will be displayed\nvar headerHeight = 0;//document.getElementById(\'options\').clientHeight;\nvar historySpacingFactor = 1.6; // will be multiplied by font size\nvar historyAlphaDelta = .25;\n// appearance\n//\nvar lineOpacity = 0.3;\nvar saturation = 0.5;\nvar lightnessBase = 0.6;\nvar lightnessMax = .8;\nvar thinLineWidth = .3;\nvar highlightLineWidth = 1.5;\nvar labelBoxBuffer = 6;\nvar labelBoxRounding = 15;\nvar labelWidthFudge = 1.05; // The width of unshortened labels are set slightly\n\t\t\t\t\t\t\t// longer than the name width so the animation\n\t\t\t\t\t\t\t// finishes faster.\nvar fontNormal;\nvar fontBold;\nvar fontFamily = \'sans-serif\';\n//var fontFaceBold = \'bold Arial\';\nvar nodeRadius;\nvar angleFactor;\nvar tickLength;\nvar compressedRadii;\n// colors\n//\nvar highlightFill = \'rgba(255, 255, 255, .3)\';\nvar colorUnclassified = \'rgb(220,220,220)\';\n// label staggering\n//\nvar labelOffsets; // will store the current offset at each depth\n//\n// This will store pointers to the last node that had a label in each offset (or "track") of a\n// each depth.  These will be used to shorten neighboring labels that would overlap.\n// The [nLabelNodes] index will store the last node with a radial label.\n// labelFirstNodes is the same, but to check for going all the way around and\n// overlapping the first labels.\n//\nvar labelLastNodes;\nvar labelFirstNodes;\n//\nvar nLabelOffsets = 3; // the number of offsets to use\nvar mouseX = -1;\nvar mouseY = -1;\nvar mouseXRel = -1;\nvar mouseYRel = -1;\n// tweening\n//\nvar progress = 0; // for tweening; goes from 0 to 1.\nvar progressLast = 0;\nvar tweenFactor = 0; // progress converted by a curve for a smoother effect.\nvar tweenLength = 850; // in ms\nvar tweenCurvature = 13;\n//\n// tweenMax is used to scale the sigmoid function so its range is [0,1] for the\n// domain [0,1]\n//\nvar tweenMax = 1 / (1 + Math.exp(-tweenCurvature / 2));\n//\nvar tweenStartTime;\n// for framerate debug\n//\nvar tweenFrames = 0;\nvar fpsDisplay = document.getElementById(\'frameRate\');\n// Arrays to translate xml attribute names into displayable attribute names\n//\nvar attributes = new Array();\n//\nvar magnitudeIndex; // the index of attribute arrays used for magnitude\nvar membersAssignedIndex;\nvar membersSummaryIndex;\n// For defining gradients\n//\nvar hueDisplayName;\nvar hueStopPositions;\nvar hueStopHues;\nvar hueStopText;\n// multiple datasets\n//\nvar currentDataset = 0;\nvar lastDataset = 0;\nvar datasets = 1;\nvar datasetNames;\nvar datasetSelectSize = 30;\nvar datasetAlpha = new Tween(0, 0);\nvar datasetWidths = new Array();\nvar datasetChanged;\nvar datasetSelectWidth = 50;\nwindow.onload = load;\nvar image;\nvar hiddenPattern;\nvar loadingImage;\nvar logoImage;\nfunction backingScale()\n{\n\tif (\'devicePixelRatio\' in window)\n\t{\n\t\tif (window.devicePixelRatio > 1)\n\t\t{\n\t\t\treturn window.devicePixelRatio;\n\t\t}\n\t}\n\t\n\treturn 1;\n}\nfunction resize()\n{\n\timageWidth = window.innerWidth;\n\timageHeight = window.innerHeight;\n\t\n\tif ( ! snapshotMode )\n\t{\n\t\tcontext.canvas.width = imageWidth * backingScale();\n\t\tcontext.canvas.height = imageHeight * backingScale();\n\t\tcontext.canvas.style.width = imageWidth + "px"\n\t\tcontext.canvas.style.height = imageHeight + "px"\n\t\tcontext.scale(backingScale(), backingScale());\n\t}\n\t\n\tif ( datasetDropDown )\n\t{\n\t\tvar ratio = \n\t\t\t(datasetDropDown.offsetTop + datasetDropDown.clientHeight) * 2 /\n\t\t\timageHeight;\n\t\t\n\t\tif ( ratio > 1 )\n\t\t{\n\t\t\tratio = 1;\n\t\t}\n\t\t\n\t\tratio = Math.sqrt(ratio);\n\t\t\n\t\tdatasetSelectWidth = \n\t\t\t(datasetDropDown.offsetLeft + datasetDropDown.clientWidth) * ratio;\n\t}\n\tvar leftMargin = datasets > 1 ? datasetSelectWidth + 30 : 0;\n\tvar minDimension = imageWidth - mapWidth - leftMargin > imageHeight ?\n\t\timageHeight :\n\t\timageWidth - mapWidth - leftMargin;\n\t\n\tmaxMapRadius = minDimension * .03;\n\tbuffer = minDimension * bufferFactor;\n\tmargin = minDimension * .015;\n\tcenterX = (imageWidth - mapWidth - leftMargin) / 2 + leftMargin;\n\tcenterY = imageHeight / 2;\n\tgRadius = minDimension / 2 - buffer;\n\t//context.font = \'11px sans-serif\';\n}\nfunction handleResize()\n{\n\tupdateViewNeeded = true;\n}\nfunction Attribute()\n{\n}\nfunction Tween(start, end)\n{\n\tthis.start = start;\n\tthis.end = end;\n\tthis.current = this.start;\n\t\n\tthis.current = function()\n\t{\n\t\tif ( progress == 1 || this.start == this.end )\n\t\t{\n\t\t\treturn this.end;\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this.start + tweenFactor * (this.end - this.start);\n\t\t}\n\t};\n\t\n\tthis.setTarget = function(target)\n\t{\n\t\tthis.start = this.current();\n\t\tthis.end = target;\n\t}\n}\nfunction Node()\n{\n\tthis.id = currentNodeID;\n\tcurrentNodeID++;\n\tnodes[this.id] = this;\n\t\n\tthis.angleStart = new Tween(Math.PI, 0);\n\tthis.angleEnd = new Tween(Math.PI, 0);\n\tthis.radiusInner = new Tween(1, 1);\n\tthis.labelRadius = new Tween(1, 1);\n\tthis.labelWidth = new Tween(0, 0);\n\tthis.scale = new Tween(1, 1); // TEMP\n\tthis.radiusOuter = new Tween(1, 1);\n\t\n\tthis.r = new Tween(255, 255);\n\tthis.g = new Tween(255, 255);\n\tthis.b = new Tween(255, 255);\n\t\n\tthis.alphaLabel = new Tween(0, 1);\n\tthis.alphaLine = new Tween(0, 1);\n\tthis.alphaArc = new Tween(0, 0);\n\tthis.alphaWedge = new Tween(0, 1);\n\tthis.alphaOther = new Tween(0, 1);\n\tthis.alphaPattern = new Tween(0, 0);\n\tthis.children = Array();\n\tthis.parent = 0;\n\t\n\tthis.attributes = new Array(attributes.length);\n\t\n\tthis.addChild = function(child)\n\t{\n\t\tthis.children.push(child);\n\t};\n\t\n\tthis.addLabelNode = function(depth, labelOffset)\n\t{\n\t\tif ( labelHeadNodes[depth][labelOffset] == 0 )\n\t\t{\n\t\t\t// this will become the head node for this list\n\t\t\t\n\t\t\tlabelHeadNodes[depth][labelOffset] = this;\n\t\t\tthis.labelPrev = this;\n\t\t}\n\t\t\n\t\tvar head = labelHeadNodes[depth][labelOffset];\n\t\t\n\t\tthis.labelNext = head;\n\t\tthis.labelPrev = head.labelPrev;\n\t\thead.labelPrev.labelNext = this;\n\t\thead.labelPrev = this;\n\t}\n\t\n\tthis.canDisplayDepth = function()\n\t{\n\t\t// whether this node is at a depth that can be displayed, according\n\t\t// to the max absolute depth\n\t\t\n\t\treturn this.depth <= maxAbsoluteDepth;\n\t}\n\t\n\tthis.canDisplayHistory = function()\n\t{\n\t\tvar radiusInner;\n\t\t\n\t\tif ( compress )\n\t\t{\n\t\t\tradiusInner = compressedRadii[0];\n\t\t}\n\t\telse\n\t\t{\n\t\t\tradiusInner = nodeRadius;\n\t\t}\n\t\t\n\t\treturn (\n\t\t\t-this.labelRadius.end * gRadius +\n\t\t\thistorySpacingFactor * fontSize / 2 <\n\t\t\tradiusInner * gRadius\n\t\t\t);\n\t}\n\t\n\tthis.canDisplayLabelCurrent = function()\n\t{\n\t\treturn (\n\t\t\t(this.angleEnd.current() - this.angleStart.current()) *\n\t\t\t(this.radiusInner.current() * gRadius + gRadius) >=\n\t\t\tminWidth());\n\t}\n\t\n\tthis.checkHighlight = function()\n\t{\n\t\tif ( this.children.length == 0 && this == focusNode )\n\t\t{\n\t\t\t//return false;\n\t\t}\n\t\t\n\t\tif ( this.hide )\n\t\t{\n\t\t\treturn false;\n\t\t}\n\t\t\n\t\tif ( this.radiusInner.end == 1 )\n\t\t{\n\t\t\t// compressed to the outside; don\'t check\n\t\t\t\n\t\t\treturn false;\n\t\t}\n\t\t\n\t\tvar highlighted = false;\n\t\t\n\t\tvar angleStartCurrent = this.angleStart.current() + rotationOffset;\n\t\tvar angleEndCurrent = this.angleEnd.current() + rotationOffset;\n\t\tvar radiusInner = this.radiusInner.current() * gRadius;\n\t\t\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\thighlighted = this.children[i].checkHighlight();\n\t\t\t\n\t\t\tif ( highlighted )\n\t\t\t{\n\t\t\t\treturn true;\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( this.radial )\n\t\t{\n\t\t\tvar angleText = (angleStartCurrent + angleEndCurrent) / 2;\n\t\t\tvar radiusText = (gRadius + radiusInner) / 2;\n\t\t\t\n\t\t\tcontext.rotate(angleText);\n\t\t\tcontext.beginPath();\n\t\t\tcontext.moveTo(radiusText, -fontSize);\n\t\t\tcontext.lineTo(radiusText, fontSize);\n\t\t\tcontext.lineTo(radiusText + centerX, fontSize);\n\t\t\tcontext.lineTo(radiusText + centerX, -fontSize);\n\t\t\tcontext.closePath();\n\t\t\tcontext.rotate(-angleText);\n\t\t\t\n\t\t\tif ( context.isPointInPath(mouseXRel, mouseYRel) )\n\t\t\t{\n\t\t\t\tvar label = String(this.getPercentage()) + \'%\' + \'   \' + this.name;\n\t\t\t\t\n\t\t\t\tif ( this.searchResultChildren() )\n\t\t\t    {\n\t\t\t\t\tlabel += searchResultString(this.searchResultChildren());\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\tMath.sqrt((mouseXRel) * (mouseXRel) + (mouseYRel) * (mouseYRel)) / backingScale() <\n\t\t\t\t\tradiusText + measureText(label)\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\thighlighted = true;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t    for ( var i = 0; i < this.hiddenLabels.length; i++ )\n\t\t    {\n\t\t        var hiddenLabel = this.hiddenLabels[i];\n\t\t        \n\t\t\t\tcontext.rotate(hiddenLabel.angle);\n\t\t\t\tcontext.beginPath();\n\t\t\t\tcontext.moveTo(gRadius, -fontSize);\n\t\t\t\tcontext.lineTo(gRadius, fontSize);\n\t\t\t\tcontext.lineTo(gRadius + centerX, fontSize);\n\t\t\t\tcontext.lineTo(gRadius + centerX, -fontSize);\n\t\t\t\tcontext.closePath();\n\t\t\t\tcontext.rotate(-hiddenLabel.angle);\n\t\t\t\t\n\t\t\t\tif ( context.isPointInPath(mouseXRel, mouseYRel) )\n\t\t\t\t{\n\t\t\t\t\tvar label = String(hiddenLabel.value) + \' more\';\n\t\t\t\t\t\n\t\t\t\t\tif ( hiddenLabel.search )\n\t\t\t\t    {\n\t\t\t\t\t\tlabel += searchResultString(hiddenLabel.search);\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tif\n\t\t\t\t\t(\n\t\t\t\t\t\tMath.sqrt((mouseXRel) * (mouseXRel) + (mouseYRel) * (mouseYRel)) / backingScale() <\n\t\t\t\t\t\tgRadius + fontSize + measureText(label)\n\t\t\t\t\t)\n\t\t\t\t\t{\n\t\t\t\t\t\thighlighted = true;\n\t\t\t\t\t\tbreak;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( ! highlighted && this != selectedNode && ! this.getCollapse() )\n\t\t{\n\t\t\tcontext.beginPath();\n\t\t\tcontext.arc(0, 0, radiusInner, angleStartCurrent, angleEndCurrent, false);\n\t\t\tcontext.arc(0, 0, gRadius, angleEndCurrent, angleStartCurrent, true);\n\t\t\tcontext.closePath();\n\t\t\t\n\t\t\tif ( context.isPointInPath(mouseXRel, mouseYRel) )\n\t\t\t{\n\t\t\t\thighlighted = true;\n\t\t\t}\n\t\t\t\n\t\t\tif\n\t\t\t(\n\t\t\t\t! highlighted &&\n\t\t\t\t(angleEndCurrent - angleStartCurrent) *\n\t\t\t\t(radiusInner + gRadius) <\n\t\t\t\tminWidth() &&\n\t\t\t\tthis.getDepth() == selectedNode.getDepth() + 1\n\t\t\t)\n\t\t\t{\n\t\t\t\tif ( showKeys && this.checkHighlightKey() )\n\t\t\t\t{\n\t\t\t\t\thighlighted = true;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( highlighted )\n\t\t{\n\t\t\tif ( this != highlightedNode )\n\t\t\t{\n\t\t\t//\tdocument.body.style.cursor=\'pointer\';\n\t\t\t}\n\t\t\t\n\t\t\thighlightedNode = this;\n\t\t}\n\t\t\n\t\treturn highlighted;\n\t}\n\t\n\tthis.checkHighlightCenter = function()\n\t{\n\t\tif ( ! this.canDisplayHistory() )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar cx = centerX;\n\t\tvar cy = centerY - this.labelRadius.end * gRadius;\n\t\t//var dim = context.measureText(this.name);\n\t\t\n\t\tvar width = this.nameWidth;\n\t\t\n\t\tif ( this.searchResultChildren() )\n\t\t{\n\t\t\tvar results = searchResultString(this.searchResultChildren());\n\t\t\tvar dim = context.measureText(results);\n\t\t\twidth += dim.width;\n\t\t}\n\t\t\n\t\tif\n\t\t(\n\t\t\tmouseX > cx - width / 2 &&\n\t\t\tmouseX < cx + width / 2 &&\n\t\t\tmouseY > cy - historySpacingFactor * fontSize / 2 &&\n\t\t\tmouseY < cy + historySpacingFactor * fontSize / 2\n\t\t)\n\t\t{\n\t\t\thighlightedNode = this;\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tif ( this.getParent() )\n\t\t{\n\t\t\tthis.getParent().checkHighlightCenter();\n\t\t}\n\t}\n\t\n\tthis.checkHighlightKey = function()\n\t{\n\t\tvar offset = keyOffset();\n\t\t\n\t\tvar xMin = imageWidth - keySize - margin - this.keyNameWidth - keyBuffer;\n\t\tvar xMax = imageWidth - margin;\n\t\tvar yMin = offset;\n\t\tvar yMax = offset + keySize;\n\t\t\n\t\tcurrentKey++;\n\t\t\n\t\treturn (\n\t\t\tmouseX > xMin &&\n\t\t\tmouseX < xMax &&\n\t\t\tmouseY > yMin &&\n\t\t\tmouseY < yMax);\n\t}\n\t\n\tthis.checkHighlightMap = function()\n\t{\n\t\tif ( this.parent )\n\t\t{\n\t\t\tthis.parent.checkHighlightMap();\n\t\t}\n\t\t\n\t\tif ( this.getCollapse() || this == focusNode )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar box = this.getMapPosition();\n\t\t\n\t\tif\n\t\t(\n\t\t\tmouseX > box.x - mapRadius &&\n\t\t\tmouseX < box.x + mapRadius &&\n\t\t\tmouseY > box.y - mapRadius &&\n\t\t\tmouseY < box.y + mapRadius\n\t\t)\n\t\t{\n\t\t\thighlightedNode = this;\n\t\t}\n\t}\n\t\n/*\tthis.collapse = function()\n\t{\n\t\tfor (var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.children[i] = this.children[i].collapse();\n\t\t}\n\t\t\n\t\tif\n\t\t(\n\t\t\tthis.children.length == 1 &&\n\t\t\tthis.children[0].magnitude == this.magnitude\n\t\t)\n\t\t{\n\t\t\tthis.children[0].parent = this.parent;\n\t\t\tthis.children[0].getDepth() = this.parent.getDepth() + 1;\n\t\t\treturn this.children[0];\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this;\n\t\t}\n\t}\n*/\t\n\tthis.draw = function(labelMode, selected, searchHighlighted)\n\t{\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n//\t\tvar hidden = false;\n\t\t\n\t\tif ( selectedNode == this )\n\t\t{\n\t\t\tselected = true;\n\t\t}\n\t\t\n\t\tvar angleStartCurrent = this.angleStart.current() + rotationOffset;\n\t\tvar angleEndCurrent = this.angleEnd.current() + rotationOffset;\n\t\tvar radiusInner = this.radiusInner.current() * gRadius;\n\t\tvar canDisplayLabelCurrent = this.canDisplayLabelCurrent();\n\t\tvar hiddenSearchResults = false;\n\t\t\n/*\t\tif ( ! this.hide )\n\t\t{\n\t\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t\t{\n\t\t\t\tif ( this.children[i].hide && this.children[i].searchResults )\n\t\t\t\t{\n\t\t\t\t\thiddenSearchResults = true;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n*/\t\t\n\t\tvar drawChildren =\n\t\t\t( ! this.hide || ! this.hidePrev && progress < 1 ) &&\n\t\t\t( ! this.hideAlone || ! this.hideAlonePrev && progress < 1 );\n\t\t\n//\t\tif ( this.alphaWedge.current() > 0 || this.alphaLabel.current() > 0 )\n\t\t{\n\t\t\tvar lastChildAngleEnd = angleStartCurrent;\n\t\t\t\n\t\t\tif ( this.hasChildren() )//canDisplayChildren )\n\t\t\t{\n\t\t\t\tlastChildAngleEnd =\n\t\t\t\t\tthis.children[this.children.length - 1].angleEnd.current()\n\t\t\t\t\t+ rotationOffset;\n\t\t\t}\n\t\t\t\n\t\t\tif ( labelMode )\n\t\t\t{\n\t\t\t\tvar drawRadial =\n\t\t\t\t!(\n\t\t\t\t\tthis.parent &&\n\t\t\t\t\tthis.parent != selectedNode &&\n\t\t\t\t\tangleEndCurrent == this.parent.angleEnd.current() + rotationOffset\n\t\t\t\t);\n\t\t\t\t\n\t\t\t\t//if ( angleStartCurrent != angleEndCurrent )\n\t\t\t\t{\n\t\t\t\t\tthis.drawLines(angleStartCurrent, angleEndCurrent, radiusInner, drawRadial, selected);\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tvar alphaOtherCurrent = this.alphaOther.current();\n\t\t\t\tvar childRadiusInner;\n\t\t\t\t\n\t\t\t\tif ( this == selectedNode || alphaOtherCurrent )\n\t\t\t\t{\n\t\t\t\t\tchildRadiusInner =\n\t\t\t\t\t\tthis.children.length ?\n\t\t\t\t\t\t\tthis.children[this.children.length - 1].radiusInner.current() * gRadius\n\t\t\t\t\t\t: radiusInner\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( this == selectedNode )\n\t\t\t\t{\n\t\t\t\t\tthis.drawReferenceRings(childRadiusInner);\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\tselected &&\n\t\t\t\t\t! searchHighlighted &&\n\t\t\t\t\tthis != selectedNode &&\n\t\t\t\t\t(\n\t\t\t\t\t\tthis.isSearchResult ||\n\t\t\t\t\t\tthis.hideAlone && this.searchResultChildren() ||\n\t\t\t\t\t\tfalse\n//\t\t\t\t\t\tthis.hide &&\n//\t\t\t\t\t\tthis.containsSearchResult\n\t\t\t\t\t)\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tcontext.globalAlpha = this.alphaWedge.current();\n\t\t\t\t\t\n\t\t\t\t\tdrawWedge\n\t\t\t\t\t(\n\t\t\t\t\t\tangleStartCurrent,\n\t\t\t\t\t\tangleEndCurrent,\n\t\t\t\t\t\tradiusInner,\n\t\t\t\t\t\tgRadius,\n\t\t\t\t\t\thighlightFill,\n\t\t\t\t\t\t0,\n\t\t\t\t\t\ttrue\n\t\t\t\t\t);\n\t\t\t\t\t\n\t\t\t\t\tif\n\t\t\t\t\t(\n\t\t\t\t\t\tthis.keyed &&\n\t\t\t\t\t\t! showKeys &&\n\t\t\t\t\t\tthis.searchResults &&\n\t\t\t\t\t\t! searchHighlighted &&\n\t\t\t\t\t\tthis != highlightedNode &&\n\t\t\t\t\t\tthis != focusNode\n\t\t\t\t\t)\n\t\t\t\t\t{\n\t\t\t\t\t\tvar angle = (angleEndCurrent + angleStartCurrent) / 2;\n\t\t\t\t\t\tthis.drawLabel(angle, true, false, true, true);\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\t//this.drawHighlight(false);\n\t\t\t\t\tsearchHighlighted = true;\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\tthis == selectedNode ||\n//\t\t\t\t\ttrue\n\t\t\t\t\t//(canDisplayLabelCurrent) &&\n\t\t\t\t\tthis != highlightedNode &&\n\t\t\t\t\tthis != focusNode\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tif ( this.radial != this.radialPrev && this.alphaLabel.end == 1 )\n\t\t\t\t\t{\n\t\t\t\t\t\tcontext.globalAlpha = tweenFactor;\n\t\t\t\t\t}\n\t\t\t\t\telse\n\t\t\t\t\t{\n\t\t\t\t\t\tcontext.globalAlpha = this.alphaLabel.current();\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tthis.drawLabel\n\t\t\t\t\t(\n\t\t\t\t\t\t(angleStartCurrent + angleEndCurrent) / 2,\n\t\t\t\t\t\tthis.hideAlone && this.searchResultChildren() ||\n\t\t\t\t\t\t(this.isSearchResult || hiddenSearchResults) && selected,\n\t\t\t\t\t\tthis == selectedNode && ! this.radial,\n\t\t\t\t\t\tselected,\n\t\t\t\t\t\tthis.radial\n\t\t\t\t\t);\n\t\t\t\t\t\n\t\t\t\t\tif ( this.radial != this.radialPrev && this.alphaLabel.start == 1 && progress < 1 )\n\t\t\t\t\t{\n\t\t\t\t\t\tcontext.globalAlpha = 1 - tweenFactor;\n\t\t\t\t\t\t\n\t\t\t\t\t\tthis.drawLabel\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\t(angleStartCurrent + angleEndCurrent) / 2,\n\t\t\t\t\t\t\t(this.isSearchResult || hiddenSearchResults) && selected,\n\t\t\t\t\t\t\tthis == selectedNodeLast && ! this.radialPrev,\n\t\t\t\t\t\t\tselected,\n\t\t\t\t\t\t\tthis.radialPrev\n\t\t\t\t\t\t);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\talphaOtherCurrent &&\n\t\t\t\t\tlastChildAngleEnd != null\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tif\n\t\t\t\t\t(\n\t\t\t\t\t\t(angleEndCurrent - lastChildAngleEnd) *\n\t\t\t\t\t\t(childRadiusInner + gRadius) >=\n\t\t\t\t\t\tminWidth()\n\t\t\t\t\t)\n\t\t\t\t\t{\n\t\t\t\t\t\t//context.font = fontNormal;\n\t\t\t\t\t\tcontext.globalAlpha = this.alphaOther.current();\n\t\t\t\t\t\t\n\t\t\t\t\t\tdrawTextPolar\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\tthis.getUnclassifiedText(),\n\t\t\t\t\t\t\tthis.getUnclassifiedPercentage(),\n\t\t\t\t\t\t\t(lastChildAngleEnd + angleEndCurrent) / 2,\n\t\t\t\t\t\t\t(childRadiusInner + gRadius) / 2,\n\t\t\t\t\t\t\ttrue,\n\t\t\t\t\t\t\tfalse,\n\t\t\t\t\t\t\tfalse,\n\t\t\t\t\t\t\t0,\n\t\t\t\t\t\t\t0\n\t\t\t\t\t\t);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( this == selectedNode && this.keyUnclassified && showKeys )\n\t\t\t\t{\n\t\t\t\t\tthis.drawKey\n\t\t\t\t\t(\n\t\t\t\t\t\t(lastChildAngleEnd + angleEndCurrent) / 2,\n\t\t\t\t\t\tfalse,\n\t\t\t\t\t\tfalse\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tvar alphaWedgeCurrent = this.alphaWedge.current();\n\t\t\t\t\n\t\t\t\tif ( alphaWedgeCurrent || this.alphaOther.current() )\n\t\t\t\t{\n\t\t\t\t\tvar currentR = this.r.current();\n\t\t\t\t\tvar currentG = this.g.current();\n\t\t\t\t\tvar currentB = this.b.current();\n\t\t\t\t\t\t\n\t\t\t\t\tvar fill = rgbText(currentR, currentG, currentB);\n\t\t\t\t\t\n\t\t\t\t\tvar radiusOuter;\n\t\t\t\t\tvar lastChildAngle;\n\t\t\t\t\tvar truncateWedge =\n\t\t\t\t\t(\n\t\t\t\t\t\t(this.hasChildren() || this == selectedNode ) &&\n\t\t\t\t\t\t! this.keyed &&\n\t\t\t\t\t\t(compress || depth < maxDisplayDepth) &&\n\t\t\t\t\t\tdrawChildren\n\t\t\t\t\t);\n\t\t\t\t\t\n\t\t\t\t\tif ( truncateWedge )\n\t\t\t\t\t{\n\t\t\t\t\t\tradiusOuter = this.children.length ? this.children[0].radiusInner.current() * gRadius : radiusInner;\n\t\t\t\t\t}\n\t\t\t\t\telse\n\t\t\t\t\t{\n\t\t\t\t\t\tradiusOuter = gRadius;\n\t\t\t\t\t}\n\t\t\t\t\t/*\n\t\t\t\t\tif ( this.hasChildren() )\n\t\t\t\t\t{\n\t\t\t\t\t\tradiusOuter = this.children[0].getUncollapsed().radiusInner.current() * gRadius + 1;\n\t\t\t\t\t}\n\t\t\t\t\telse\n\t\t\t\t\t{ // TEMP\n\t\t\t\t\t\tradiusOuter = radiusInner + nodeRadius * gRadius;\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( radiusOuter > gRadius )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tradiusOuter = gRadius;\n\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t\t*/\n\t\t\t\t\tcontext.globalAlpha = alphaWedgeCurrent;\n\t\t\t\t\t\n\t\t\t\t\tif ( radiusInner != radiusOuter || truncateWedge )\n\t\t\t\t\t{\n\t\t\t\t\t\tdrawWedge\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\tangleStartCurrent,\n\t\t\t\t\t\t\tangleEndCurrent,\n\t\t\t\t\t\t\tradiusInner,\n\t\t\t\t\t\t\tradiusOuter,//this.radiusOuter.current() * gRadius,\n\t\t\t\t\t\t\t//\'rgba(0, 200, 0, .1)\',\n\t\t\t\t\t\t\tfill,\n\t\t\t\t\t\t\tthis.alphaPattern.current()\n\t\t\t\t\t\t);\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( truncateWedge )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\t// fill in the extra space if the sum of our childrens\'\n\t\t\t\t\t\t\t// magnitudes is less than ours\n\t\t\t\t\t\t\t\n\t\t\t\t\t\t\tif ( lastChildAngleEnd < angleEndCurrent )//&& false) // TEMP\n\t\t\t\t\t\t\t{\n\t\t\t\t\t\t\t\tif ( radiusOuter > 1 )\n\t\t\t\t\t\t\t\t{\n\t\t\t\t\t\t\t\t\t// overlap slightly to hide the seam\n\t\t\t\t\t\t\t\t\t\n\t//\t\t\t\t\t\t\t\tradiusOuter -= 1;\n\t\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\t\t\n\t\t\t\t\t\t\t\tif ( alphaWedgeCurrent < 1 )\n\t\t\t\t\t\t\t\t{\n\t\t\t\t\t\t\t\t\tcontext.globalAlpha = this.alphaOther.current();\n\t\t\t\t\t\t\t\t\tdrawWedge\n\t\t\t\t\t\t\t\t\t(\n\t\t\t\t\t\t\t\t\t\tlastChildAngleEnd,\n\t\t\t\t\t\t\t\t\t\tangleEndCurrent,\n\t\t\t\t\t\t\t\t\t\tradiusOuter,\n\t\t\t\t\t\t\t\t\t\tgRadius,\n\t\t\t\t\t\t\t\t\t\tcolorUnclassified,\n\t\t\t\t\t\t\t\t\t\t0\n\t\t\t\t\t\t\t\t\t);\n\t\t\t\t\t\t\t\t\tcontext.globalAlpha = alphaWedgeCurrent;\n\t\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\t\t\n\t\t\t\t\t\t\t\tdrawWedge\n\t\t\t\t\t\t\t\t(\n\t\t\t\t\t\t\t\t\tlastChildAngleEnd,\n\t\t\t\t\t\t\t\t\tangleEndCurrent,\n\t\t\t\t\t\t\t\t\tradiusOuter,\n\t\t\t\t\t\t\t\t\tgRadius,//this.radiusOuter.current() * gRadius,\n\t\t\t\t\t\t\t\t\t//\'rgba(200, 0, 0, .1)\',\n\t\t\t\t\t\t\t\t\tfill,\n\t\t\t\t\t\t\t\t\tthis.alphaPattern.current()\n\t\t\t\t\t\t\t\t);\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( radiusOuter < gRadius )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\t// patch up the seam\n\t\t\t\t\t\t\t//\n\t\t\t\t\t\t\tcontext.beginPath();\n\t\t\t\t\t\t\tcontext.arc(0, 0, radiusOuter, angleStartCurrent/*lastChildAngleEnd*/, angleEndCurrent, false);\n\t\t\t\t\t\t\tcontext.strokeStyle = fill;\n\t\t\t\t\t\t\tcontext.lineWidth = 1;\n\t\t\t\t\t\t\tcontext.stroke();\n\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tif ( this.keyed && selected && showKeys )//&& progress == 1 )\n\t\t\t\t\t{\n\t\t\t\t\t\tthis.drawKey\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\t(angleStartCurrent + angleEndCurrent) / 2,\n\t\t\t\t\t\t\t(\n\t\t\t\t\t\t\t\tthis == highlightedNode ||\n\t\t\t\t\t\t\t\tthis == focusNode ||\n\t\t\t\t\t\t\t\tthis.searchResults\n\t\t\t\t\t\t\t),\n\t\t\t\t\t\t\tthis == highlightedNode || this == focusNode\n\t\t\t\t\t\t);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n\t\tthis.hiddenLabels = Array();\n\t\t\n\t\tif ( drawChildren )\n\t\t{\n\t\t\t// draw children\n\t\t\t//\n\t\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t\t{\n\t\t\t\tif ( this.drawHiddenChildren(i, selected, labelMode, searchHighlighted) )\n\t\t\t\t{\n\t\t\t\t\ti = this.children[i].hiddenEnd;\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.children[i].draw(labelMode, selected, searchHighlighted);\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t};\n\t\n\tthis.drawHiddenChildren = function\n\t(\n\t\tfirstHiddenChild,\n\t\tselected,\n\t\tlabelMode,\n\t\tsearchHighlighted\n\t)\n\t{\n\t\tvar firstChild = this.children[firstHiddenChild];\n\t\t\n\t\tif ( firstChild.hiddenEnd == null || firstChild.radiusInner.current() == 1 )\n\t\t{\n\t\t\treturn false;\n\t\t}\n\t\t\n\t\tfor ( var i = firstHiddenChild; i < firstChild.hiddenEnd; i++ )\n\t\t{\n\t\t\tif ( ! this.children[i].hide || ! this.children[i].hidePrev && progress < 1 )\n\t\t\t{\n\t\t\t\treturn false;\n\t\t\t}\n\t\t}\n\t\t\n\t\tvar angleStart = firstChild.angleStart.current() + rotationOffset;\n\t\tvar lastChild = this.children[firstChild.hiddenEnd];\n\t\tvar angleEnd = lastChild.angleEnd.current() + rotationOffset;\n\t\tvar radiusInner = gRadius * firstChild.radiusInner.current();\n\t\tvar hiddenChildren = firstChild.hiddenEnd - firstHiddenChild + 1;\n\t\t\n\t\tif ( labelMode )\n\t\t{\n\t\t\tvar hiddenSearchResults = 0;\n\t\t\t\n\t\t\tfor ( var i = firstHiddenChild; i <= firstChild.hiddenEnd; i++ )\n\t\t\t{\n\t\t\t\thiddenSearchResults += this.children[i].searchResults;\n\t\t\t\t\n\t\t\t\tif ( this.children[i].magnitude == 0 )\n\t\t\t\t{\n\t\t\t\t\thiddenChildren--;\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\tif\n\t\t\t(\n\t\t\t\tselected &&\n\t\t\t\t(angleEnd - angleStart) * \n\t\t\t\t(gRadius + gRadius) >=\n\t\t\t\tminWidth() ||\n\t\t\t\tthis == highlightedNode &&\n\t\t\t\thiddenChildren ||\n\t\t\t\thiddenSearchResults\n\t\t\t)\n\t\t\t{\n\t\t\t\tcontext.globalAlpha = this.alphaWedge.current();\n\t\t\t\t\n\t\t\t\tthis.drawHiddenLabel\n\t\t\t\t(\n\t\t\t\t\tangleStart,\n\t\t\t\t\tangleEnd,\n\t\t\t\t\thiddenChildren,\n\t\t\t\t\thiddenSearchResults\n\t\t\t\t);\n\t\t\t}\n\t\t}\n\t\t\n\t\tvar drawWedges = true;\n\t\t\n\t\tfor ( var i = firstHiddenChild; i <= firstChild.hiddenEnd; i++ )\n\t\t{\n\t\t\t// all hidden children must be completely hidden to draw together\n\t\t\t\n\t\t\tif ( this.children[i].alphaPattern.current() != this.children[i].alphaWedge.current() )\n\t\t\t{\n\t\t\t\tdrawWedges = false;\n\t\t\t\tbreak;\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( labelMode )\n\t\t{\n\t\t\tif ( drawWedges )\n\t\t\t{\n\t\t\t\tvar drawRadial = (angleEnd < this.angleEnd.current() + rotationOffset);\n\t\t\t\tthis.drawLines(angleStart, angleEnd, radiusInner, drawRadial);\n\t\t\t}\n\t\t\t\n\t\t\tif ( hiddenSearchResults && ! searchHighlighted )\n\t\t\t{\n\t\t\t\tdrawWedge\n\t\t\t\t(\n\t\t\t\t\tangleStart,\n\t\t\t\t\tangleEnd,\n\t\t\t\t\tradiusInner,\n\t\t\t\t\tgRadius,//this.radiusOuter.current() * gRadius,\n\t\t\t\t\thighlightFill,\n\t\t\t\t\t0,\n\t\t\t\t\ttrue\n\t\t\t\t);\n\t\t\t}\n\t\t}\n\t\telse if ( drawWedges )\n\t\t{\n\t\t\tcontext.globalAlpha = this.alphaWedge.current();\n\t\t\t\n\t\t\tvar fill = rgbText\n\t\t\t(\n\t\t\t\tfirstChild.r.current(),\n\t\t\t\tfirstChild.g.current(),\n\t\t\t\tfirstChild.b.current()\n\t\t\t);\n\t\t\t\n\t\t\tdrawWedge\n\t\t\t(\n\t\t\t\tangleStart,\n\t\t\t\tangleEnd,\n\t\t\t\tradiusInner,\n\t\t\t\tgRadius,//this.radiusOuter.current() * gRadius,\n\t\t\t\tfill,\n\t\t\t\tcontext.globalAlpha,\n\t\t\t\tfalse\n\t\t\t);\n\t\t}\n\t\t\n\t\treturn drawWedges;\n\t}\n\t\n\tthis.drawHiddenLabel = function(angleStart, angleEnd, value, hiddenSearchResults)\n\t{\n\t\tvar textAngle = (angleStart + angleEnd) / 2;\n\t\tvar labelRadius = gRadius + fontSize;//(radiusInner + radius) / 2;\n\t\t\n\t\tvar hiddenLabel = Array();\n\t\t\n\t\thiddenLabel.value = value;\n\t\thiddenLabel.angle = textAngle;\n\t\thiddenLabel.search = hiddenSearchResults;\n\t\t\n\t\tthis.hiddenLabels.push(hiddenLabel);\n\t\t\n\t\tdrawTick(gRadius - fontSize * .75, fontSize * 1.5, textAngle);\n\t\tdrawTextPolar\n\t\t(\n\t\t\tvalue.toString() + \' more\',\n\t\t\t0, // inner text\n\t\t\ttextAngle,\n\t\t\tlabelRadius,\n\t\t\ttrue, // radial\n\t\t\thiddenSearchResults, // bubble\n\t\t\tthis == highlightedNode || this == focusNode, // bold\n\t\t\tfalse,\n\t\t\thiddenSearchResults\n\t\t);\n\t}\n\t\n\tthis.drawHighlight = function(bold)\n\t{\n\t\tvar angleStartCurrent = this.angleStart.current() + rotationOffset;\n\t\tvar angleEndCurrent = this.angleEnd.current() + rotationOffset;\n\t\tvar radiusInner = this.radiusInner.current() * gRadius;\n\t\t\n\t\t//this.setHighlightStyle();\n\t\t\n\t\tif ( this == focusNode && this == highlightedNode && this.hasChildren() )\n\t\t{\n//\t\t\tcontext.fillStyle = "rgba(255, 255, 255, .3)";\n\t\t\tarrow\n\t\t\t(\n\t\t\t\tangleStartCurrent,\n\t\t\t\tangleEndCurrent,\n\t\t\t\tradiusInner\n\t\t\t);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tdrawWedge\n\t\t\t(\n\t\t\t\tangleStartCurrent,\n\t\t\t\tangleEndCurrent,\n\t\t\t\tradiusInner,\n\t\t\t\tgRadius,\n\t\t\t\thighlightFill,\n\t\t\t\t0,\n\t\t\t\ttrue\n\t\t\t);\n\t\t}\n\t\t\n\t\t// check if hidden children should be highlighted\n\t\t//\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tif\n\t\t\t(\n\t\t\t\tthis.children[i].getDepth() - selectedNode.getDepth() + 1 <=\n\t\t\t\tmaxDisplayDepth &&\n\t\t\t\tthis.children[i].hiddenEnd != null\n\t\t\t)\n\t\t\t{\n\t\t\t\tvar firstChild = this.children[i];\n\t\t\t\tvar lastChild = this.children[firstChild.hiddenEnd];\n\t\t\t\tvar hiddenAngleStart = firstChild.angleStart.current() + rotationOffset;\n\t\t\t\tvar hiddenAngleEnd = lastChild.angleEnd.current() + rotationOffset;\n\t\t\t\tvar hiddenRadiusInner = gRadius * firstChild.radiusInner.current();\n\t\t\t\t\n\t\t\t\tdrawWedge\n\t\t\t\t(\n\t\t\t\t\thiddenAngleStart,\n\t\t\t\t\thiddenAngleEnd,\n\t\t\t\t\thiddenRadiusInner,\n\t\t\t\t\tgRadius,\n\t\t\t\t\t\'rgba(255, 255, 255, .3)\',\n\t\t\t\t\t0,\n\t\t\t\t\ttrue\n\t\t\t\t);\n\t\t\t\t\n\t\t\t\tif ( false && ! this.searchResults )\n\t\t\t\t{\n\t\t\t\t\tthis.drawHiddenLabel\n\t\t\t\t\t(\n\t\t\t\t\t\thiddenAngleStart,\n\t\t\t\t\t\thiddenAngleEnd,\n\t\t\t\t\t\tfirstChild.hiddenEnd - i + 1\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\ti = firstChild.hiddenEnd;\n\t\t\t}\n\t\t}\n\t\t\n//\t\t\tcontext.strokeStyle = \'black\';\n\t\tcontext.fillStyle = \'black\';\n\t\t\n\t\tvar highlight = ! ( progress < 1 && zoomOut && this == selectedNodeLast );\n\t\t\n\t\tvar angle = (angleEndCurrent + angleStartCurrent) / 2;\n\t\t\n\t\tif ( ! (this.keyed && showKeys) )\n\t\t{\n\t\t\tthis.drawLabel(angle, true, bold, true, this.radial);\n\t\t}\n\t}\n\t\n\tthis.drawHighlightCenter = function()\n\t{\n\t\tif ( ! this.canDisplayHistory() )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tcontext.lineWidth = highlightLineWidth;\n\t\tcontext.strokeStyle = \'black\';\n\t\tcontext.fillStyle = "rgba(255, 255, 255, .6)";\n\t\t\n\t\tcontext.fillStyle = \'black\';\n\t\tthis.drawLabel(3 * Math.PI / 2, true, true, false);\n\t\tcontext.font = fontNormal;\n\t}\n\t\n\tthis.drawKey = function(angle, highlight, bold)\n\t{\n\t\tvar offset = keyOffset();\n\t\tvar color;\n\t\tvar colorText = this.magnitude == 0 ? \'gray\' : \'black\';\n\t\tvar patternAlpha = this.alphaPattern.end;\n\t\tvar boxLeft = imageWidth - keySize - margin;\n\t\tvar textY = offset + keySize / 2;\n\t\t\n\t\tvar label;\n\t\tvar keyNameWidth;\n\t\t\n\t\tif ( this == selectedNode )\n\t\t{\n\t\t\tcolor = colorUnclassified;\n\t\t\tlabel =\n\t\t\t\tthis.getUnclassifiedText() +\n\t\t\t\t\'   \' +\n\t\t\t\tthis.getUnclassifiedPercentage();\n\t\t\tkeyNameWidth = measureText(label, false);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tlabel = this.keyLabel;\n\t\t\tcolor = rgbText(this.r.end, this.g.end, this.b.end);\n\t\t\t\n\t\t\tif ( highlight )\n\t\t\t{\n\t\t\t\tif ( this.searchResultChildren() )\n\t\t\t\t{\n\t\t\t\t\tlabel = label + searchResultString(this.searchResultChildren());\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tkeyNameWidth = measureText(label, bold);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tkeyNameWidth = this.keyNameWidth;\n\t\t\t}\n\t\t}\n\t\t\n\t\tvar textLeft = boxLeft - keyBuffer - keyNameWidth - fontSize / 2;\n\t\tvar labelLeft = textLeft;\n\t\t\n\t\tif ( labelLeft > keyMinTextLeft - fontSize / 2 )\n\t\t{\n\t\t\tkeyMinTextLeft -= fontSize / 2;\n\t\t\t\n\t\t\tif ( keyMinTextLeft < centerX - gRadius + fontSize / 2 )\n\t\t\t{\n\t\t\t\tkeyMinTextLeft = centerX - gRadius + fontSize / 2;\n\t\t\t}\n\t\t\t\n\t\t\tlabelLeft = keyMinTextLeft;\n\t\t}\n\t\t\n\t\tvar lineX = new Array();\n\t\tvar lineY = new Array();\n\t\t\n\t\tvar bendRadius;\n\t\tvar keyAngle = Math.atan((textY - centerY) / (labelLeft - centerX));\n\t\tvar arcAngle;\n\t\t\n\t\tif ( keyAngle < 0 )\n\t\t{\n\t\t\tkeyAngle += Math.PI;\n\t\t}\n\t\t\n\t\tif ( keyMinAngle == 0 || angle < keyMinAngle )\n\t\t{\n\t\t\tkeyMinAngle = angle;\n\t\t}\n\t\t\n\t\tif ( angle > Math.PI && keyMinAngle > Math.PI )\n\t\t{\n\t\t\t// allow lines to come underneath the chart\n\t\t\t\n\t\t\tangle -= Math.PI * 2;\n\t\t}\n\t\t\n\t\tlineX.push(Math.cos(angle) * gRadius);\n\t\tlineY.push(Math.sin(angle) * gRadius);\n\t\t\n\t\tif ( angle < keyAngle && textY > centerY + Math.sin(angle) * (gRadius + buffer * (currentKey - 1) / (keys + 1) / 2 + buffer / 2) )\n\t\t{\n\t\t\tbendRadius = gRadius + buffer - buffer * currentKey / (keys + 1) / 2;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tbendRadius = gRadius + buffer * currentKey / (keys + 1) / 2 + buffer / 2;\n\t\t}\n\t\t\n\t\tvar outside =\n\t\t\tMath.sqrt\n\t\t\t(\n\t\t\t\tMath.pow(labelLeft - centerX, 2) +\n\t\t\t\tMath.pow(textY - centerY, 2)\n\t\t\t) > bendRadius;\n\t\t\n\t\tif ( ! outside )\n\t\t{\n\t\t\tarcAngle = Math.asin((textY - centerY) / bendRadius);\n\t\t\t\n\t\t\tkeyMinTextLeft = min(keyMinTextLeft, centerX + bendRadius * Math.cos(arcAngle) - fontSize / 2);\n\t\t\t\n\t\t\tif ( labelLeft < textLeft && textLeft > centerX + bendRadius * Math.cos(arcAngle) )\n\t\t\t{\n\t\t\t\tlineX.push(textLeft - centerX);\n\t\t\t\tlineY.push(textY - centerY);\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tkeyMinTextLeft = min(keyMinTextLeft, labelLeft - fontSize / 2);\n\t\t\t\n\t\t\tif ( angle < keyAngle )\n\t\t\t{\n\t\t\t\t// flip everything over y = x\n\t\t\t\t//\n\t\t\t\tarcAngle = Math.PI / 2 - keyLineAngle\n\t\t\t\t(\n\t\t\t\t\tMath.PI / 2 - angle,\n\t\t\t\t\tMath.PI / 2 - keyAngle,\n\t\t\t\t\tbendRadius,\n\t\t\t\t\ttextY - centerY,\n\t\t\t\t\tlabelLeft - centerX,\n\t\t\t\t\tlineY,\n\t\t\t\t\tlineX\n\t\t\t\t);\n\t\t\t\t\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tarcAngle = keyLineAngle\n\t\t\t\t(\n\t\t\t\t\tangle,\n\t\t\t\t\tkeyAngle,\n\t\t\t\t\tbendRadius,\n\t\t\t\t\tlabelLeft - centerX,\n\t\t\t\t\ttextY - centerY,\n\t\t\t\t\tlineX,\n\t\t\t\t\tlineY\n\t\t\t\t);\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( labelLeft > centerX + bendRadius * Math.cos(arcAngle) ||\n\t\ttextY > centerY + bendRadius * Math.sin(arcAngle) + .01)\n//\t\tif ( outside ||  )\n\t\t{\n\t\t\tlineX.push(labelLeft - centerX);\n\t\t\tlineY.push(textY - centerY);\n\t\t\t\n\t\t\tif ( textLeft != labelLeft )\n\t\t\t{\n\t\t\t\tlineX.push(textLeft - centerX);\n\t\t\t\tlineY.push(textY - centerY);\n\t\t\t}\n\t\t}\n\t\t\n\t\tcontext.globalAlpha = this.alphaWedge.current();\n\t\t\n\t\tif ( snapshotMode )\n\t\t{\n\t\t\tvar labelSVG;\n\t\t\t\n\t\t\tif ( this == selectedNode )\n\t\t\t{\n\t\t\t\tlabelSVG =\n\t\t\t\t\tthis.getUnclassifiedText() +\n\t\t\t\t\tspacer() +\n\t\t\t\t\tthis.getUnclassifiedPercentage();\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tlabelSVG = this.name + spacer() + this.getPercentage() + \'%\';\n\t\t\t}\n\t\t\t\n\t\t\tsvg +=\n\t\t\t\t\'<rect fill="\' + color + \'" \' +\n\t\t\t\t\'x="\' + boxLeft + \'" y="\' + offset +\n\t\t\t\t\'" width="\' + keySize + \'" height="\' + keySize + \'"/>\';\n\t\t\t\n\t\t\tif ( patternAlpha )\n\t\t\t{\n\t\t\t\tsvg +=\n\t\t\t\t\t\'<rect fill="url(#hiddenPattern)" style="stroke:none" \' +\n\t\t\t\t\t\'x="\' + boxLeft + \'" y="\' + offset +\n\t\t\t\t\t\'" width="\' + keySize + \'" height="\' + keySize + \'"/>\';\n\t\t\t}\n\t\t\t\n\t\t\tsvg +=\n\t\t\t\t\'<path class="line\' +\n\t\t\t\t(highlight ? \' highlight\' : \'\') +\n\t\t\t\t\'" d="M \' + (lineX[0] + centerX) + \',\' +\n\t\t\t\t(lineY[0] + centerY);\n\t\t\t\n\t\t\tif ( angle != arcAngle )\n\t\t\t{\n\t\t\t\tsvg +=\n\t\t\t\t\t\' L \' + (centerX + bendRadius * Math.cos(angle)) + \',\' +\n\t\t\t\t\t(centerY + bendRadius * Math.sin(angle)) +\n\t\t\t\t\t\' A \' + bendRadius + \',\' + bendRadius + \' 0 \' +\n\t\t\t\t\t\'0,\' + (angle > arcAngle ? \'0\' : \'1\') + \' \' +\n\t\t\t\t\t(centerX + bendRadius * Math.cos(arcAngle)) + \',\' +\n\t\t\t\t\t(centerY + bendRadius * Math.sin(arcAngle));\n\t\t\t}\n\t\t\t\n\t\t\tfor ( var i = 1; i < lineX.length; i++ )\n\t\t\t{\n\t\t\t\tsvg +=\n\t\t\t\t\t\' L \' + (centerX + lineX[i]) + \',\' +\n\t\t\t\t\t(centerY + lineY[i]);\n\t\t\t}\n\t\t\t\n\t\t\tsvg += \'"/>\';\n\t\t\t\n\t\t\tif ( highlight )\n\t\t\t{\n\t\t\t\tif ( this.searchResultChildren() )\n\t\t\t\t{\n\t\t\t\t\tlabelSVG = labelSVG + searchResultString(this.searchResultChildren());\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tdrawBubbleSVG\n\t\t\t\t(\n\t\t\t\t\tboxLeft - keyBuffer - keyNameWidth - fontSize / 2,\n\t\t\t\t\ttextY - fontSize,\n\t\t\t\t\tkeyNameWidth + fontSize,\n\t\t\t\t\tfontSize * 2,\n\t\t\t\t\tfontSize,\n\t\t\t\t\t0\n\t\t\t\t);\n\t\t\t\t\n\t\t\t\tif ( this.isSearchResult )\n\t\t\t\t{\n\t\t\t\t\tdrawSearchHighlights\n\t\t\t\t\t(\n\t\t\t\t\t\tlabel,\n\t\t\t\t\t\tboxLeft - keyBuffer - keyNameWidth,\n\t\t\t\t\t\ttextY,\n\t\t\t\t\t\t0\n\t\t\t\t\t)\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\tsvg += svgText(labelSVG, boxLeft - keyBuffer, textY, \'end\', bold, colorText);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.fillStyle = color;\n\t\t\tcontext.translate(-centerX, -centerY);\n\t\t\tcontext.strokeStyle = \'black\';\n\t\t\t\tcontext.globalAlpha = 1;//this.alphaWedge.current();\n\t\t\t\n\t\t\tcontext.fillRect(boxLeft, offset, keySize, keySize);\n\t\t\t\n\t\t\tif ( patternAlpha )\n\t\t\t{\n\t\t\t\tcontext.globalAlpha = patternAlpha;\n\t\t\t\tcontext.fillStyle = hiddenPattern;\n\t\t\t\t\n\t\t\t\t// make clipping box for Firefox performance\n\t\t\t\tcontext.beginPath();\n\t\t\t\tcontext.moveTo(boxLeft, offset);\n\t\t\t\tcontext.lineTo(boxLeft + keySize, offset);\n\t\t\t\tcontext.lineTo(boxLeft + keySize, offset + keySize);\n\t\t\t\tcontext.lineTo(boxLeft, offset + keySize);\n\t\t\t\tcontext.closePath();\n\t\t\t\tcontext.save();\n\t\t\t\tcontext.clip();\n\t\t\t\t\n\t\t\t\tcontext.fillRect(boxLeft, offset, keySize, keySize);\n\t\t\t\tcontext.fillRect(boxLeft, offset, keySize, keySize);\n\t\t\t\t\n\t\t\t\tcontext.restore(); // remove clipping region\n\t\t\t}\n\t\t\t\n\t\t\tif ( highlight )\n\t\t\t{\n\t\t\t\tthis.setHighlightStyle();\n\t\t\t\tcontext.fillRect(boxLeft, offset, keySize, keySize);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tcontext.lineWidth = thinLineWidth;\n\t\t\t}\n\t\t\t\n\t\t\tcontext.strokeRect(boxLeft, offset, keySize, keySize);\n\t\t\t\n\t\t\tif ( lineX.length )\n\t\t\t{\n\t\t\t\tcontext.beginPath();\n\t\t\t\tcontext.moveTo(lineX[0] + centerX, lineY[0] + centerY);\n\t\t\t\t\n\t\t\t\tcontext.arc(centerX, centerY, bendRadius, angle, arcAngle, angle > arcAngle);\n\t\t\t\t\n\t\t\t\tfor ( var i = 1; i < lineX.length; i++ )\n\t\t\t\t{\n\t\t\t\t\tcontext.lineTo(lineX[i] + centerX, lineY[i] + centerY);\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tcontext.globalAlpha = this == selectedNode ?\n\t\t\t\t\tthis.children[0].alphaWedge.current() :\n\t\t\t\t\tthis.alphaWedge.current();\n\t\t\t\tcontext.lineWidth = highlight ? highlightLineWidth : thinLineWidth;\n\t\t\t\tcontext.stroke();\n\t\t\t\tcontext.globalAlpha = 1;\n\t\t\t}\n\t\t\t\n\t\t\tif ( highlight )\n\t\t\t{\n\t\t\t\tdrawBubbleCanvas\n\t\t\t\t(\n\t\t\t\t\tboxLeft - keyBuffer - keyNameWidth - fontSize / 2,\n\t\t\t\t\ttextY - fontSize,\n\t\t\t\t\tkeyNameWidth + fontSize,\n\t\t\t\t\tfontSize * 2,\n\t\t\t\t\tfontSize,\n\t\t\t\t\t0\n\t\t\t\t);\n\t\t\t\t\n\t\t\t\tif ( this.isSearchResult )\n\t\t\t\t{\n\t\t\t\t\tdrawSearchHighlights\n\t\t\t\t\t(\n\t\t\t\t\t\tlabel,\n\t\t\t\t\t\tboxLeft - keyBuffer - keyNameWidth,\n\t\t\t\t\t\ttextY,\n\t\t\t\t\t\t0\n\t\t\t\t\t)\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\tdrawText(label, boxLeft - keyBuffer, offset + keySize / 2, 0, \'end\', bold, colorText);\n\t\t\t\n\t\t\tcontext.translate(centerX, centerY);\n\t\t}\n\t\t\n\t\tcurrentKey++;\n\t}\n\t\n\tthis.drawLabel = function(angle, bubble, bold, selected, radial)\n\t{\n\t\tif ( context.globalAlpha == 0 )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar innerText;\n\t\tvar label;\n\t\tvar radius;\n\t\t\n\t\tif ( radial )\n\t\t{\n\t\t\tradius = (this.radiusInner.current() + 1) * gRadius / 2;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tradius = this.labelRadius.current() * gRadius;\n\t\t}\n\t\t\n\t\tif ( radial && (selected || bubble ) )\n\t\t{\n\t\t\tvar percentage = this.getPercentage();\n\t\t\tinnerText = percentage + \'%\';\n\t\t}\n\t\t\n\t\tif\n\t\t(\n\t\t\t! radial &&\n\t\t\tthis != selectedNode &&\n\t\t\t! bubble &&\n\t\t\t( !zoomOut || this != selectedNodeLast)\n\t\t)\n\t\t{\n\t\t\tlabel = this.shortenLabel();\n\t\t}\n\t\telse\n\t\t{\n\t\t\tlabel = this.name;\n\t\t}\n\t\t\n\t\tvar flipped = drawTextPolar\n\t\t(\n\t\t\tlabel,\n\t\t\tinnerText,\n\t\t\tangle,\n\t\t\tradius,\n\t\t\tradial,\n\t\t\tbubble,\n\t\t\tbold,\n//\t\t\tthis.isSearchResult && this.shouldAddSearchResultsString() && (!selected || this == selectedNode || highlight),\n\t\t\tthis.isSearchResult && (!selected || this == selectedNode || bubble),\n\t\t\t(this.hideAlone || !selected || this == selectedNode ) ? this.searchResultChildren() : 0\n\t\t);\n\t\t\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n\t\t\n\t\tif\n\t\t(\n\t\t\t! radial &&\n\t\t\t! bubble &&\n\t\t\tthis != selectedNode &&\n\t\t\tthis.angleEnd.end != this.angleStart.end &&\n\t\t\tnLabelOffsets[depth - 2] > 2 &&\n\t\t\tthis.labelWidth.current() > (this.angleEnd.end - this.angleStart.end) * Math.abs(radius) &&\n\t\t\t! ( zoomOut && this == selectedNodeLast ) &&\n\t\t\tthis.labelRadius.end > 0\n\t\t)\n\t\t{\n\t\t\t// name extends beyond wedge; draw tick mark towards the central\n\t\t\t// radius for easier identification\n\t\t\t\n\t\t\tvar radiusCenter = compress ?\n\t\t\t\t(compressedRadii[depth - 1] + compressedRadii[depth - 2]) / 2 :\n\t\t\t\t(depth - .5) * nodeRadius;\n\t\t\t\n\t\t\tif ( this.labelRadius.end > radiusCenter )\n\t\t\t{\n\t\t\t\tif ( flipped )\n\t\t\t\t{\n\t\t\t\t\tdrawTick(radius - tickLength * 1.4 , tickLength, angle);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tdrawTick(radius - tickLength * 1.7, tickLength, angle);\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tif ( flipped )\n\t\t\t\t{\n\t\t\t\t\tdrawTick(radius + tickLength * .7, tickLength, angle);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tdrawTick(radius + tickLength * .4, tickLength, angle);\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t}\n\t\n\tthis.drawLines = function(angleStart, angleEnd, radiusInner, drawRadial, selected)\n\t{\n\t\tif ( snapshotMode )\n\t\t{\n\t\t\tif ( this != selectedNode)\n\t\t\t{\n\t\t\t\tif ( angleEnd == angleStart + Math.PI * 2 )\n\t\t\t\t{\n\t\t\t\t\t// fudge to prevent overlap, which causes arc ambiguity\n\t\t\t\t\t//\n\t\t\t\t\tangleEnd -= .1 / gRadius;\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tvar longArc = angleEnd - angleStart > Math.PI ? 1 : 0;\n\t\t\t\t\n\t\t\t\tvar x1 = centerX + radiusInner * Math.cos(angleStart);\n\t\t\t\tvar y1 = centerY + radiusInner * Math.sin(angleStart);\n\t\t\t\t\n\t\t\t\tvar x2 = centerX + gRadius * Math.cos(angleStart);\n\t\t\t\tvar y2 = centerY + gRadius * Math.sin(angleStart);\n\t\t\t\t\n\t\t\t\tvar x3 = centerX + gRadius * Math.cos(angleEnd);\n\t\t\t\tvar y3 = centerY + gRadius * Math.sin(angleEnd);\n\t\t\t\t\n\t\t\t\tvar x4 = centerX + radiusInner * Math.cos(angleEnd);\n\t\t\t\tvar y4 = centerY + radiusInner * Math.sin(angleEnd);\n\t\t\t\t\n\t\t\t\tif ( this.alphaArc.end )\n\t\t\t\t{\n\t\t\t\t\tvar dArray =\n\t\t\t\t\t[\n\t\t\t\t\t\t" M ", x4, ",", y4,\n\t\t\t\t\t\t" A ", radiusInner, ",", radiusInner, " 0 ", longArc,\n\t\t\t\t\t\t\t" 0 ", x1, ",", y1\n\t\t\t\t\t];\n\t\t\t\t\t\n\t\t\t\t\tsvg += \'<path class="line" d="\' + dArray.join(\'\') + \'"/>\';\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( drawRadial && this.alphaLine.end )\n\t\t\t\t{\n\t\t\t\t\tsvg += \'<line x1="\' + x3 + \'" y1="\' + y3 + \'" x2="\' + x4 + \'" y2="\' + y4 + \'"/>\';\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.lineWidth = thinLineWidth;\n\t\t\tcontext.strokeStyle = \'black\';\n\t\t\tcontext.beginPath();\n\t\t\tcontext.arc(0, 0, radiusInner, angleStart, angleEnd, false);\n\t\t\tcontext.globalAlpha = this.alphaArc.current();\n\t\t\tcontext.stroke();\n\t\t\t\n\t\t\tif ( drawRadial )\n\t\t\t{\n\t\t\t\tvar x1 = radiusInner * Math.cos(angleEnd);\n\t\t\t\tvar y1 = radiusInner * Math.sin(angleEnd);\n\t\t\t\tvar x2 = gRadius * Math.cos(angleEnd);\n\t\t\t\tvar y2 = gRadius * Math.sin(angleEnd);\n\t\t\t\t\n\t\t\t\tcontext.beginPath();\n\t\t\t\tcontext.moveTo(x1, y1);\n\t\t\t\tcontext.lineTo(x2, y2);\n\t\t\t\t\n//\t\t\t\tif ( this.getCollapse() )//( selected && this != selectedNode )\n\t\t\t\t{\n\t\t\t\t\tcontext.globalAlpha = this.alphaLine.current();\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tcontext.stroke();\n\t\t\t}\n\t\t}\n\t}\n\t\n\tthis.drawMap = function(child)\n\t{\n\t\tif ( this.parent )\n\t\t{\n\t\t\tthis.parent.drawMap(child);\n\t\t}\n\t\t\n\t\tif ( this.getCollapse() && this != child || this == focusNode )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar angleStart =\n\t\t\t(child.baseMagnitude - this.baseMagnitude) / this.magnitude * Math.PI * 2 +\n\t\t\trotationOffset;\n\t\tvar angleEnd =\n\t\t\t(child.baseMagnitude - this.baseMagnitude + child.magnitude) /\n\t\t\tthis.magnitude * Math.PI * 2 +\n\t\t\trotationOffset;\n\t\t\n\t\tvar box = this.getMapPosition();\n\t\t\n\t\tcontext.save();\n\t\tcontext.fillStyle = \'black\';\n\t\tcontext.textAlign = \'end\';\n\t\tcontext.textBaseline = \'middle\';\n\t\t\n\t\tvar textX = box.x - mapRadius - mapBuffer;\n\t\tvar percentage = getPercentage(child.magnitude / this.magnitude);\n\t\t\n\t\tvar highlight = this == selectedNode || this == highlightedNode;\n\t\t\n\t\tif ( highlight )\n\t\t{\n\t\t\tcontext.font = fontBold;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.font = fontNormal;\n\t\t}\n\t\t\n\t\tcontext.fillText(percentage + \'% of\', textX, box.y - mapRadius / 3);\n\t\tcontext.fillText(this.name, textX, box.y + mapRadius / 3);\n\t\t\n\t\tif ( highlight )\n\t\t{\n\t\t\tcontext.font = fontNormal;\n\t\t}\n\t\t\n\t\tif ( this == highlightedNode && this != selectedNode )\n\t\t{\n\t\t\tcontext.fillStyle = \'rgb(245, 245, 245)\';\n//\t\t\tcontext.fillStyle = \'rgb(200, 200, 200)\';\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.fillStyle = \'rgb(255, 255, 255)\';\n\t\t}\n\t\t\n\t\tcontext.beginPath();\n\t\tcontext.arc(box.x, box.y, mapRadius, 0, Math.PI * 2, true);\n\t\tcontext.closePath();\n\t\tcontext.fill();\n\t\t\n\t\tif ( this == selectedNode )\n\t\t{\n\t\t\tcontext.lineWidth = 1;\n\t\t\tcontext.fillStyle = \'rgb(100, 100, 100)\';\n\t\t}\n\t\telse\n\t\t{\n\t\t\tif ( this == highlightedNode )\n\t\t\t{\n\t\t\t\tcontext.lineWidth = .2;\n\t\t\t\tcontext.fillStyle = \'rgb(190, 190, 190)\';\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tcontext.lineWidth = .2;\n\t\t\t\tcontext.fillStyle = \'rgb(200, 200, 200)\';\n\t\t\t}\n\t\t}\n\t\t\n\t\tvar maxDepth = this.getMaxDepth();\n\t\t\n\t\tif ( ! compress && maxDepth > maxPossibleDepth + this.getDepth() - 1 )\n\t\t{\n\t\t\tmaxDepth = maxPossibleDepth + this.getDepth() - 1;\n\t\t}\n\t\t\n\t\tif ( this.getDepth() < selectedNode.getDepth() )\n\t\t{\n\t\t\tif ( child.getDepth() - 1 >= maxDepth )\n\t\t\t{\n\t\t\t\tmaxDepth = child.getDepth();\n\t\t\t}\n\t\t}\n\t\t\n\t\tvar radiusInner;\n\t\t\n\t\tif ( compress )\n\t\t{\n\t\t\tradiusInner = 0;\n//\t\t\t\tMath.atan(child.getDepth() - this.getDepth()) /\n//\t\t\t\tMath.PI * 2 * .9;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tradiusInner =\n\t\t\t\t(child.getDepth() - this.getDepth()) /\n\t\t\t\t(maxDepth - this.getDepth() + 1);\n\t\t}\n\t\t\n\t\tcontext.stroke();\n\t\tcontext.beginPath();\n\t\t\n\t\tif ( radiusInner == 0 )\n\t\t{\n\t\t\tcontext.moveTo(box.x, box.y);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.arc(box.x, box.y, mapRadius * radiusInner, angleEnd, angleStart, true);\n\t\t}\n\t\t\n\t\tcontext.arc(box.x, box.y, mapRadius, angleStart, angleEnd, false);\n\t\tcontext.closePath();\n\t\tcontext.fill();\n\t\t\n\t\tif ( this == highlightedNode && this != selectedNode )\n\t\t{\n\t\t\tcontext.lineWidth = 1;\n\t\t\tcontext.stroke();\n\t\t}\n\t\t\n\t\tcontext.restore();\n\t}\n\t\n\tthis.drawReferenceRings = function(childRadiusInner)\n\t{\n\t\tif ( snapshotMode )\n\t\t{\n\t\t\tsvg +=\n\t\t\t\t\'<circle cx="\' + centerX + \'" cy="\' + centerY +\n\t\t\t\t\'" r="\' + childRadiusInner + \'"/>\';\n\t\t\tsvg +=\n\t\t\t\t\'<circle cx="\' + centerX + \'" cy="\' + centerY +\n\t\t\t\t\'" r="\' + gRadius + \'"/>\';\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.globalAlpha = 1 - this.alphaLine.current();//this.getUncollapsed().alphaLine.current();\n\t\t\tcontext.beginPath();\n\t\t\tcontext.arc(0, 0, childRadiusInner, 0, Math.PI * 2, false);\n\t\t\tcontext.stroke();\n\t\t\tcontext.beginPath();\n\t\t\tcontext.arc(0, 0, gRadius, 0, Math.PI * 2, false);\n\t\t\tcontext.stroke();\n\t\t}\n\t}\n\t\n\tthis.getCollapse = function()\n\t{\n\t\treturn (\n\t\t\tcollapse &&\n\t\t\tthis.collapse &&\n\t\t\tthis.depth != maxAbsoluteDepth\n\t\t\t);\n\t}\n\t\n\tthis.getDepth = function()\n\t{\n\t\tif ( collapse )\n\t\t{\n\t\t\treturn this.depthCollapsed;\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this.depth;\n\t\t}\n\t}\n\t\n\tthis.getMagnitude = function()\n\t{\n\t\treturn this.attributes[magnitudeIndex][currentDataset];\n\t}\n\t\n\tthis.getMapPosition = function()\n\t{\n\t\treturn {\n\t\t\tx : (details.offsetLeft + details.clientWidth - mapRadius),\n\t\t\ty : ((focusNode.getDepth() - this.getDepth()) *\n\t\t\t\t(mapBuffer + mapRadius * 2) - mapRadius) +\n\t\t\t\tdetails.clientHeight + details.offsetTop\n\t\t};\n\t}\n\t\n\tthis.getMaxDepth = function(limit)\n\t{\n\t\tvar max;\n\t\t\n\t\tif ( collapse )\n\t\t{\n\t\t\treturn this.maxDepthCollapsed;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tif ( this.maxDepth > maxAbsoluteDepth )\n\t\t\t{\n\t\t\t\treturn maxAbsoluteDepth;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\treturn this.maxDepth;\n\t\t\t}\n\t\t}\n\t}\n\t\n\tthis.getData = function(index, summary)\n\t{\n\t\tvar files = new Array();\n\t\t\n\t\tif\n\t\t(\n\t\t\tthis.attributes[index] != null &&\n\t\t\tthis.attributes[index][currentDataset] != null &&\n\t\t\tthis.attributes[index][currentDataset] != \'\'\n\t\t)\n\t\t{\n\t\t\tfiles.push\n\t\t\t(\n\t\t\t\tdocument.location +\n\t\t\t\t\'.files/\' +\n\t\t\t\tthis.attributes[index][currentDataset]\n\t\t\t);\n\t\t}\n\t\t\n\t\tif ( summary )\n\t\t{\n\t\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t\t{\n\t\t\t\tfiles = files.concat(this.children[i].getData(index, true));\n\t\t\t}\n\t\t}\n\t\t\n\t\treturn files;\n\t}\n\t\n\tthis.getList = function(index, summary)\n\t{\n\t\tvar list;\n\t\t\n\t\tif\n\t\t(\n\t\t\tthis.attributes[index] != null &&\n\t\t\tthis.attributes[index][currentDataset] != null\n\t\t)\n\t\t{\n\t\t\tlist = this.attributes[index][currentDataset];\n\t\t}\n\t\telse\n\t\t{\n\t\t\tlist = new Array();\n\t\t}\n\t\t\n\t\tif ( summary )\n\t\t{\n\t\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t\t{\n\t\t\t\tlist = list.concat(this.children[i].getList(index, true));\n\t\t\t}\n\t\t}\n\t\t\n\t\treturn list;\n\t}\n\t\n\tthis.getParent = function()\n\t{\n\t\t// returns parent, accounting for collapsing or 0 if doesn\'t exist\n\t\t\n\t\tvar parent = this.parent;\n\t\t\n\t\twhile ( parent != 0 && parent.getCollapse() )\n\t\t{\n\t\t\tparent = parent.parent;\n\t\t}\n\t\t\n\t\treturn parent;\n\t}\n\t\n\tthis.getPercentage = function()\n\t{\n\t\treturn getPercentage(this.magnitude / selectedNode.magnitude);\n\t}\n\t\n\tthis.getUnclassifiedPercentage = function()\n\t{\n\t\tif ( this.children.length )\n\t\t{\n\t\t\tvar lastChild = this.children[this.children.length - 1];\n\t\t\n\t\t\treturn getPercentage\n\t\t\t(\n\t\t\t\t(\n\t\t\t\t\tthis.baseMagnitude +\n\t\t\t\t\tthis.magnitude -\n\t\t\t\t\tlastChild.magnitude -\n\t\t\t\t\tlastChild.baseMagnitude\n\t\t\t\t) / this.magnitude\n\t\t\t) + \'%\';\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn \'100%\';\n\t\t}\n\t}\n\t\n\tthis.getUnclassifiedText = function()\n\t{\n\t\treturn \'[other \'+ this.name + \']\';\n\t}\n\t\n\tthis.getUncollapsed = function()\n\t{\n\t\t// recurse through collapsed children until uncollapsed node is found\n\t\t\n\t\tif ( this.getCollapse() )\n\t\t{\n\t\t\treturn this.children[0].getUncollapsed();\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this;\n\t\t}\n\t}\n\t\n\tthis.hasChildren = function()\n\t{\n\t\treturn this.children.length && this.depth < maxAbsoluteDepth && this.magnitude;\n\t}\n\t\n\tthis.hasParent = function(parent)\n\t{\n\t\tif ( this.parent )\n\t\t{\n\t\t\tif ( this.parent == parent )\n\t\t\t{\n\t\t\t\treturn true;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\treturn this.parent.hasParent(parent);\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn false;\n\t\t}\n\t}\n\t\n\tthis.maxVisibleDepth = function(maxDepth)\n\t{\n\t\tvar childInnerRadius;\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n\t\tvar currentMaxDepth = depth;\n\t\t\n\t\tif ( this.hasChildren() && depth < maxDepth)\n\t\t{\n\t\t\tvar lastChild = this.children[this.children.length - 1];\n\t\t\t\n\t\t\tif ( this.name == \'Pseudomonadaceae\' )\n\t\t\t{\n\t\t\t\tvar x = 3;\n\t\t\t}\n\t\t\t\n\t\t\tif\n\t\t\t(\n\t\t\t\tlastChild.baseMagnitude + lastChild.magnitude <\n\t\t\t\tthis.baseMagnitude + this.magnitude\n\t\t\t)\n\t\t\t{\n\t\t\t\tcurrentMaxDepth++;\n\t\t\t}\n\t\t\t\n\t\t\tif ( compress )\n\t\t\t{\n\t\t\t\tchildInnerRadius = compressedRadii[depth - 1];\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tchildInnerRadius = (depth) / maxDepth;\n\t\t\t}\n\t\t\t\n\t\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t\t{\n\t\t\t\tif\n\t\t\t\t(//true ||\n\t\t\t\t\tthis.children[i].magnitude *\n\t\t\t\t\tangleFactor *\n\t\t\t\t\t(childInnerRadius + 1) *\n\t\t\t\t\tgRadius >=\n\t\t\t\t\tminWidth()\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tvar childMaxDepth = this.children[i].maxVisibleDepth(maxDepth);\n\t\t\t\t\t\n\t\t\t\t\tif ( childMaxDepth > currentMaxDepth )\n\t\t\t\t\t{\n\t\t\t\t\t\tcurrentMaxDepth = childMaxDepth;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n\t\treturn currentMaxDepth;\n\t}\n\t\n\tthis.resetLabelWidth = function()\n\t{\n\t\tvar nameWidthOld = this.nameWidth;\n\t\t\n\t\tif ( true || ! this.radial )//&& fontSize != fontSizeLast )\n\t\t{\n\t\t\tvar dim = context.measureText(this.name);\n\t\t\tthis.nameWidth = dim.width;\n\t\t}\n\t\t\n\t\tif ( fontSize != fontSizeLast && this.labelWidth.end == nameWidthOld * labelWidthFudge )\n\t\t{\n\t\t\t// font size changed; adjust start of tween to match\n\t\t\t\n\t\t\tthis.labelWidth.start = this.nameWidth * labelWidthFudge;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.labelWidth.start = this.labelWidth.current();\n\t\t}\n\t\t\n\t\tthis.labelWidth.end = this.nameWidth * labelWidthFudge;\n\t}\n\t\n\tthis.restrictLabelWidth = function(width)\n\t{\n\t\tif ( width < this.labelWidth.end )\n\t\t{\n\t\t\tthis.labelWidth.end = width;\n\t\t}\n\t}\n\t\n\tthis.search = function()\n\t{\n\t\tthis.isSearchResult = false;\n\t\tthis.searchResults = 0;\n\t\t\n\t\tif\n\t\t(\n\t\t\t! this.getCollapse() &&\n\t\t\tsearch.value != \'\' &&\n\t\t\tthis.name.toLowerCase().indexOf(search.value.toLowerCase()) != -1\n\t\t)\n\t\t{\n\t\t\tthis.isSearchResult = true;\n\t\t\tthis.searchResults = 1;\n\t\t\tnSearchResults++;\n\t\t}\n\t\t\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.searchResults += this.children[i].search();\n\t\t}\n\t\t\n\t\treturn this.searchResults;\n\t}\n\t\n\tthis.searchResultChildren = function()\n\t{\n\t\tif ( this.isSearchResult )\n\t\t{\n\t\t\treturn this.searchResults - 1;\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this.searchResults;\n\t\t}\n\t}\n\t\n\tthis.setDepth = function(depth, depthCollapsed)\n\t{\n\t\tthis.depth = depth;\n\t\tthis.depthCollapsed = depthCollapsed;\n\t\t\n\t\tif\n\t\t(\n\t\t\tthis.children.length == 1 &&\n//\t\t\tthis.magnitude > 0 &&\n\t\t\tthis.children[0].magnitude == this.magnitude &&\n\t\t\t( head.children.length > 1 || this.children[0].children.length )\n\t\t)\n\t\t{\n\t\t\tthis.collapse = true;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.collapse = false;\n\t\t\tdepthCollapsed++;\n\t\t}\n\t\t\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.children[i].setDepth(depth + 1, depthCollapsed);\n\t\t}\n\t}\n\t\n\tthis.setHighlightStyle = function()\n\t{\n\t\tcontext.lineWidth = highlightLineWidth;\n\t\t\n\t\tif ( this.hasChildren() || this != focusNode || this != highlightedNode )\n\t\t{\n\t\t\tcontext.strokeStyle = \'black\';\n\t\t\tcontext.fillStyle = "rgba(255, 255, 255, .3)";\n\t\t}\n\t\telse\n\t\t{\n\t\t\tcontext.strokeStyle = \'rgb(90,90,90)\';\n\t\t\tcontext.fillStyle = "rgba(155, 155, 155, .3)";\n\t\t}\n\t}\n\t\n\tthis.setLabelWidth = function(node)\n\t{\n\t\tif ( ! shorten || this.radial )\n\t\t{\n\t\t\treturn; // don\'t need to set width\n\t\t}\n\t\t\n\t\tif ( node.hide )\n\t\t{\n\t\t\talert(\'wtf\');\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar angle = (this.angleStart.end + this.angleEnd.end) / 2;\n\t\tvar a; // angle difference\n\t\t\n\t\tif ( node == selectedNode )\n\t\t{\n\t\t\ta = Math.abs(angle - node.angleOther);\n\t\t}\n\t\telse\n\t\t{\n\t\t\ta = Math.abs(angle - (node.angleStart.end + node.angleEnd.end) / 2);\n\t\t}\n\t\t\n\t\tif ( a == 0 )\n\t\t{\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tif ( a > Math.PI )\n\t\t{\n\t\t\ta = 2 * Math.PI - a;\n\t\t}\n\t\t\n\t\tif ( node.radial || node == selectedNode )\n\t\t{\n\t\t\tvar nodeLabelRadius;\n\t\t\t\n\t\t\tif ( node == selectedNode )\n\t\t\t{\n\t\t\t\t// radial \'other\' label\n\t\t\t\t\n\t\t\t\tnodeLabelRadius = (node.children[0].radiusInner.end + 1) / 2;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tnodeLabelRadius = (node.radiusInner.end + 1) / 2;\n\t\t\t}\n\t\t\t\n\t\t\tif ( a < Math.PI / 2 )\n\t\t\t{\n\t\t\t\tvar r = this.labelRadius.end * gRadius + .5 * fontSize\n\t\t\t\tvar hypotenuse = r / Math.cos(a);\n\t\t\t\tvar opposite = r * Math.tan(a);\n\t\t\t\tvar fontRadius = .8 * fontSize;\n\t\t\t\t\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\tnodeLabelRadius * gRadius < hypotenuse &&\n\t\t\t\t\tthis.labelWidth.end / 2 + fontRadius > opposite\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tthis.labelWidth.end = 2 * (opposite - fontRadius);\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\telse if\n\t\t(\n\t\t\tthis.labelRadius.end == node.labelRadius.end &&\n\t\t\ta < Math.PI / 4\n\t\t)\n\t\t{\n\t\t\t// same radius with small angle; use circumferential approximation\n\t\t\t\n\t\t\tvar dist = a * this.labelRadius.end * gRadius - fontSize * (1 - a * 4 / Math.PI) * 1.3;\n\t\t\t\n\t\t\tif ( this.labelWidth.end < dist )\n\t\t\t{\n\t\t\t\tnode.restrictLabelWidth((dist - this.labelWidth.end / 2) * 2);\n\t\t\t}\n\t\t\telse if ( node.labelWidth.end < dist )\n\t\t\t{\n\t\t\t\tthis.restrictLabelWidth((dist - node.labelWidth.end / 2) * 2);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\t// both labels reach halfway point; restrict both\n\t\t\t\t\n\t\t\t\tthis.labelWidth.end = dist;\n\t\t\t\tnode.labelWidth.end = dist\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tvar r1 = this.labelRadius.end * gRadius;\n\t\t\tvar r2 = node.labelRadius.end * gRadius;\n\t\t\t\n\t\t\t// first adjust the radii to account for the height of the font by shifting them\n\t\t\t// toward each other\n\t\t\t//\n\t\t\tvar fontFudge = .35 * fontSize;\n\t\t\t//\n\t\t\tif ( this.labelRadius.end < node.labelRadius.end )\n\t\t\t{\n\t\t\t\tr1 += fontFudge;\n\t\t\t\tr2 -= fontFudge;\n\t\t\t}\n\t\t\telse if ( this.labelRadius.end > node.labelRadius.end )\n\t\t\t{\n\t\t\t\tr1 -= fontFudge;\n\t\t\t\tr2 += fontFudge;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tr1 -= fontFudge;\n\t\t\t\tr2 -= fontFudge;\n\t\t\t}\n\t\t\t\n\t\t\tvar r1s = r1 * r1;\n\t\t\tvar r2s = r2 * r2;\n\t\t\t\n\t\t\t// distance between the centers of the two labels\n\t\t\t//\n\t\t\tvar dist = Math.sqrt(r1s + r2s - 2 * r1 * r2 * Math.cos(a));\n\t\t\t\n\t\t\t// angle at our label center between our radius and the line to the other label center\n\t\t\t//\n\t\t\tvar b = Math.acos((r1s + dist * dist - r2s) / (2 * r1 * dist));\n\t\t\t\n\t\t\t// distance from our label center to the intersection of the two tangents\n\t\t\t//\n\t\t\tvar l1 = Math.sin(a + b - Math.PI / 2) * dist / Math.sin(Math.PI - a);\n\t\t\t\n\t\t\t// distance from other label center the the intersection of the two tangents\n\t\t\t//\n\t\t\tvar l2 = Math.sin(Math.PI / 2 - b) * dist / Math.sin(Math.PI - a);\n\t\t\t\n\t\t\tl1 = Math.abs(l1) - .4 * fontSize;\n\t\t\tl2 = Math.abs(l2) - .4 * fontSize;\n/*\t\t\t\n\t\t\t// amount to shorten the distances because of the height of the font\n\t\t\t//\n\t\t\tvar l3 = 0;\n\t\t\tvar fontRadius = fontSize * .55;\n\t\t\t//\n\t\t\tif ( l1 < 0 || l2 < 0 )\n\t\t\t{\n\t\t\t\tvar l4 = fontRadius / Math.tan(a);\n\t\t\tl1 = Math.abs(l1);\n\t\t\tl2 = Math.abs(l2);\n\t\t\t\n\t\t\t\tl1 -= l4;\n\t\t\t\tl2 -= l4;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tvar c = Math.PI - a;\n\t\t\t\t\n\t\t\t\tl3 = fontRadius * Math.tan(c / 2);\n\t\t\t}\n*/\t\t\t\n\t\t\tif ( this.labelWidth.end / 2 > l1 && node.labelWidth.end / 2 > l2 )\n\t\t\t{\n\t\t\t\t// shorten the farthest one from the intersection\n\t\t\t\t\n\t\t\t\tif ( l1 > l2 )\n\t\t\t\t{\n\t\t\t\t\tthis.restrictLabelWidth(2 * (l1));// - l3 - fontRadius));\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tnode.restrictLabelWidth(2 * (l2));// - l3 - fontRadius));\n\t\t\t\t}\n\t\t\t}/*\n\t\t\telse if ( this.labelWidth.end / 2 > l1 + l3 && node.labelWidth.end / 2 > l2 - l3 )\n\t\t\t{\n\t\t\t\tnode.restrictLabelWidth(2 * (l2 - l3));\n\t\t\t}\n\t\t\telse if ( this.labelWidth.end / 2 > l1 - l3 && node.labelWidth.end / 2 > l2 + l3 )\n\t\t\t{\n\t\t\t\tthis.restrictLabelWidth(2 * (l1 - l3));\n\t\t\t}*/\n\t\t}\n\t}\n\t\n\tthis.setMagnitudes = function(baseMagnitude)\n\t{\n\t\tthis.magnitude = this.getMagnitude();\n\t\tthis.baseMagnitude = baseMagnitude;\n\t\t\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.children[i].setMagnitudes(baseMagnitude);\n\t\t\tbaseMagnitude += this.children[i].magnitude;\n\t\t}\n\t\t\n\t\tthis.maxChildMagnitude = baseMagnitude;\n\t}\n\t\n\tthis.setMaxDepths = function()\n\t{\n\t\tthis.maxDepth = this.depth;\n\t\tthis.maxDepthCollapsed = this.depthCollapsed;\n\t\t\n\t\tfor ( i in this.children )\n\t\t{\n\t\t\tvar child = this.children[i];\n\t\t\t\n\t\t\tchild.setMaxDepths();\n\t\t\t\n\t\t\tif ( child.maxDepth > this.maxDepth )\n\t\t\t{\n\t\t\t\tthis.maxDepth = child.maxDepth;\n\t\t\t}\n\t\t\t\n\t\t\tif\n\t\t\t(\n\t\t\t\tchild.maxDepthCollapsed > this.maxDepthCollapsed &&\n\t\t\t\t(child.depth <= maxAbsoluteDepth || maxAbsoluteDepth == 0)\n\t\t\t)\n\t\t\t{\n\t\t\t\tthis.maxDepthCollapsed = child.maxDepthCollapsed;\n\t\t\t}\n\t\t}\n\t}\n\t\n\tthis.setTargetLabelRadius = function()\n\t{\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n\t\tvar index = depth - 2;\n\t\tvar labelOffset = labelOffsets[index];\n\t\t\n\t\tif ( this.radial )\n\t\t{\n\t\t\t//this.labelRadius.setTarget((this.radiusInner.end + 1) / 2);\n\t\t\tvar max =\n\t\t\t\tdepth == maxDisplayDepth ?\n\t\t\t\t1 :\n\t\t\t\tcompressedRadii[index + 1];\n\t\t\t\n\t\t\tthis.labelRadius.setTarget((compressedRadii[index] + max) / 2);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tvar radiusCenter;\n\t\t\tvar width;\n\t\t\t\n\t\t\tif ( compress )\n\t\t\t{\n\t\t\t\tif ( nLabelOffsets[index] > 1 )\n\t\t\t\t{\n\t\t\t\t\tthis.labelRadius.setTarget\n\t\t\t\t\t(\n\t\t\t\t\t\tlerp\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\tlabelOffset + .75,\n\t\t\t\t\t\t\t0,\n\t\t\t\t\t\t\tnLabelOffsets[index] + .5,\n\t\t\t\t\t\t\tcompressedRadii[index],\n\t\t\t\t\t\t\tcompressedRadii[index + 1]\n\t\t\t\t\t\t)\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.labelRadius.setTarget((compressedRadii[index] + compressedRadii[index + 1]) / 2);\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tradiusCenter =\n\t\t\t\t\tnodeRadius * (depth - 1) +\n\t\t\t\t\tnodeRadius / 2;\n\t\t\t\twidth = nodeRadius;\n\t\t\t\t\n\t\t\t\tthis.labelRadius.setTarget\n\t\t\t\t(\n\t\t\t\t\tradiusCenter + width * ((labelOffset + 1) / (nLabelOffsets[index] + 1) - .5)\n\t\t\t\t);\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( ! this.hide && ! this.keyed && nLabelOffsets[index] )\n\t\t{\n\t\t\t// check last and first labels in each track for overlap\n\t\t\t\n\t\t\tfor ( var i = 0; i < maxDisplayDepth - 1; i++ )\n\t\t\t{\n\t\t\t\tfor ( var j = 0; j <= nLabelOffsets[i]; j++ )\n\t\t\t\t{\n\t\t\t\t\tvar last = labelLastNodes[i][j];\n\t\t\t\t\tvar first = labelFirstNodes[i][j];\n\t\t\t\t\t\n\t\t\t\t\tif ( last )\n\t\t\t\t\t{\n\t\t\t\t\t\tif ( j == nLabelOffsets[i] )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\t// last is radial\n\t\t\t\t\t\t\tthis.setLabelWidth(last);\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tlast.setLabelWidth(this);\n\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tif ( first )\n\t\t\t\t\t{\n\t\t\t\t\t\tif ( j == nLabelOffsets[i] )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tthis.setLabelWidth(first);\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tfirst.setLabelWidth(this);\n\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\tif ( selectedNode.canDisplayLabelOther )\n\t\t\t{\n\t\t\t\tthis.setLabelWidth(selectedNode); // in case there is an \'other\' label\n\t\t\t}\n\t\t\t\n\t\t\tif ( this.radial )\n\t\t\t{\n\t\t\t\t// use the last \'track\' of this depth for radial\n\t\t\t\t\n\t\t\t\tlabelLastNodes[index][nLabelOffsets[index]] = this;\n\t\t\t\t\n\t\t\t\tif ( labelFirstNodes[index][nLabelOffsets[index]] == 0 )\n\t\t\t\t{\n\t\t\t\t\tlabelFirstNodes[index][nLabelOffsets[index]] = this;\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tlabelLastNodes[index][labelOffset] = this;\n\t\t\t\t\n\t\t\t\t// update offset\n\t\t\t\t\n\t\t\t\tlabelOffsets[index] += 1;\n\t\t\t\t\n\t\t\t\tif ( labelOffsets[index] > nLabelOffsets[index] )\n\t\t\t\t{\n\t\t\t\t\tlabelOffsets[index] -= nLabelOffsets[index];\n\t\t\t\t\t\n\t\t\t\t\tif ( !(nLabelOffsets[index] & 1) )\n\t\t\t\t\t{\n\t\t\t\t\t\tlabelOffsets[index]--;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\telse if ( labelOffsets[index] == nLabelOffsets[index] )\n\t\t\t\t{\n\t\t\t\t\tlabelOffsets[index] -= nLabelOffsets[index];\n\t\t\t\t\t\n\t\t\t\t\tif ( false && !(nLabelOffsets[index] & 1) )\n\t\t\t\t\t{\n\t\t\t\t\t\tlabelOffsets[index]++;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( labelFirstNodes[index][labelOffset] == 0 )\n\t\t\t\t{\n\t\t\t\t\tlabelFirstNodes[index][labelOffset] = this;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\telse if ( this.hide )\n\t\t{\n\t\t\tthis.labelWidth.end = 0;\n\t\t}\n\t}\n\t\n\tthis.setTargets = function()\n\t{\n\t\tif ( this == selectedNode )\n\t\t{\n\t\t\tthis.setTargetsSelected\n\t\t\t(\n\t\t\t\t0,\n\t\t\t\t1,\n\t\t\t\tlightnessBase,\n\t\t\t\tfalse,\n\t\t\t\tfalse\n\t\t\t);\n\t\t\treturn;\n\t\t}\n\t\t\n\t\tvar depthRelative = this.getDepth() - selectedNode.getDepth();\n\t\t\n\t\tvar parentOfSelected = selectedNode.hasParent(this);\n/*\t\t(\n//\t\t\t! this.getCollapse() &&\n\t\t\tthis.baseMagnitude <= selectedNode.baseMagnitude &&\n\t\t\tthis.baseMagnitude + this.magnitude >=\n\t\t\tselectedNode.baseMagnitude + selectedNode.magnitude\n\t\t);\n*/\t\t\n\t\tif ( parentOfSelected )\n\t\t{\n\t\t\tthis.resetLabelWidth();\n\t\t}\n\t\telse\n\t\t{\n\t\t\t//context.font = fontNormal;\n\t\t\tvar dim = context.measureText(this.name);\n\t\t\tthis.nameWidth = dim.width;\n\t\t\t//this.labelWidth.setTarget(this.labelWidth.end);\n\t\t\tthis.labelWidth.setTarget(0);\n\t\t}\n\t\t\n\t\t// set angles\n\t\t//\n\t\tif ( this.baseMagnitude <= selectedNode.baseMagnitude )\n\t\t{\n\t\t\tthis.angleStart.setTarget(0);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.angleStart.setTarget(Math.PI * 2);\n\t\t}\n\t\t//\n\t\tif\n\t\t(\n\t\t\tparentOfSelected ||\n\t\t\tthis.baseMagnitude + this.magnitude >=\n\t\t\tselectedNode.baseMagnitude + selectedNode.magnitude\n\t\t)\n\t\t{\n\t\t\tthis.angleEnd.setTarget(Math.PI * 2);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.angleEnd.setTarget(0);\n\t\t}\n\t\t\n\t\t// children\n\t\t//\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.children[i].setTargets();\n\t\t}\n\t\t\n\t\tif ( this.getDepth() <= selectedNode.getDepth() )\n\t\t{\n\t\t\t// collapse in\n\t\t\t\n\t\t\tthis.radiusInner.setTarget(0);\n\t\t\t\n\t\t\tif ( parentOfSelected )\n\t\t\t{\n\t\t\t\tthis.labelRadius.setTarget\n\t\t\t\t(\n\t\t\t\t\t(depthRelative) *\n\t\t\t\t\thistorySpacingFactor * fontSize / gRadius\n\t\t\t\t);\n\t\t\t\t//this.scale.setTarget(1 - (selectedNode.getDepth() - this.getDepth()) / 18); // TEMP\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.labelRadius.setTarget(0);\n\t\t\t\t//this.scale.setTarget(1); // TEMP\n\t\t\t}\n\t\t}\n\t\telse if ( depthRelative + 1 > maxDisplayDepth )\n\t\t{\n\t\t\t// collapse out\n\t\t\t\n\t\t\tthis.radiusInner.setTarget(1);\n\t\t\tthis.labelRadius.setTarget(1);\n\t\t\t//this.scale.setTarget(1); // TEMP\n\t\t}\n\t\telse\n\t\t{\n\t\t\t// don\'t collapse\n\t\t\t\n\t\t\tif ( compress )\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(compressedRadii[depthRelative - 1]);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(nodeRadius * (depthRelative));\n\t\t\t}\n\t\t\t\n\t\t\t//this.scale.setTarget(1); // TEMP\n\t\t\t\n\t\t\tif ( this == selectedNode )\n\t\t\t{\n\t\t\t\tthis.labelRadius.setTarget(0);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tif ( compress )\n\t\t\t\t{\n\t\t\t\t\tthis.labelRadius.setTarget\n\t\t\t\t\t(\n\t\t\t\t\t\t(compressedRadii[depthRelative - 1] + compressedRadii[depthRelative]) / 2\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.labelRadius.setTarget(nodeRadius * (depthRelative) + nodeRadius / 2);\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n//\t\tthis.r.start = this.r.end;\n//\t\tthis.g.start = this.g.end;\n//\t\tthis.b.start = this.b.end;\n\t\t\n\t\tthis.r.setTarget(255);\n\t\tthis.g.setTarget(255);\n\t\tthis.b.setTarget(255);\n\t\tthis.alphaLine.setTarget(0);\n\t\tthis.alphaArc.setTarget(0);\n\t\tthis.alphaWedge.setTarget(0);\n\t\tthis.alphaPattern.setTarget(0);\n\t\tthis.alphaOther.setTarget(0);\n\t\t\n\t\tif ( parentOfSelected && ! this.getCollapse() )\n\t\t{\n\t\t\tvar alpha =\n\t\t\t(\n\t\t\t\t1 -\n\t\t\t\t(selectedNode.getDepth() - this.getDepth()) /\n\t\t\t\t(Math.floor((compress ? compressedRadii[0] : nodeRadius) * gRadius / (historySpacingFactor * fontSize) - .5) + 1)\n\t\t\t);\n\t\t\t\n\t\t\tif ( alpha < 0 )\n\t\t\t{\n\t\t\t\talpha = 0;\n\t\t\t}\n\t\t\t\n\t\t\tthis.alphaLabel.setTarget(alpha);\n\t\t\tthis.radial = false;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.alphaLabel.setTarget(0);\n\t\t}\n\t\t\n\t\tthis.hideAlonePrev = this.hideAlone;\n\t\tthis.hidePrev = this.hide;\n\t\t\n\t\tif ( parentOfSelected )\n\t\t{\n\t\t\tthis.hideAlone = false;\n\t\t\tthis.hide = false;\n\t\t}\n\t\t\n\t\tif ( this.getParent() == selectedNode.getParent() )\n\t\t{\n\t\t\tthis.hiddenEnd = null;\n\t\t}\n\t\t\n\t\tthis.radialPrev = this.radial;\n\t}\n\t\n\tthis.setTargetsSelected = function(hueMin, hueMax, lightness, hide, nextSiblingHidden)\n\t{\n\t\tvar collapse = this.getCollapse();\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n\t\tvar canDisplayChildLabels = false;\n\t\tvar lastChild;\n\t\t\n\t\tif ( this.hasChildren() )//&& ! hide )\n\t\t{\n\t\t\tlastChild = this.children[this.children.length - 1];\n\t\t\tthis.hideAlone = true;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.hideAlone = false;\n\t\t}\n\t\t\n\t\t// set child wedges\n\t\t//\n\t\tfor ( var i = 0; i < this.children.length; i++ )\n\t\t{\n\t\t\tthis.children[i].setTargetWedge();\n\t\t\t\n\t\t\tif\n\t\t\t(\n\t\t\t\t! this.children[i].hide &&\n\t\t\t\t( collapse || depth < maxDisplayDepth ) &&\n\t\t\t\tthis.depth < maxAbsoluteDepth\n\t\t\t)\n\t\t\t{\n\t\t\t\tcanDisplayChildLabels = true;\n\t\t\t\tthis.hideAlone = false;\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( this == selectedNode || lastChild && lastChild.angleEnd.end < this.angleEnd.end - .01)\n\t\t{\n\t\t\tthis.hideAlone = false;\n\t\t}\n\t\t\n\t\tif ( this.hideAlonePrev == undefined )\n\t\t{\n\t\t\tthis.hideAlonePrev = this.hideAlone;\n\t\t}\n\t\t\n\t\tif ( this == selectedNode )\n\t\t{\n\t\t\tvar otherArc = \n\t\t\t\tthis.children.length ?\n\t\t\t\t\tangleFactor *\n\t\t\t\t\t(\n\t\t\t\t\t\tthis.baseMagnitude + this.magnitude -\n\t\t\t\t\t\tlastChild.baseMagnitude - lastChild.magnitude\n\t\t\t\t\t)\n\t\t\t\t: this.baseMagnitude + this.magnitude;\n\t\t\tthis.canDisplayLabelOther =\n\t\t\t\tthis.children.length ?\n\t\t\t\t\totherArc *\n\t\t\t\t\t(this.children[0].radiusInner.end + 1) * gRadius >=\n\t\t\t\t\tminWidth()\n\t\t\t\t: true;\n\t\t\t\n\t\t\tthis.keyUnclassified = false;\n\t\t\t\n\t\t\tif ( this.canDisplayLabelOther )\n\t\t\t{\n\t\t\t\tthis.angleOther = Math.PI * 2 - otherArc / 2;\n\t\t\t}\n\t\t\telse if ( otherArc > 0.0000000001 )\n\t\t\t{\n\t\t\t\tthis.keyUnclassified = true;\n\t\t\t\tkeys++;\n\t\t\t}\n\t\t\t\n\t\t\tthis.angleStart.setTarget(0);\n\t\t\tthis.angleEnd.setTarget(Math.PI * 2);\n\t\t\t\n\t\t\tif ( this.children.length )\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(0);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(compressedRadii[0]);\n\t\t\t}\n\t\t\t\n\t\t\tthis.hidePrev = this.hide;\n\t\t\tthis.hide = false;\n\t\t\tthis.hideAlonePrev = this.hideAlone;\n\t\t\tthis.hideAlone = false;\n\t\t\tthis.keyed = false;\n\t\t}\n\t\t\n\t\tif ( hueMax - hueMin > 1 / 12 )\n\t\t{\n\t\t\thueMax = hueMin + 1 / 12;\n\t\t}\n\t\t\n\t\t// set lightness\n\t\t//\n\t\tif ( ! ( hide || this.hideAlone ) )\n\t\t{\n\t\t\tif ( useHue() )\n\t\t\t{\n\t\t\t\tlightness = (lightnessBase + lightnessMax) / 2;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tlightness = lightnessBase + (depth - 1) * lightnessFactor;\n\t\t\t\t\n\t\t\t\tif ( lightness > lightnessMax )\n\t\t\t\t{\n\t\t\t\t\tlightness = lightnessMax;\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( hide )\n\t\t{\n\t\t\tthis.hide = true;\n\t\t}\n\t\t\n\t\tif ( this.hidePrev == undefined )\n\t\t{\n\t\t\tthis.hidePrev = this.hide;\n\t\t}\n\t\t\n\t\tvar hiddenStart = -1;\n\t\tvar hiddenHueNumer = 0;\n\t\tvar hiddenHueDenom = 0;\n\t\tvar i = 0;\n\t\t\n\t\tif ( ! this.hide )\n\t\t{\n\t\t\tthis.hiddenEnd = null;\n\t\t}\n\t\t\n\t\twhile ( true )\n\t\t{\n\t\t\tif ( ! this.hideAlone && ! hide && ( i == this.children.length || ! this.children[i].hide ) )\n\t\t\t{\n\t\t\t\t// reached a non-hidden child or the end; set targets for\n\t\t\t\t// previous group of hidden children (if any) using their\n\t\t\t\t// average hue\n\t\t\t\t\n\t\t\t\tif ( hiddenStart != -1 )\n\t\t\t\t{\n\t\t\t\t\tvar hiddenHue = hiddenHueDenom ? hiddenHueNumer / hiddenHueDenom : hueMin;\n\t\t\t\t\t\n\t\t\t\t\tfor ( var j = hiddenStart; j < i; j++ )\n\t\t\t\t\t{\n\t\t\t\t\t\tthis.children[j].setTargetsSelected\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\thiddenHue,\n\t\t\t\t\t\t\tnull,\n\t\t\t\t\t\t\tlightness,\n\t\t\t\t\t\t\tfalse,\n\t\t\t\t\t\t\tj < i - 1\n\t\t\t\t\t\t);\n\t\t\t\t\t\t\n\t\t\t\t\t\tthis.children[j].hiddenEnd = null;\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tthis.children[hiddenStart].hiddenEnd = i - 1;\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\tif ( i == this.children.length )\n\t\t\t{\n\t\t\t\tbreak;\n\t\t\t}\n\t\t\t\n\t\t\tvar child = this.children[i];\n\t\t\tvar childHueMin;\n\t\t\tvar childHueMax;\n\t\t\t\n\t\t\tif ( this.magnitude > 0 && ! this.hide && ! this.hideAlone )\n\t\t\t{\n\t\t\t\tif ( useHue() )\n\t\t\t\t{\n\t\t\t\t\tchildHueMin = child.hues[currentDataset];\n\t\t\t\t}\n\t\t\t\telse if ( this == selectedNode )\n\t\t\t\t{\n\t\t\t\t\tvar min = 0.0;\n\t\t\t\t\tvar max = 1.0;\n\t\t\t\t\t\n\t\t\t\t\tif ( this.children.length > 6 )\n\t\t\t\t\t{\n\t\t\t\t\t\tchildHueMin = lerp((1 - Math.pow(1 - i / this.children.length, 1.4)) * .95, 0, 1, min, max);\n\t\t\t\t\t\tchildHueMax = lerp((1 - Math.pow(1 - (i + .55) / this.children.length, 1.4)) * .95, 0, 1, min, max);\n\t\t\t\t\t}\n\t\t\t\t\telse\n\t\t\t\t\t{\n\t\t\t\t\t\tchildHueMin = lerp(i / this.children.length, 0, 1, min, max);\n\t\t\t\t\t\tchildHueMax = lerp((i + .55) / this.children.length, 0, 1, min, max);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tchildHueMin = lerp\n\t\t\t\t\t(\n\t\t\t\t\t\tchild.baseMagnitude,\n\t\t\t\t\t\tthis.baseMagnitude, \n\t\t\t\t\t\tthis.baseMagnitude + this.magnitude,\n\t\t\t\t\t\thueMin,\n\t\t\t\t\t\thueMax\n\t\t\t\t\t);\n\t\t\t\t\tchildHueMax = lerp\n\t\t\t\t\t(\n\t\t\t\t\t\tchild.baseMagnitude + child.magnitude * .99,\n\t\t\t\t\t\tthis.baseMagnitude,\n\t\t\t\t\t\tthis.baseMagnitude + this.magnitude,\n\t\t\t\t\t\thueMin,\n\t\t\t\t\t\thueMax\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tchildHueMin = hueMin;\n\t\t\t\tchildHueMax = hueMax;\n\t\t\t}\n\t\t\t\n\t\t\tif ( ! this.hideAlone && ! hide && ! this.hide && child.hide )\n\t\t\t{\n\t\t\t\tif ( hiddenStart == -1 )\n\t\t\t\t{\n\t\t\t\t\thiddenStart = i;\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( useHue() )\n\t\t\t\t{\n\t\t\t\t\thiddenHueNumer += childHueMin * child.magnitude;\n\t\t\t\t\thiddenHueDenom += child.magnitude;\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\thiddenHueNumer += childHueMin;\n\t\t\t\t\thiddenHueDenom++;\n\t\t\t\t}\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\thiddenStart = -1;\n\t\t\t\t\n\t\t\t\tthis.children[i].setTargetsSelected\n\t\t\t\t(\n\t\t\t\t\tchildHueMin,\n\t\t\t\t\tchildHueMax,\n\t\t\t\t\tlightness,\n\t\t\t\t\thide || this.keyed || this.hideAlone || this.hide && ! collapse,\n\t\t\t\t\tfalse\n\t\t\t\t);\n\t\t\t}\n\t\t\t\n\t\t\ti++;\n\t\t}\n\t\t\n\t \tif ( this.hue && this.magnitude )\n\t \t{\n\t\t \tthis.hue.setTarget(this.hues[currentDataset]);\n\t\t\t\n\t\t\tif ( this.attributes[magnitudeIndex][lastDataset] == 0 )\n\t\t\t{\n\t\t\t\tthis.hue.start = this.hue.end;\n\t\t\t}\n\t\t}\n\t \t\n\t\tthis.radialPrev = this.radial;\n\t\t\n\t\tif ( this == selectedNode )\n\t\t{\n\t\t\tthis.resetLabelWidth();\n\t\t\tthis.labelWidth.setTarget(this.nameWidth * labelWidthFudge);\n\t\t\tthis.alphaWedge.setTarget(0);\n\t\t\tthis.alphaLabel.setTarget(1);\n\t\t\tthis.alphaOther.setTarget(1);\n\t\t\tthis.alphaArc.setTarget(0);\n\t\t\tthis.alphaLine.setTarget(0);\n\t\t\tthis.alphaPattern.setTarget(0);\n\t\t\tthis.r.setTarget(255);\n\t\t\tthis.g.setTarget(255);\n\t\t\tthis.b.setTarget(255);\n\t\t\tthis.radial = false;\n\t\t\tthis.labelRadius.setTarget(0);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tvar rgb = hslToRgb\n\t\t\t(\n\t\t\t\thueMin,\n\t\t\t\tsaturation,\n\t\t\t\tlightness\n\t\t\t);\n\t\t\t\n\t\t\tthis.r.setTarget(rgb.r);\n\t\t\tthis.g.setTarget(rgb.g);\n\t\t\tthis.b.setTarget(rgb.b);\n\t\t\tthis.alphaOther.setTarget(0);\n\t\t\t\n\t\t\tthis.alphaWedge.setTarget(1);\n\t\t\t\n\t\t\tif ( this.hide || this.hideAlone )\n\t\t\t{\n\t\t\t\tthis.alphaPattern.setTarget(1);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.alphaPattern.setTarget(0);\n\t\t\t}\n\t\t\t\n\t\t\t// set radial\n\t\t\t//\n\t\t\tif ( ! ( hide || this.hide ) )//&& ! this.keyed )\n\t\t\t{\n\t\t\t\tif ( this.hideAlone )\n\t\t\t\t{\n\t\t\t\t\tthis.radial = true;\n\t\t\t\t}\n\t\t\t\telse if ( false && canDisplayChildLabels )\n\t\t\t\t{\n\t\t\t\t\tthis.radial = false;\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.radial = true;\n\t\t\t\t\t\n\t\t\t\t\tif ( this.hasChildren() && depth < maxDisplayDepth )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar lastChild = this.children[this.children.length - 1];\n\t\t\t\t\t\t\n\t\t\t\t\t\tif\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\tlastChild.angleEnd.end == this.angleEnd.end ||\n\t\t\t\t\t\t\t(\n\t\t\t\t\t\t\t\t(this.angleStart.end + this.angleEnd.end) / 2 -\n\t\t\t\t\t\t\t\tlastChild.angleEnd.end\n\t\t\t\t\t\t\t) * (this.radiusInner.end + 1) * gRadius * 2 <\n\t\t\t\t\t\t\tminWidth()\n\t\t\t\t\t\t)\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tthis.radial = false;\n\t\t\t\t\t\t}\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\t// set alphaLabel\n\t\t\t//\n\t\t\tif\n\t\t\t(\n\t\t\t\tcollapse ||\n\t\t\t\thide ||\n\t\t\t\tthis.hide ||\n\t\t\t\tthis.keyed ||\n\t\t\t\tdepth > maxDisplayDepth ||\n\t\t\t\t! this.canDisplayDepth()\n\t\t\t)\n\t\t\t{\n\t\t\t\tthis.alphaLabel.setTarget(0);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tif\n\t\t\t\t(\n\t\t\t\t\t(this.radial || nLabelOffsets[depth - 2])\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tthis.alphaLabel.setTarget(1);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.alphaLabel.setTarget(0);\n\t\t\t\t\t\n\t\t\t\t\tif ( this.radialPrev )\n\t\t\t\t\t{\n\t\t\t\t\t\tthis.alphaLabel.start = 0;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t}\n\t\t\t\n\t\t\t// set alphaArc\n\t\t\t//\n\t\t\tif\n\t\t\t(\n\t\t\t\tcollapse ||\n\t\t\t\thide ||\n\t\t\t\tdepth > maxDisplayDepth ||\n\t\t\t\t! this.canDisplayDepth()\n\t\t\t)\n\t\t\t{\n\t\t\t\tthis.alphaArc.setTarget(0);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.alphaArc.setTarget(1);\n\t\t\t}\n\t\t\t\n\t\t\t// set alphaLine\n\t\t\t//\n\t\t\tif\n\t\t\t(\n\t\t\t\thide ||\n\t\t\t\tthis.hide && nextSiblingHidden ||\n\t\t\t\tdepth > maxDisplayDepth ||\n\t\t\t\t! this.canDisplayDepth()\n\t\t\t)\n\t\t\t{\n\t\t\t\tthis.alphaLine.setTarget(0);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.alphaLine.setTarget(1);\n\t\t\t}\n\t\t\t\n\t\t\t//if (  ! this.radial )\n\t\t\t{\n\t\t\t\tthis.resetLabelWidth();\n\t\t\t}\n\t\t\t\n\t\t\t// set labelRadius target\n\t\t\t//\n\t\t\tif ( collapse )\n\t\t\t{\n\t\t\t\tthis.labelRadius.setTarget(this.radiusInner.end);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tif ( depth > maxDisplayDepth || ! this.canDisplayDepth() )\n\t\t\t\t{\n\t\t\t\t\tthis.labelRadius.setTarget(1);\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tthis.setTargetLabelRadius();\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t}\n\t\n\tthis.setTargetWedge = function()\n\t{\n\t\tvar depth = this.getDepth() - selectedNode.getDepth() + 1;\n\t\t\n\t\t// set angles\n\t\t//\n\t\tvar baseMagnitudeRelative = this.baseMagnitude - selectedNode.baseMagnitude;\n\t\t//\n\t\tthis.angleStart.setTarget(baseMagnitudeRelative * angleFactor);\n\t\tthis.angleEnd.setTarget((baseMagnitudeRelative + this.magnitude) * angleFactor);\n\t\t\n\t\t// set radiusInner\n\t\t//\n\t\tif ( depth > maxDisplayDepth || ! this.canDisplayDepth() )\n\t\t{\n\t\t\tthis.radiusInner.setTarget(1);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tif ( compress )\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(compressedRadii[depth - 2]);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.radiusInner.setTarget(nodeRadius * (depth - 1));\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( this.hide != undefined )\n\t\t{\n\t\t\tthis.hidePrev = this.hide;\n\t\t}\n\t\t\n\t\tif ( this.hideAlone != undefined )\n\t\t{\n\t\t\tthis.hideAlonePrev = this.hideAlone;\n\t\t}\n\t\t\n\t\t// set hide\n\t\t//\n\t\tif\n\t\t(\n\t\t\t(this.angleEnd.end - this.angleStart.end) *\n\t\t\t(this.radiusInner.end * gRadius + gRadius) <\n\t\t\tminWidth()\n\t\t)\n\t\t{\n\t\t\tif ( depth == 2 && ! this.getCollapse() && this.depth <= maxAbsoluteDepth )\n\t\t\t{\n\t\t\t\tthis.keyed = true;\n\t\t\t\tkeys++;\n\t\t\t\tthis.hide = false;\n\t\t\t\t\n\t\t\t\tvar percentage = this.getPercentage();\n\t\t\t\tthis.keyLabel = this.name + \'   \' + percentage + \'%\';\n\t\t\t\tvar dim = context.measureText(this.keyLabel);\n\t\t\t\tthis.keyNameWidth = dim.width;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tthis.keyed = false;\n\t\t\t\tthis.hide = depth > 2;\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tthis.keyed = false;\n\t\t\tthis.hide = false;\n\t\t}\n\t}\n\t\n\tthis.shortenLabel = function()\n\t{\n\t\tvar label = this.name;\n\t\t\n\t\tvar labelWidth = this.nameWidth;\n\t\tvar maxWidth = this.labelWidth.current();\n\t\tvar minEndLength = 0;\n\t\t\n\t\tif ( labelWidth > maxWidth && label.length > minEndLength * 2 )\n\t\t{\n\t\t\tvar endLength =\n\t\t\t\tMath.floor((label.length - 1) * maxWidth / labelWidth / 2);\n\t\t\t\n\t\t\tif ( endLength < minEndLength )\n\t\t\t{\n\t\t\t\tendLength = minEndLength;\n\t\t\t}\n\t\t\t\n\t\t\treturn (\n\t\t\t\tlabel.substring(0, endLength) +\n\t\t\t\t\'...\' +\n\t\t\t\tlabel.substring(label.length - endLength));\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn label;\n\t\t}\n\t}\n\t\n/*\tthis.shouldAddSearchResultsString = function()\n\t{\n\t\tif ( this.isSearchResult )\n\t\t{\n\t\t\treturn this.searchResults > 1;\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn this.searchResults > 0;\n\t\t}\n\t}\n*/\t\n\tthis.sort = function()\n\t{\n\t\tthis.children.sort(function(a, b){return b.getMagnitude() - a.getMagnitude()});\n\t\t\n\t\tfor (var i = 0; i < this.children.length; i++)\n\t\t{\n\t\t\tthis.children[i].sort();\n\t\t}\n\t}\n}\nvar options;\nfunction addOptionElement(position, innerHTML, title)\n{\n\tvar div = document.createElement("div");\n//\tdiv.style.position = \'absolute\';\n//\tdiv.style.top = position + \'px\';\n\tdiv.innerHTML = innerHTML;\n//\tdiv.style.display = \'block\';\n\tdiv.style.padding = \'2px\';\n\t\n\tif ( title )\n\t{\n\t\tdiv.title = title;\n\t}\n\t\n\toptions.appendChild(div);\n\tvar height = 0;//div.clientHeight;\n\treturn position + height;\n}\nfunction addOptionElements(hueName, hueDefault)\n{\n\toptions = document.createElement(\'div\');\n\toptions.style.position = \'absolute\';\n\toptions.style.top = \'0px\';\n\toptions.addEventListener(\'mousedown\', function(e) {mouseClick(e)}, false);\n//\toptions.onmouseup = function(e) {mouseUp(e)}\n\tdocument.body.appendChild(options);\n\t\n\tdocument.body.style.font = \'11px sans-serif\';\n\tvar position = 5;\n\t\n\tdetails = document.createElement(\'div\');\n\tdetails.style.position = \'absolute\';\n\tdetails.style.top = \'1%\';\n\tdetails.style.right = \'2%\';\n\tdetails.style.textAlign = \'right\';\n\tdocument.body.insertBefore(details, canvas);\n//\t\t<div id="details" style="position:absolute;top:1%;right:2%;text-align:right;">\n\tdetails.innerHTML = \'\\\n<span id="detailsName" style="font-weight:bold"></span>&nbsp;\\\n<input type="button" id="detailsExpand" onclick="expand(focusNode);"\\\nvalue="&harr;" title="Expand this wedge to become the new focus of the chart"/><br/>\\\n<div id="detailsInfo" style="float:right"></div>\';\n\tkeyControl = document.createElement(\'input\');\n\tkeyControl.type = \'button\';\n\tkeyControl.value = showKeys ? \'x\' : \'...\';\n\tkeyControl.style.position = \'\';\n\tkeyControl.style.position = \'fixed\';\n\tkeyControl.style.visibility = \'hidden\';\n\t\n\tdocument.body.insertBefore(keyControl, canvas);\n\t\n\tvar logoElement = document.getElementById(\'logo\');\n\t\n\tif ( logoElement )\n\t{\n\t\tlogoImage = logoElement.src;\n\t}\n\telse\n\t{\n\t\tlogoImage = \'http://marbl.github.io/Krona/img/logo-med.png\';\n\t}\n\t\n//\tdocument.getElementById(\'options\').style.fontSize = \'9pt\';\n\tposition = addOptionElement\n\t(\n\t\tposition,\n\'<a style="margin:2px" target="_blank" href="https://github.com/marbl/Krona/wiki"><img style="vertical-align:middle;width:108px;height:30px;" src="\' + logoImage + \'"/></a><input type="button" id="back" value="&larr;" title="Go back (Shortcut: &larr;)"/>\\\n<input type="button" id="forward" value="&rarr;" title="Go forward (Shortcut: &rarr;)"/> \\\n&nbsp;Search: <input type="text" id="search"/>\\\n<input id="searchClear" type="button" value="x" onclick="clearSearch()"/> \\\n<span id="searchResults"></span>\'\n\t);\n\t\n\tif ( datasets > 1 )\n\t{\n\t\tvar size = datasets < datasetSelectSize ? datasets : datasetSelectSize;\n\t\t\n\t\tvar select =\n\t\t\t\'<table style="border-collapse:collapse;padding:0px"><tr><td style="padding:0px">\' +\n\t\t\t\'<select id="datasets" style="min-width:100px" size="\' + size + \'" onchange="onDatasetChange()">\';\n\t\t\n\t\tfor ( var i = 0; i < datasetNames.length; i++ )\n\t\t{\n\t\t\tselect += \'<option>\' + datasetNames[i] + \'</option>\';\n\t\t}\n\t\t\n\t\tselect +=\n\t\t\t\'</select></td><td style="vertical-align:top;padding:1px;">\' +\n\t\t\t\'<input style="display:block" title="Previous dataset (Shortcut: &uarr;)" id="prevDataset" type="button" value="&uarr;" onclick="prevDataset()" disabled="true"/>\' +\n\t\t\t\'<input title="Next dataset (Shortcut: &darr;)" id="nextDataset" type="button" value="&darr;" onclick="nextDataset()"/><br/></td>\' +\n\t\t\t\'<td style="padding-top:1px;vertical-align:top"><input title="Switch to the last dataset that was viewed (Shortcut: TAB)" id="lastDataset" type="button" style="font:11px Times new roman" value="last" onclick="selectLastDataset()"/></td></tr></table>\';\n\t\t\n\t\tposition = addOptionElement(position + 5, select);\n\t\t\n\t\tdatasetDropDown = document.getElementById(\'datasets\');\n\t\tdatasetButtonLast = document.getElementById(\'lastDataset\');\n\t\tdatasetButtonPrev = document.getElementById(\'prevDataset\');\n\t\tdatasetButtonNext = document.getElementById(\'nextDataset\');\n\t\t\n\t\tposition += datasetDropDown.clientHeight;\n\t}\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition + 5,\n\'<input type="button" id="maxAbsoluteDepthDecrease" value="-"/>\\\n<span id="maxAbsoluteDepth"></span>\\\n&nbsp;<input type="button" id="maxAbsoluteDepthIncrease" value="+"/> Max depth\',\n\'Maximum depth to display, counted from the top level \\\nand including collapsed wedges.\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition,\n\'<input type="button" id="fontSizeDecrease" value="-"/>\\\n<span id="fontSize"></span>\\\n&nbsp;<input type="button" id="fontSizeIncrease" value="+"/> Font size\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition,\n\'<input type="button" id="radiusDecrease" value="-"/>\\\n<input type="button" id="radiusIncrease" value="+"/> Chart size\'\n\t);\n\t\n\tif ( hueName )\n\t{\n\t\thueDisplayName = attributes[attributeIndex(hueName)].displayName;\n\t\t\n\t\tposition = addOptionElement\n\t\t(\n\t\t\tposition + 5,\n\t\t\t\'<input type="checkbox" id="useHue" style="float:left" \' +\n\t\t\t\'/><div>Color by<br/>\' + hueDisplayName +\n\t\t\t\'</div>\'\n\t\t);\n\t\t\n\t\tuseHueCheckBox = document.getElementById(\'useHue\');\n\t\tuseHueCheckBox.checked = hueDefault;\n\t\tuseHueCheckBox.onclick = handleResize;\n\t\tuseHueCheckBox.onmousedown = suppressEvent;\n\t}\n\t/*\n\tposition = addOptionElement\n\t(\n\t\tposition + 5,\n\t\t\'&nbsp;<input type="checkbox" id="shorten" checked="checked" />Shorten labels</div>\',\n\t\t\'Prevent labels from overlapping by shortening them\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition,\n\t\t\'&nbsp;<input type="checkbox" id="compress" checked="checked" />Compress\',\n\t\t\'Compress wedges if needed to show the entire depth\'\n\t);\n\t*/\n\tposition = addOptionElement\n\t(\n\t\tposition,\n\t\t\'<input type="checkbox" id="collapse" checked="checked" />Collapse\',\n\t\t\'Collapse wedges that are redundant (entirely composed of another wedge)\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition + 5,\n\t\t\'<input type="button" id="snapshot" value="Snapshot"/>\',\n\'Render the current view as SVG (Scalable Vector Graphics), a publication-\\\nquality format that can be printed and saved (see Help for browser compatibility)\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition + 5,\n\'<input type="button" id="linkButton" value="Link"/>\\\n<input type="text" size="30" id="linkText"/>\',\n\'Show a link to this view that can be copied for bookmarking or sharing\'\n\t);\n\t\n\tposition = addOptionElement\n\t(\n\t\tposition + 5,\n\'<input type="button" id="help" value="?"\\\nonclick="window.open(\\\'https://github.com/marbl/Krona/wiki/Browsing%20Krona%20charts\\\', \\\'help\\\')"/>\',\n\'Help\'\n\t);\n}\nfunction arrow(angleStart, angleEnd, radiusInner)\n{\n\tif ( context.globalAlpha == 0 )\n\t{\n\t\treturn;\n\t}\n\t\n\tvar angleCenter = (angleStart + angleEnd) / 2;\n\tvar radiusArrowInner = radiusInner - gRadius / 10;//nodeRadius * gRadius;\n\tvar radiusArrowOuter = gRadius * 1.1;//(1 + nodeRadius);\n\tvar radiusArrowCenter = (radiusArrowInner + radiusArrowOuter) / 2;\n\tvar pointLength = (radiusArrowOuter - radiusArrowInner) / 5;\n\t\n\tcontext.fillStyle = highlightFill;\n\tcontext.lineWidth = highlightLineWidth;\n\t\n\t// First, mask out the first half of the arrow.  This will prevent the tips\n\t// from superimposing if the arrow goes most of the way around the circle.\n\t// Masking is done by setting the clipping region to the inverse of the\n\t// half-arrow, which is defined by cutting the half-arrow out of a large\n\t// rectangle\n\t//\n\tcontext.beginPath();\n\tcontext.arc(0, 0, radiusInner, angleCenter, angleEnd, false);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowInner * Math.cos(angleEnd),\n\t\tradiusArrowInner * Math.sin(angleEnd)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowCenter * Math.cos(angleEnd) - pointLength * Math.sin(angleEnd),\n\t\tradiusArrowCenter * Math.sin(angleEnd) + pointLength * Math.cos(angleEnd)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowOuter * Math.cos(angleEnd),\n\t\tradiusArrowOuter * Math.sin(angleEnd)\n\t);\n\tcontext.arc(0, 0, gRadius, angleEnd, angleCenter, true);\n\tcontext.closePath();\n\tcontext.moveTo(-imageWidth, -imageHeight);\n\tcontext.lineTo(imageWidth, -imageHeight);\n\tcontext.lineTo(imageWidth, imageHeight);\n\tcontext.lineTo(-imageWidth, imageHeight);\n\tcontext.closePath();\n\tcontext.save();\n\tcontext.clip();\n\t\n\t// Next, draw the other half-arrow with the first half masked out\n\t//\n\tcontext.beginPath();\n\tcontext.arc(0, 0, radiusInner, angleCenter, angleStart, true);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowInner * Math.cos(angleStart),\n\t\tradiusArrowInner * Math.sin(angleStart)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowCenter * Math.cos(angleStart) + pointLength * Math.sin(angleStart),\n\t\tradiusArrowCenter * Math.sin(angleStart) - pointLength * Math.cos(angleStart)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowOuter * Math.cos(angleStart),\n\t\tradiusArrowOuter * Math.sin(angleStart)\n\t);\n\tcontext.arc(0, 0, gRadius, angleStart, angleCenter, false);\n\tcontext.fill();\n\tcontext.stroke();\n\t\n\t// Finally, remove the clipping region and draw the first half-arrow.  This\n\t// half is extended slightly to fill the seam.\n\t//\n\tcontext.restore();\n\tcontext.beginPath();\n\tcontext.arc(0, 0, radiusInner, angleCenter - 2 / (2 * Math.PI * radiusInner), angleEnd, false);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowInner * Math.cos(angleEnd),\n\t\tradiusArrowInner * Math.sin(angleEnd)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowCenter * Math.cos(angleEnd) - pointLength * Math.sin(angleEnd),\n\t\tradiusArrowCenter * Math.sin(angleEnd) + pointLength * Math.cos(angleEnd)\n\t);\n\tcontext.lineTo\n\t(\n\t\tradiusArrowOuter * Math.cos(angleEnd),\n\t\tradiusArrowOuter * Math.sin(angleEnd)\n\t);\n\tcontext.arc(0, 0, gRadius, angleEnd, angleCenter - 2 / (2 * Math.PI * gRadius), true);\n\tcontext.fill();\n\tcontext.stroke();\n}\nfunction attributeIndex(aname)\n{\n\tfor ( var i = 0 ; i < attributes.length; i++ )\n\t{\n\t\tif ( aname == attributes[i].name )\n\t\t{\n\t\t\treturn i;\n\t\t}\n\t}\n\t\n\treturn null;\n}\nfunction checkHighlight()\n{\n\tvar lastHighlightedNode = highlightedNode;\n\tvar lastHighlightingHidden = highlightingHidden;\n\t\n\thighlightedNode = selectedNode;\n\tresetKeyOffset();\n\t\n\tif ( progress == 1 )\n\t{\n\t\tselectedNode.checkHighlight();\n\t\tif ( selectedNode.getParent() )\n\t\t{\n\t\t\tselectedNode.getParent().checkHighlightCenter();\n\t\t}\n\t\t\n\t\tfocusNode.checkHighlightMap();\n\t}\n\t\n\tif ( highlightedNode != selectedNode )\n\t{\n\t\tif ( highlightedNode == focusNode )\n\t\t{\n//\t\t\tcanvas.style.display=\'none\';\n//\t\t\twindow.resizeBy(1,0);\n//\t\t\tcanvas.style.cursor=\'ew-resize\';\n//\t\t\twindow.resizeBy(-1,0);\n//\t\t\tcanvas.style.display=\'inline\';\n\t\t}\n\t\telse\n\t\t{\n//\t\t\tcanvas.style.cursor=\'pointer\';\n\t\t}\n\t}\n\telse\n\t{\n//\t\tcanvas.style.cursor=\'auto\';\n\t}\n\t\n\tif\n\t(\n\t\t(\n\t\t\ttrue ||\n\t\t\thighlightedNode != lastHighlightedNode ||\n\t\t\thighlightingHidden != highlightingHiddenLast\n\t\t) &&\n\t\tprogress == 1\n\t)\n\t{\n\t\tdraw(); // TODO: handle in update()\n\t}\n}\nfunction checkSelectedCollapse()\n{\n\tvar newNode = selectedNode;\n\t\n\twhile ( newNode.getCollapse() )\n\t{\n\t\tnewNode = newNode.children[0];\n\t}\n\t\n\tif ( newNode.children.length == 0 && newNode.getParent() )\n\t{\n\t\tnewNode = newNode.getParent();\n\t}\n\t\n\tif ( newNode != selectedNode )\n\t{\n\t\tselectNode(newNode);\n\t}\n}\nfunction clearSearch()\n{\n\tif ( search.value != \'\' )\n\t{\n\t\tsearch.value = \'\';\n\t\tonSearchChange();\n\t}\n}\nfunction createSVG()\n{\n\tsvgNS = "http://www.w3.org/2000/svg";\n\tvar SVG = {};\n\tSVG.xlinkns = "http://www.w3.org/1999/xlink";\n\t\n\tvar newSVG = document.createElementNS(svgNS, "svg:svg");\n\t\n\tnewSVG.setAttribute("id", "canvas");\n\t// How big is the canvas in pixels\n\tnewSVG.setAttribute("width", \'100%\');\n\tnewSVG.setAttribute("height", \'100%\');\n\t// Set the coordinates used by drawings in the canvas\n//\tnewSVG.setAttribute("viewBox", "0 0 " + imageWidth + " " + imageHeight);\n\t// Define the XLink namespace that SVG uses\n\tnewSVG.setAttributeNS\n\t(\n\t\t"http://www.w3.org/2000/xmlns/",\n\t\t"xmlns:xlink",\n\t\tSVG.xlinkns\n\t);\n\t\n\treturn newSVG;\n}\nfunction degrees(radians)\n{\n\treturn radians * 180 / Math.PI;\n}\nfunction draw()\n{\n\ttweenFrames++;\n\t//resize();\n//\tcontext.fillRect(0, 0, imageWidth, imageHeight);\n\tcontext.clearRect(0, 0, imageWidth, imageHeight);\n\t\n\tcontext.font = fontNormal;\n\tcontext.textBaseline = \'middle\';\n\t\n\t//context.strokeStyle = \'rgba(0, 0, 0, 0.3)\';\n\tcontext.translate(centerX, centerY);\n\t\n\tresetKeyOffset();\n\t\n\thead.draw(false, false); // draw pie slices\n\thead.draw(true, false); // draw labels\n\t\n\tvar pathRoot = selectedNode;\n\t\n\tif ( focusNode != 0 && focusNode != selectedNode )\n\t{\n\t\tcontext.globalAlpha = 1;\n\t\tfocusNode.drawHighlight(true);\n\t\tpathRoot = focusNode;\n\t}\n\t\n\tif\n\t(\n\t\thighlightedNode &&\n\t\thighlightedNode.getDepth() >= selectedNode.getDepth() &&\n\t\thighlightedNode != focusNode\n\t)\n\t{\n\t\tif\n\t\t(\n\t\t\tprogress == 1 &&\n\t\t\thighlightedNode != selectedNode &&\n\t\t\t(\n\t\t\t\thighlightedNode != focusNode ||\n\t\t\t\tfocusNode.children.length > 0\n\t\t\t)\n\t\t)\n\t\t{\n\t\t\tcontext.globalAlpha = 1;\n\t\t\thighlightedNode.drawHighlight(true);\n\t\t}\n\t\t\n\t\t//pathRoot = highlightedNode;\n\t}\n\telse if\n\t(\n\t\tprogress == 1 &&\n\t\thighlightedNode.getDepth() < selectedNode.getDepth()\n\t)\n\t{\n\t\tcontext.globalAlpha = 1;\n\t\thighlightedNode.drawHighlightCenter();\n\t}\n\t\n\tif ( quickLook && false) // TEMP\n\t{\n\t\tcontext.globalAlpha = 1 - progress / 2;\n\t\tselectedNode.drawHighlight(true);\n\t}\n\telse if ( progress < 1 )//&& zoomOut() )\n\t{\n\t\tif ( !zoomOut)//() )\n\t\t{\n\t\t\tcontext.globalAlpha = selectedNode.alphaLine.current();\n\t\t\tselectedNode.drawHighlight(true);\n\t\t}\n\t\telse if ( selectedNodeLast )\n\t\t{\n\t\t\tcontext.globalAlpha = 1 - 4 * Math.pow(progress - .5, 2);\n\t\t\tselectedNodeLast.drawHighlight(false);\n\t\t}\n\t}\n\t\n\tdrawDatasetName();\n\t\n\t//drawHistory();\n\t\n\tcontext.translate(-centerX, -centerY);\n\tcontext.globalAlpha = 1;\n\t\n\tmapRadius =\n\t\t(imageHeight / 2 - details.clientHeight - details.offsetTop) /\n\t\t(pathRoot.getDepth() - 1) * 3 / 4 / 2;\n\t\n\tif ( mapRadius > maxMapRadius )\n\t{\n\t\tmapRadius = maxMapRadius;\n\t}\n\t\n\tmapBuffer = mapRadius / 2;\n\t\n\t//context.font = fontNormal;\n\tpathRoot.drawMap(pathRoot);\n\t\n\tif ( hueDisplayName && useHue() )\n\t{\n\t\tdrawLegend();\n\t}\n}\nfunction drawBubble(angle, radius, width, radial, flip)\n{\n\tvar height = fontSize * 2;\n\tvar x;\n\tvar y;\n\t\n\twidth = width + fontSize;\n\t\n\tif ( radial )\n\t{\n\t\ty = -fontSize;\n\t\t\n\t\tif ( flip )\n\t\t{\n\t\t\tx = radius - width + fontSize / 2;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tx = radius - fontSize / 2;\n\t\t}\n\t}\n\telse\n\t{\n\t\tx = -width / 2;\n\t\ty = -radius - fontSize;\n\t}\n\t\n\tif ( snapshotMode )\n\t{\n\t\tdrawBubbleSVG(x + centerX, y + centerY, width, height, fontSize, angle);\n\t}\n\telse\n\t{\n\t\tdrawBubbleCanvas(x, y, width, height, fontSize, angle);\n\t}\n}\nfunction drawBubbleCanvas(x, y, width, height, radius, rotation)\n{\n\tcontext.strokeStyle = \'black\';\n\tcontext.lineWidth = highlightLineWidth;\n\tcontext.fillStyle = \'rgba(255, 255, 255, .75)\';\n\tcontext.rotate(rotation);\n\troundedRectangle(x, y, width, fontSize * 2, fontSize);\n\tcontext.fill();\n\tcontext.stroke();\n\tcontext.rotate(-rotation);\n}\nfunction drawBubbleSVG(x, y, width, height, radius, rotation)\n{\n\tsvg +=\n\t\t\'<rect x="\' + x + \'" y="\' + y +\n\t\t\'" width="\' + width +\n\t\t\'" height="\' + height +\n\t\t\'" rx="\' + radius +\n\t\t\'" ry="\' + radius +\n\t\t\'" fill="rgba(255, 255, 255, .75)\' +\n\t\t\'" class="highlight" \' +\n\t\t\'transform="rotate(\' +\n\t\tdegrees(rotation) + \',\' + centerX + \',\' + centerY +\n\t\t\')"/>\';\n}\nfunction drawDatasetName()\n{\n\tvar alpha = datasetAlpha.current();\n\t\n\tif ( alpha > 0 )\n\t{\n\t\tvar radius = gRadius * compressedRadii[0] / -2;\n\t\t\n\t\tif ( alpha > 1 )\n\t\t{\n\t\t\talpha = 1;\n\t\t}\n\t\t\n\t\tcontext.globalAlpha = alpha;\n\t\t\n\t\tdrawBubble(0, -radius, datasetWidths[currentDataset], false, false);\n\t\tdrawText(datasetNames[currentDataset], 0, radius, 0, \'center\', true);\n\t}\n}\nfunction drawHistory()\n{\n\tvar alpha = 1;\n\tcontext.textAlign = \'center\';\n\t\n\tfor ( var i = 0; i < nodeHistoryPosition && alpha > 0; i++ )\n\t{\n\t\t\n\t\tcontext.globalAlpha = alpha - historyAlphaDelta * tweenFactor;\n\t\tcontext.fillText\n\t\t(\n\t\t\tnodeHistory[nodeHistoryPosition - i - 1].name,\n\t\t\t0,\n\t\t\t(i + tweenFactor) * historySpacingFactor * fontSize - 1\n\t\t);\n\t\t\n\t\tif ( alpha > 0 )\n\t\t{\n\t\t\talpha -= historyAlphaDelta;\n\t\t}\n\t}\n\t\n\tcontext.globalAlpha = 1;\n}\nfunction drawLegend()\n{\n\tvar left = imageWidth * .01;\n\tvar width = imageHeight * .0265;\n\tvar height = imageHeight * .15;\n\tvar top = imageHeight - fontSize * 3.5 - height;\n\tvar textLeft = left + width + fontSize / 2;\n\t\n\tcontext.fillStyle = \'black\';\n\tcontext.textAlign = \'start\';\n\tcontext.font = fontNormal;\n//\tcontext.fillText(valueStartText, textLeft, top + height);\n//\tcontext.fillText(valueEndText, textLeft, top);\n\tcontext.fillText(hueDisplayName, left, imageHeight - fontSize * 1.5);\n\t\n\tvar gradient = context.createLinearGradient(0, top + height, 0, top);\n\t\n\tfor ( var i = 0; i < hueStopPositions.length; i++ )\n\t{\n\t\tgradient.addColorStop(hueStopPositions[i], hueStopHsl[i]);\n\t\t\n\t\tvar textY = top + (1 - hueStopPositions[i]) * height;\n\t\t\n\t\tif\n\t\t(\n\t\t\ti == 0 ||\n\t\t\ti == hueStopPositions.length - 1 ||\n\t\t\ttextY > top + fontSize && textY < top + height - fontSize\n\t\t)\n\t\t{\n\t\t\tcontext.fillText(hueStopText[i], textLeft, textY);\n\t\t}\n\t}\n\t\n\tcontext.fillStyle = gradient;\n\tcontext.fillRect(left, top, width, height);\n\tcontext.lineWidth = thinLineWidth;\n\tcontext.strokeRect(left, top, width, height);\n}\nfunction drawLegendSVG()\n{\n\tvar left = imageWidth * .01;\n\tvar width = imageHeight * .0265;\n\tvar height = imageHeight * .15;\n\tvar top = imageHeight - fontSize * 3.5 - height;\n\tvar textLeft = left + width + fontSize / 2;\n\tvar text = \'\';\n\t\n\ttext += svgText(hueDisplayName, left, imageHeight - fontSize * 1.5);\n\t\n\tvar svgtest = \'<linearGradient id="gradient" x1="0%" y1="100%" x2="0%" y2="0%">\';\n\t\n\tfor ( var i = 0; i < hueStopPositions.length; i++ )\n\t{\n\t\tsvgtest +=\n\t\t\t\'<stop offset="\' + round(hueStopPositions[i] * 100) +\n\t\t\t\'%" style="stop-color:\' + hueStopHsl[i] + \'"/>\';\n\t\t\n\t\tvar textY = top + (1 - hueStopPositions[i]) * height;\n\t\t\n\t\tif\n\t\t(\n\t\t\ti == 0 ||\n\t\t\ti == hueStopPositions.length - 1 ||\n\t\t\ttextY > top + fontSize && textY < top + height - fontSize\n\t\t)\n\t\t{\n\t\t\ttext += svgText(hueStopText[i], textLeft, textY);\n\t\t}\n\t}\n\t\n\tsvgtest += \'</linearGradient>\';\n\t//alert(svgtest);\n\tsvg += svgtest;\n\tsvg +=\n\t\t\'<rect style="fill:url(#gradient)" x="\' + left + \'" y="\' + top +\n\t\t\'" width="\' + width + \'" height="\' + height + \'"/>\';\n\t\n\tsvg += text;\n}\nfunction drawSearchHighlights(label, bubbleX, bubbleY, rotation, center)\n{\n\tvar index = -1;\n\tvar labelLength = label.length;\n\t\n\tbubbleX -= fontSize / 4;\n\t\n\tdo\n\t{\n\t\tindex = label.toLowerCase().indexOf(search.value.toLowerCase(), index + 1);\n\t\t\n\t\tif ( index != -1 && index < labelLength )\n\t\t{\n\t\t\tvar dim = context.measureText(label.substr(0, index));\n\t\t\tvar x = bubbleX + dim.width;\n\t\t\t\n\t\t\tdim = context.measureText(label.substr(index, search.value.length));\n\t\t\t\n\t\t\tvar y = bubbleY - fontSize * 3 / 4;\n\t\t\tvar width = dim.width + fontSize / 2;\n\t\t\tvar height = fontSize * 3 / 2;\n\t\t\tvar radius = fontSize / 2;\n\t\t\t\n\t\t\tif ( snapshotMode )\n\t\t\t{\n\t\t\t\tif ( center )\n\t\t\t\t{\n\t\t\t\t\tx += centerX;\n\t\t\t\t\ty += centerY;\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tsvg +=\n\t\t\t\t\t\'<rect x="\' + x + \'" y="\' + y +\n\t\t\t\t\t\'" width="\' + width +\n\t\t\t\t\t\'" height="\' + height +\n\t\t\t\t\t\'" rx="\' + radius +\n\t\t\t\t\t\'" ry="\' + radius +\n\t\t\t\t\t\'" class="searchHighlight\' +\n\t\t\t\t\t\'" transform="rotate(\' +\n\t\t\t\t\tdegrees(rotation) + \',\' + centerX + \',\' + centerY +\n\t\t\t\t\t\')"/>\';\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tcontext.fillStyle = \'rgb(255, 255, 100)\';\n\t\t\t\tcontext.rotate(rotation);\n\t\t\t\troundedRectangle(x, y, width, height, radius);\n\t\t\t\tcontext.fill();\n\t\t\t\tcontext.rotate(-rotation);\n\t\t\t}\n\t\t}\n\t}\n\twhile ( index != -1 && index < labelLength );\n}\nfunction drawText(text, x, y, angle, anchor, bold, color)\n{\n\tif ( color == undefined )\n\t{\n\t\tcolor = \'black\';\n\t}\n\t\n\tif ( snapshotMode )\n\t{\n\t\tsvg +=\n\t\t\t\'<text x="\' + (centerX + x) + \'" y="\' + (centerY + y) +\n\t\t\t\'" text-anchor="\' + anchor + \'" style="font-color:\' + color + \';font-weight:\' + (bold ? \'bold\' : \'normal\') +\n\t\t\t\'" transform="rotate(\' + degrees(angle) + \',\' + centerX + \',\' + centerY + \')">\' +\n\t\t\ttext + \'</text>\';\n\t}\n\telse\n\t{\n\t\tcontext.fillStyle = color;\n\t\tcontext.textAlign = anchor;\n\t\tcontext.font = bold ? fontBold : fontNormal;\n\t\tcontext.rotate(angle);\n\t\tcontext.fillText(text, x, y);\n\t\tcontext.rotate(-angle);\n\t}\n}\nfunction drawTextPolar\n(\n\ttext,\n\tinnerText,\n\tangle,\n\tradius,\n\tradial,\n\tbubble,\n\tbold, \n\tsearchResult,\n\tsearchResults\n)\n{\n\tvar anchor;\n\tvar textX;\n\tvar textY;\n\tvar spacer;\n\tvar totalText = text;\n\tvar flip;\n\t\n\tif ( snapshotMode )\n\t{\n\t\tspacer = \'&#160;&#160;&#160;\';\n\t}\n\telse\n\t{\n\t\tspacer = \'   \';\n\t}\n\t\n\tif ( radial )\n\t{\n\t\tflip = angle < 3 * Math.PI / 2;\n\t\t\n\t\tif ( flip )\n\t\t{\n\t\t\tangle -= Math.PI;\n\t\t\tradius = -radius;\n\t\t\tanchor = \'end\';\n\t\t\t\n\t\t\tif ( innerText )\n\t\t\t{\n\t\t\t\ttotalText = text + spacer + innerText;\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tanchor = \'start\';\n\t\t\t\n\t\t\tif ( innerText )\n\t\t\t{\n\t\t\t\ttotalText = innerText + spacer + text;\n\t\t\t}\n\t\t}\n\t\t\n\t\ttextX = radius;\n\t\ttextY = 0;\n\t}\n\telse\n\t{\n\t\tflip = angle < Math.PI || angle > 2 * Math.PI;\n\t\tvar label;\n\t\t\n\t\tanchor = snapshotMode ? \'middle\' : \'center\';\n\t\t\n\t\tif ( flip )\n\t\t{\n\t\t\tangle -= Math.PI;\n\t\t\tradius = -radius;\n\t\t}\n\t\t\n\t\tangle += Math.PI / 2;\n\t\ttextX = 0;\n\t\ttextY = -radius;\n\t}\n\t\n\tif ( bubble )\n\t{\n\t\tvar textActual = totalText;\n\t\t\n\t\tif ( innerText && snapshotMode )\n\t\t{\n\t\t\tif ( flip )\n\t\t\t{\n\t\t\t\ttextActual = text + \'   \' + innerText;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\ttextActual = innerText + \'   \' + text;\n\t\t\t}\n\t\t}\n\t\t\n\t\tif ( searchResults )\n\t\t{\n\t\t\ttextActual = textActual + searchResultString(searchResults);\n\t\t}\n\t\t\n\t\tvar textWidth = measureText(textActual, bold);\n\t\t\n\t\tvar x = textX;\n\t\t\n\t\tif ( anchor == \'end\' )\n\t\t{\n\t\t\tx -= textWidth;\n\t\t}\n\t\telse if ( anchor != \'start\' )\n\t\t{\n\t\t\t// centered\n\t\t\tx -= textWidth / 2;\n\t\t}\n\t\t\n\t\tdrawBubble(angle, radius, textWidth, radial, flip);\n\t\t\n\t\tif ( searchResult )\n\t\t{\n\t\t\tdrawSearchHighlights\n\t\t\t(\n\t\t\t\ttextActual,\n\t\t\t\tx,\n\t\t\t\ttextY,\n\t\t\t\tangle,\n\t\t\t\ttrue\n\t\t\t)\n\t\t}\n\t}\n\t\n\tif ( searchResults )\n\t{\n\t\ttotalText = totalText + searchResultString(searchResults);\n\t}\n\t\n\tdrawText(totalText, textX, textY, angle, anchor, bold);\n\t\n\treturn flip;\n}\nfunction drawTick(start, length, angle)\n{\n\tif ( snapshotMode )\n\t{\n\t\tsvg +=\n\t\t\t\'<line x1="\' + (centerX + start) +\n\t\t\t\'" y1="\' + centerY +\n\t\t\t\'" x2="\' + (centerX + start + length) +\n\t\t\t\'" y2="\' + centerY +\n\t\t\t\'" class="tick" transform="rotate(\' +\n\t\t\tdegrees(angle) + \',\' + centerX + \',\' + centerY +\n\t\t\t\')"/>\';\n\t}\n\telse\n\t{\n\t\tcontext.rotate(angle);\n\t\tcontext.beginPath();\n\t\tcontext.moveTo(start, 0);\n\t\tcontext.lineTo(start + length, 0);\n\t\tcontext.lineWidth = thinLineWidth * 2;\n\t\tcontext.stroke();\n\t\tcontext.rotate(-angle);\n\t}\n}\nfunction drawWedge\n(\n\tangleStart,\n\tangleEnd,\n\tradiusInner,\n\tradiusOuter,\n\tcolor,\n\tpatternAlpha,\n\thighlight\n)\n{\n\tif ( context.globalAlpha == 0 )\n\t{\n\t\treturn;\n\t}\n\t\n\tif ( snapshotMode )\n\t{\n\t\tif ( angleEnd == angleStart + Math.PI * 2 )\n\t\t{\n\t\t\t// fudge to prevent overlap, which causes arc ambiguity\n\t\t\t//\n\t\t\tangleEnd -= .1 / gRadius;\n\t\t}\n\t\t\n\t\tvar longArc = angleEnd - angleStart > Math.PI ? 1 : 0;\n\t\t\n\t\tvar x1 = centerX + radiusInner * Math.cos(angleStart);\n\t\tvar y1 = centerY + radiusInner * Math.sin(angleStart);\n\t\t\n\t\tvar x2 = centerX + gRadius * Math.cos(angleStart);\n\t\tvar y2 = centerY + gRadius * Math.sin(angleStart);\n\t\t\n\t\tvar x3 = centerX + gRadius * Math.cos(angleEnd);\n\t\tvar y3 = centerY + gRadius * Math.sin(angleEnd);\n\t\t\n\t\tvar x4 = centerX + radiusInner * Math.cos(angleEnd);\n\t\tvar y4 = centerY + radiusInner * Math.sin(angleEnd);\n\t\t\n\t\tvar dArray =\n\t\t[\n\t\t\t" M ", x1, ",", y1,\n\t\t\t" L ", x2, ",", y2,\n\t\t\t" A ", gRadius, ",", gRadius, " 0 ", longArc, ",1 ", x3, ",", y3,\n\t\t\t" L ", x4, ",", y4,\n\t\t\t" A ", radiusInner, ",", radiusInner, " 0 ", longArc, " 0 ", x1, ",", y1,\n\t\t\t" Z "\n\t\t];\n\t\t\n\t\tsvg +=\n\t\t\t\'<path class="\'+ (highlight ? \'highlight\' : \'wedge\') + \'" fill="\' + color +\n\t\t\t\'" d="\' + dArray.join(\'\') + \'"/>\';\n\t\t\n\t\tif ( patternAlpha > 0 )\n\t\t{\n\t\t\tsvg +=\n\t\t\t\t\'<path class="wedge" fill="url(#hiddenPattern)" d="\' +\n\t\t\t\tdArray.join(\'\') + \'"/>\';\n\t\t}\n\t}\n\telse\n\t{\n\t\t// fudge to prevent seams during animation\n\t\t//\n\t\tangleEnd += 1 / gRadius;\n\t\t\n\t\tcontext.fillStyle = color;\n\t\tcontext.beginPath();\n\t\tcontext.arc(0, 0, radiusInner, angleStart, angleEnd, false);\n\t\tcontext.arc(0, 0, radiusOuter, angleEnd, angleStart, true);\n\t\tcontext.closePath();\n\t\tcontext.fill();\n\t\t\n\t\tif ( patternAlpha > 0 )\n\t\t{\n\t\t\tcontext.save();\n\t\t\tcontext.clip();\n\t\t\tcontext.globalAlpha = patternAlpha;\n\t\t\tcontext.fillStyle = hiddenPattern;\n\t\t\tcontext.fill();\n\t\t\tcontext.restore();\n\t\t}\n\t\t\n\t\tif ( highlight )\n\t\t{\n\t\t\tcontext.lineWidth = highlight ? highlightLineWidth : thinLineWidth;\n\t\t\tcontext.strokeStyle = \'black\';\n\t\t\tcontext.stroke();\n\t\t}\n\t}\n}\nfunction expand(node)\n{\n\tselectNode(node);\n\tupdateView();\n}\nfunction focusLost()\n{\n\tmouseX = -1;\n\tmouseY = -1;\n\tcheckHighlight();\n\tdocument.body.style.cursor = \'auto\';\n}\nfunction fontSizeDecrease()\n{\n\tif ( fontSize > 1 )\n\t{\n\t\tfontSize--;\n\t\tupdateViewNeeded = true;\n\t}\n}\nfunction fontSizeIncrease()\n{\n\tfontSize++;\n\tupdateViewNeeded = true;\n}\nfunction getGetString(name, value, bool)\n{\n\treturn name + \'=\' + (bool ? value ? \'true\' : \'false\' : value);\n}\nfunction hideLink()\n{\n\thide(linkText);\n\tshow(linkButton);\n}\nfunction show(object)\n{\n\tobject.style.display = \'inline\';\n}\nfunction hide(object)\n{\n\tobject.style.display = \'none\';\n}\nfunction showLink()\n{\n\tvar urlHalves = String(document.location).split(\'?\');\n\tvar newGetVariables = new Array();\n\t\n\tnewGetVariables.push\n\t(\n\t\tgetGetString(\'dataset\', currentDataset, false),\n\t\tgetGetString(\'node\', selectedNode.id, false),\n\t\tgetGetString(\'collapse\', collapse, true),\n\t\tgetGetString(\'color\', useHue(), true),\n\t\tgetGetString(\'depth\', maxAbsoluteDepth - 1, false),\n\t\tgetGetString(\'font\', fontSize, false),\n\t\tgetGetString(\'key\', showKeys, true)\n\t);\n\t\n\thide(linkButton);\n\tshow(linkText);\n\tlinkText.value = urlHalves[0] + \'?\' + getVariables.concat(newGetVariables).join(\'&\');\n\t//linkText.disabled = false;\n\tlinkText.focus();\n\tlinkText.select();\n\t//linkText.disabled = true;\n//\tdocument.location = urlHalves[0] + \'?\' + getVariables.join(\'&\');\n}\nfunction getFirstChild(element)\n{\n\telement = element.firstChild;\n\t\n\tif ( element && element.nodeType != 1 )\n\t{\n\t\telement = getNextSibling(element);\n\t}\n\t\n\treturn element;\n}\nfunction getNextSibling(element)\n{\n\tdo\n\t{\n\t\telement = element.nextSibling;\n\t}\n\twhile ( element && element.nodeType != 1 );\n\t\n\treturn element;\n}\nfunction getPercentage(fraction)\n{\n\treturn round(fraction * 100);\n}\nfunction hslText(hue)\n{\n\tif ( 1 || snapshotMode )\n\t{\n\t\t// Safari doesn\'t seem to allow hsl() in SVG\n\t\t\n\t\tvar rgb = hslToRgb(hue, saturation, (lightnessBase + lightnessMax) / 2);\n\t\t\n\t\treturn rgbText(rgb.r, rgb.g, rgb.b);\n\t}\n\telse\n\t{\n\t\tvar hslArray =\n\t\t[\n\t\t\t\'hsl(\',\n\t\t\tMath.floor(hue * 360),\n\t\t\t\',\',\n\t\t\tMath.floor(saturation * 100),\n\t\t\t\'%,\',\n\t\t\tMath.floor((lightnessBase + lightnessMax) * 50),\n\t\t\t\'%)\'\n\t\t];\n\t\t\n\t\treturn hslArray.join(\'\');\n\t}\n}\nfunction hslToRgb(h, s, l)\n{\n\tvar m1, m2;\n\tvar r, g, b;\n\t\n\tif (s == 0)\n\t{\n\t\tr = g = b = Math.floor((l * 255));\n\t}\n\telse\n\t{\n\t\tif (l <= 0.5)\n\t\t{\n\t\t\tm2 = l * (s + 1);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tm2 = l + s - l * s;\n\t\t}\n\t\t\n\t\tm1 = l * 2 - m2;\n\t\t\n\t\tr = Math.floor(hueToRgb(m1, m2, h + 1 / 3));\n\t\tg = Math.floor(hueToRgb(m1, m2, h));\n\t\tb = Math.floor(hueToRgb(m1, m2, h - 1/3));\n\t}\n\t\n\treturn {r: r, g: g, b: b};\n}\nfunction hueToRgb(m1, m2, hue)\n{\n\tvar v;\n\t\n\twhile (hue < 0)\n\t{\n\t\thue += 1;\n\t}\n\t\n\twhile (hue > 1)\n\t{\n\t\thue -= 1;\n\t}\n\t\n\tif (6 * hue < 1)\n\t\tv = m1 + (m2 - m1) * hue * 6;\n\telse if (2 * hue < 1)\n\t\tv = m2;\n\telse if (3 * hue < 2)\n\t\tv = m1 + (m2 - m1) * (2/3 - hue) * 6;\n\telse\n\t\tv = m1;\n\treturn 255 * v;\n}\nfunction interpolateHue(hueStart, hueEnd, valueStart, valueEnd)\n{\n\t// since the gradient will be RGB based, we need to add stops to hit all the\n\t// colors in the hue spectrum\n\t\n\thueStopPositions = new Array();\n\thueStopHsl = new Array();\n\thueStopText = new Array();\n\t\n\thueStopPositions.push(0);\n\thueStopHsl.push(hslText(hueStart));\n\thueStopText.push(round(valueStart));\n\t\n\tfor\n\t(\n\t\tvar i = (hueStart > hueEnd ? 5 / 6 : 1 / 6);\n\t\t(hueStart > hueEnd ? i > 0 : i < 1);\n\t\ti += (hueStart > hueEnd ? -1 : 1) / 6\n\t)\n\t{\n\t\tif\n\t\t(\n\t\t\thueStart > hueEnd ?\n\t\t\t\ti > hueEnd && i < hueStart :\n\t\t\t\ti > hueStart && i < hueEnd\n\t\t)\n\t\t{\n\t\t\thueStopPositions.push(lerp(i, hueStart, hueEnd, 0, 1));\n\t\t\thueStopHsl.push(hslText(i));\n\t\t\thueStopText.push(round(lerp\n\t\t\t(\n\t\t\t\ti,\n\t\t\t\thueStart,\n\t\t\t\thueEnd,\n\t\t\t\tvalueStart,\n\t\t\t\tvalueEnd\n\t\t\t)));\n\t\t}\n\t}\n\t\n\thueStopPositions.push(1);\n\thueStopHsl.push(hslText(hueEnd));\n\thueStopText.push(round(valueEnd));\n}\nfunction keyLineAngle(angle, keyAngle, bendRadius, keyX, keyY, pointsX, pointsY)\n{\n\tif ( angle < Math.PI / 2 && keyY < bendRadius * Math.sin(angle) \n\t|| angle > Math.PI / 2 && keyY < bendRadius)\n\t{\n\t\treturn Math.asin(keyY / bendRadius);\n\t}\n\telse\n\t{\n\t\t// find the angle of the normal to a tangent line that goes to\n\t\t// the label\n\t\t\n\t\tvar textDist = Math.sqrt\n\t\t(\n\t\t\tMath.pow(keyX, 2) +\n\t\t\tMath.pow(keyY, 2)\n\t\t);\n\t\t\n\t\tvar tanAngle = Math.acos(bendRadius / textDist) + keyAngle;\n\t\t\n\t\tif ( angle < tanAngle || angle < Math.PI / 2 )//|| labelLeft < centerX )\n\t\t{\n\t\t\t// angle doesn\'t reach far enough for tangent; collapse and\n\t\t\t// connect directly to label\n\t\t\t\n\t\t\tif ( keyY / Math.tan(angle) > 0 )\n\t\t\t{\n\t\t\t\tpointsX.push(keyY / Math.tan(angle));\n\t\t\t\tpointsY.push(keyY);\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tpointsX.push(bendRadius * Math.cos(angle));\n\t\t\t\tpointsY.push(bendRadius * Math.sin(angle));\n\t\t\t}\n\t\t\t\n\t\t\treturn angle;\n\t\t}\n\t\telse\n\t\t{\n\t\t\treturn tanAngle;\n\t\t}\n\t}\n}\nfunction keyOffset()\n{\n\treturn imageHeight - (keys - currentKey + 1) * (keySize + keyBuffer) + keyBuffer - margin;\n}\nfunction lerp(value, fromStart, fromEnd, toStart, toEnd)\n{\n\treturn (value - fromStart) *\n\t\t(toEnd - toStart) /\n\t\t(fromEnd - fromStart) +\n\t\ttoStart;\n}\nfunction createCanvas()\n{\n\tcanvas = document.createElement(\'canvas\');\n\tdocument.body.appendChild(canvas);\n\tcontext = canvas.getContext(\'2d\');\n}\nfunction load()\n{\n\tdocument.body.style.overflow = "hidden";\n\tdocument.body.style.margin = 0;\n\t\n\tcreateCanvas();\n\t\n\tif ( context == undefined )\n\t{\n\t\tdocument.body.innerHTML = \'\\\n<br/>This browser does not support HTML5 (see \\\n<a href="https://github.com/marbl/Krona/wiki/Browser%20support">Browser support</a>).\\\n\t\';\n\t\treturn;\n\t}\n\tif ( typeof context.fillText != \'function\' )\n\t{\n\t\tdocument.body.innerHTML = \'\\\n<br/>This browser does not support HTML5 canvas text (see \\\n<a href="https://github.com/marbl/Krona/wiki/Browser%20support">Browser support</a>).\\\n\t\';\n\t\treturn;\n\t}\n\t\n\tresize();\n\t\n\tvar kronaElement = document.getElementsByTagName(\'krona\')[0];\n\t\n\tvar magnitudeName;\n\tvar hueName;\n\tvar hueDefault;\n\tvar hueStart;\n\tvar hueEnd;\n\tvar valueStart;\n\tvar valueEnd;\n\t\n\tif ( kronaElement.getAttribute(\'collapse\') != undefined )\n\t{\n\t\tcollapse = kronaElement.getAttribute(\'collapse\') == \'true\';\n\t}\n\t\n\tif ( kronaElement.getAttribute(\'key\') != undefined )\n\t{\n\t\tshowKeys = kronaElement.getAttribute(\'key\') == \'true\';\n\t}\n\t\n\tfor\n\t(\n\t\tvar element = getFirstChild(kronaElement);\n\t\telement;\n\t\telement = getNextSibling(element)\n\t)\n\t{\n\t\tswitch ( element.tagName.toLowerCase() )\n\t\t{\n\t\t\tcase \'attributes\':\n\t\t\t\tmagnitudeName = element.getAttribute(\'magnitude\');\n\t\t\t\t//\n\t\t\t\tfor\n\t\t\t\t(\n\t\t\t\t\tvar attributeElement = getFirstChild(element);\n\t\t\t\t\tattributeElement;\n\t\t\t\t\tattributeElement = getNextSibling(attributeElement)\n\t\t\t\t)\n\t\t\t\t{\n\t\t\t\t\tvar tag = attributeElement.tagName.toLowerCase();\n\t\t\t\t\t\n\t\t\t\t\tif ( tag == \'attribute\' )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar attribute = new Attribute();\n\t\t\t\t\t\tattribute.name = attributeElement.firstChild.nodeValue.toLowerCase();\n\t\t\t\t\t\tattribute.displayName = attributeElement.getAttribute(\'display\');\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'hrefBase\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.hrefBase = attributeElement.getAttribute(\'hrefBase\');\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'target\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.target = attributeElement.getAttribute(\'target\');\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attribute.name == magnitudeName )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tmagnitudeIndex = attributes.length;\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'listAll\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.listAll = attributeElement.getAttribute(\'listAll\').toLowerCase();\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse if ( attributeElement.getAttribute(\'listNode\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.listNode = attributeElement.getAttribute(\'listNode\').toLowerCase();\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse if ( attributeElement.getAttribute(\'dataAll\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.dataAll = attributeElement.getAttribute(\'dataAll\').toLowerCase();\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse if ( attributeElement.getAttribute(\'dataNode\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.dataNode = attributeElement.getAttribute(\'dataNode\').toLowerCase();\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'postUrl\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.postUrl = attributeElement.getAttribute(\'postUrl\');\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'postVar\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.postVar = attributeElement.getAttribute(\'postVar\');\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributeElement.getAttribute(\'mono\') )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\tattribute.mono = true;\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tattributes.push(attribute);\n\t\t\t\t\t}\n\t\t\t\t\telse if ( tag == \'list\' )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar attribute = new Attribute();\n\t\t\t\t\t\t\n\t\t\t\t\t\tattribute.name = attributeElement.firstChild.nodeValue;\n\t\t\t\t\t\tattribute.list = true;\n\t\t\t\t\t\tattributes.push(attribute);\n\t\t\t\t\t}\n\t\t\t\t\telse if ( tag == \'data\' )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar attribute = new Attribute();\n\t\t\t\t\t\t\n\t\t\t\t\t\tattribute.name = attributeElement.firstChild.nodeValue;\n\t\t\t\t\t\tattribute.data = true;\n\t\t\t\t\t\tattributes.push(attribute);\n\t\t\t\t\t\t\n\t\t\t\t\t\tvar enableScript = document.createElement(\'script\');\n\t\t\t\t\t\tvar date = new Date();\n\t\t\t\t\t\tenableScript.src =\n\t\t\t\t\t\t\tattributeElement.getAttribute(\'enable\') + \'?\' +\n\t\t\t\t\t\t\tdate.getTime();\n\t\t\t\t\t\tdocument.body.appendChild(enableScript);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\tbreak;\n\t\t\t\n\t\t\tcase \'color\':\n\t\t\t\thueName = element.getAttribute(\'attribute\');\n\t\t\t\thueStart = Number(element.getAttribute(\'hueStart\')) / 360;\n\t\t\t\thueEnd = Number(element.getAttribute(\'hueEnd\')) / 360;\n\t\t\t\tvalueStart = Number(element.getAttribute(\'valueStart\'));\n\t\t\t\tvalueEnd = Number(element.getAttribute(\'valueEnd\'));\n\t\t\t\t//\n\t\t\t\tinterpolateHue(hueStart, hueEnd, valueStart, valueEnd);\n\t\t\t\t//\n\t\t\t\tif ( element.getAttribute(\'default\') == \'true\' )\n\t\t\t\t{\n\t\t\t\t\thueDefault = true;\n\t\t\t\t}\n\t\t\t\tbreak;\n\t\t\t\n\t\t\tcase \'datasets\':\n\t\t\t\tdatasetNames = new Array();\n\t\t\t\t//\n\t\t\t\tfor ( j = getFirstChild(element); j; j = getNextSibling(j) )\n\t\t\t\t{\n\t\t\t\t\tdatasetNames.push(j.firstChild.nodeValue);\n\t\t\t\t}\n\t\t\t\tdatasets = datasetNames.length;\n\t\t\t\tbreak;\n\t\t\t\n\t\t\tcase \'node\':\n\t\t\t\thead = loadTreeDOM\n\t\t\t\t(\n\t\t\t\t\telement,\n\t\t\t\t\tmagnitudeName,\n\t\t\t\t\thueName,\n\t\t\t\t\thueStart,\n\t\t\t\t\thueEnd,\n\t\t\t\t\tvalueStart,\n\t\t\t\t\tvalueEnd\n\t\t\t\t);\n\t\t\t\tbreak;\n\t\t}\n\t}\n\t\n\t// get GET options\n\t//\n\tvar urlHalves = String(document.location).split(\'?\');\n\tvar datasetDefault = 0;\n\tvar maxDepthDefault;\n\tvar nodeDefault = 0;\n\t//\n\tif ( urlHalves[1] )\n\t{\n\t\tvar vars = urlHalves[1].split(\'&\');\n\t\t\n\t\tfor ( i = 0; i < vars.length; i++ )\n\t\t{\n\t\t\tvar pair = vars[i].split(\'=\');\n\t\t\t\n\t\t\tswitch ( pair[0] )\n\t\t\t{\n\t\t\t\tcase \'collapse\':\n\t\t\t\t\tcollapse = pair[1] == \'true\';\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tcase \'color\':\n\t\t\t\t\thueDefault = pair[1] == \'true\';\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tcase \'dataset\':\n\t\t\t\t\tdatasetDefault = Number(pair[1]);\n\t\t\t\t\tbreak;\n\t\t\t\t\t\n\t\t\t\tcase \'depth\':\n\t\t\t\t\tmaxDepthDefault = Number(pair[1]) + 1;\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tcase \'key\':\n\t\t\t\t\tshowKeys = pair[1] == \'true\';\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tcase \'font\':\n\t\t\t\t\tfontSize = Number(pair[1]);\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tcase \'node\':\n\t\t\t\t\tnodeDefault = Number(pair[1]);\n\t\t\t\t\tbreak;\n\t\t\t\t\n\t\t\t\tdefault:\n\t\t\t\t\tgetVariables.push(pair[0] + \'=\' + pair[1]);\n\t\t\t\t\tbreak;\n\t\t\t}\n\t\t}\n\t}\n\t\n\taddOptionElements(hueName, hueDefault);\n\tsetCallBacks();\n\t\n\thead.sort();\n\tmaxAbsoluteDepth = 0;\n\tselectDataset(datasetDefault);\n\t\n\tif ( maxDepthDefault && maxDepthDefault < head.maxDepth )\n\t{\n\t\tmaxAbsoluteDepth = maxDepthDefault;\n\t}\n\telse\n\t{\n\t\tmaxAbsoluteDepth = head.maxDepth;\n\t}\n\t\n\tselectNode(nodes[nodeDefault]);\n\t\n\tsetInterval(update, 20);\n\t\n\twindow.onresize = handleResize;\n\tupdateMaxAbsoluteDepth();\n\tupdateViewNeeded = true;\n}\nfunction loadTreeDOM\n(\n\tdomNode,\n\tmagnitudeName,\n\thueName,\n\thueStart,\n\thueEnd,\n\tvalueStart,\n\tvalueEnd\n)\n{\n\tvar newNode = new Node();\n\t\n\tnewNode.name = domNode.getAttribute(\'name\');\n\t\n\tif ( domNode.getAttribute(\'href\') )\n\t{\n\t\tnewNode.href = domNode.getAttribute(\'href\');\n\t}\n\t\n\tif ( hueName )\n\t{\n\t\tnewNode.hues = new Array();\n\t}\n\t\n\tfor ( var i = getFirstChild(domNode); i; i = getNextSibling(i) )\n\t{\n\t\tswitch ( i.tagName.toLowerCase() )\n\t\t{\n\t\tcase \'node\': \n\t\t\tvar newChild = loadTreeDOM\n\t\t\t(\n\t\t\t\ti,\n\t\t\t\tmagnitudeName,\n\t\t\t\thueName,\n\t\t\t\thueStart,\n\t\t\t\thueEnd,\n\t\t\t\tvalueStart,\n\t\t\t\tvalueEnd\n\t\t\t);\n\t\t\tnewChild.parent = newNode;\n\t\t\tnewNode.children.push(newChild);\n\t\t\tbreak;\n\t\t\t\n\t\tdefault:\n\t\t\tvar attributeName = i.tagName.toLowerCase();\n\t\t\tvar index = attributeIndex(attributeName);\n\t\t\t//\n\t\t\tnewNode.attributes[index] = new Array();\n\t\t\t//\n\t\t\tfor ( var j = getFirstChild(i); j; j = getNextSibling(j) )\n\t\t\t{\n\t\t\t\tif ( attributes[index] == undefined )\n\t\t\t\t{\n\t\t\t\t\tvar x = 5;\n\t\t\t\t}\n\t\t\t\tif ( attributes[index].list )\n\t\t\t\t{\n\t\t\t\t\tnewNode.attributes[index].push(new Array());\n\t\t\t\t\t\n\t\t\t\t\tfor ( var k = getFirstChild(j); k; k = getNextSibling(k) )\n\t\t\t\t\t{\n\t\t\t\t\t\tnewNode.attributes[index][newNode.attributes[index].length - 1].push(k.firstChild.nodeValue);\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\telse\n\t\t\t\t{\n\t\t\t\t\tvar value = j.firstChild ? j.firstChild.nodeValue : \'\';\n\t\t\t\t\t\n\t\t\t\t\tif ( j.getAttribute(\'href\') )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar target;\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( attributes[index].target )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\ttarget = \' target="\' + attributes[index].target + \'"\';\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tvalue = \'<a href="\' + attributes[index].hrefBase + j.getAttribute(\'href\') + \'"\' + target + \'>\' + value + \'</a>\';\n\t\t\t\t\t}\n\t\t\t\t\t\n\t\t\t\t\tnewNode.attributes[index].push(value);\n\t\t\t\t}\n\t\t\t}\n\t\t\t//\n\t\t\tif ( attributeName == magnitudeName || attributeName == hueName )\n\t\t\t{\n\t\t\t\tfor ( j = 0; j < datasets; j++ )\n\t\t\t\t{\n\t\t\t\t\tvar value = newNode.attributes[index][j] == undefined ? 0 : Number(newNode.attributes[index][j]);\n\t\t\t\t\t\n\t\t\t\t\tnewNode.attributes[index][j] = value;\n\t\t\t\t\t\n\t\t\t\t\tif ( attributeName == hueName )\n\t\t\t\t\t{\n\t\t\t\t\t\tvar hue = lerp\n\t\t\t\t\t\t(\n\t\t\t\t\t\t\tvalue,\n\t\t\t\t\t\t\tvalueStart,\n\t\t\t\t\t\t\tvalueEnd,\n\t\t\t\t\t\t\thueStart,\n\t\t\t\t\t\t\thueEnd\n\t\t\t\t\t\t);\n\t\t\t\t\t\t\n\t\t\t\t\t\tif ( hue < hueStart == hueStart < hueEnd )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\thue = hueStart;\n\t\t\t\t\t\t}\n\t\t\t\t\t\telse if ( hue > hueEnd == hueStart < hueEnd )\n\t\t\t\t\t\t{\n\t\t\t\t\t\t\thue = hueEnd;\n\t\t\t\t\t\t}\n\t\t\t\t\t\t\n\t\t\t\t\t\tnewNode.hues[j] = hue;\n\t\t\t\t\t}\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\tif ( attributeName == hueName )\n\t\t\t\t{\n\t\t\t\t\tnewNode.hue = new Tween(newNode.hues[0], newNode.hues[0]);\n\t\t\t\t}\n\t\t\t}\n\t\t\tbreak;\n\t\t}\n\t}\n\t\n\treturn newNode;\n}\nfunction maxAbsoluteDepthDecrease()\n{\n\tif ( maxAbsoluteDepth > 2 )\n\t{\n\t\tmaxAbsoluteDepth--;\n\t\thead.setMaxDepths();\n\t\thandleResize();\n\t}\n}\nfunction maxAbsoluteDepthIncrease()\n{\n\tif ( maxAbsoluteDepth < head.maxDepth )\n\t{\n\t\tmaxAbsoluteDepth++;\n\t\thead.setMaxDepths();\n\t\thandleResize();\n\t}\n}\nfunction measureText(text, bold)\n{\n\tcontext.font = bold ? fontBold : fontNormal;\n\tvar dim = context.measureText(text);\n\treturn dim.width;\n}\nfunction min(a, b)\n{\n\treturn a < b ? a : b;\n}\nfunction minWidth()\n{\n\t// Min wedge width (at center) for displaying a node (or for displaying a\n\t// label if it\'s at the highest level being viewed, multiplied by 2 to make\n\t// further calculations simpler\n\t\n\treturn (fontSize * 2.3);\n}\nfunction mouseMove(e)\n{\n\tmouseX = e.pageX;\n\tmouseY = e.pageY - headerHeight;\n\tmouseXRel = (mouseX - centerX) * backingScale()\n\tmouseYRel = (mouseY - centerY) * backingScale()\n\t\n\tif ( head && ! quickLook )\n\t{\n\t\tcheckHighlight();\n\t}\n}\nfunction mouseClick(e)\n{\n\tif ( highlightedNode == focusNode && focusNode != selectedNode || selectedNode.hasParent(highlightedNode) )\n\t{\n\t\tif ( highlightedNode.hasChildren() )\n\t\t{\n\t\t\texpand(highlightedNode);\n\t\t}\n\t}\n\telse if ( progress == 1 )//( highlightedNode != selectedNode )\n\t{\n\t\tsetFocus(highlightedNode);\n//\t\tdocument.body.style.cursor=\'ew-resize\';\n\t\tdraw();\n\t\tcheckHighlight();\n\t\tvar date = new Date();\n\t\tmouseDownTime = date.getTime();\n\t\tmouseDown = true;\n\t}\n}\nfunction mouseUp(e)\n{\n\tif ( quickLook )\n\t{\n\t\tnavigateBack();\n\t\tquickLook = false;\n\t}\n\t\n\tmouseDown = false;\n}\nfunction navigateBack()\n{\n\tif ( nodeHistoryPosition > 0 )\n\t{\n\t\tnodeHistory[nodeHistoryPosition] = selectedNode;\n\t\tnodeHistoryPosition--;\n\t\t\n\t\tif ( nodeHistory[nodeHistoryPosition].collapse )\n\t\t{\n\t\t\tcollapseCheckBox.checked = collapse = false;\n\t\t}\n\t\t\n\t\tsetSelectedNode(nodeHistory[nodeHistoryPosition]);\n\t\tupdateDatasetButtons();\n\t\tupdateView();\n\t}\n}\nfunction navigateUp()\n{\n\tif ( selectedNode.getParent() )\n\t{\n\t\tselectNode(selectedNode.getParent());\n\t\tupdateView();\n\t}\n}\nfunction navigateForward()\n{\n\tif ( nodeHistoryPosition < nodeHistory.length - 1 )\n\t{\n\t\tnodeHistoryPosition++;\n\t\tvar newNode = nodeHistory[nodeHistoryPosition];\n\t\t\n\t\tif ( newNode.collapse )\n\t\t{\n\t\t\tcollapseCheckBox.checked = collapse = false;\n\t\t}\n\t\t\n\t\tif ( nodeHistoryPosition == nodeHistory.length - 1 )\n\t\t{\n\t\t\t// this will ensure the forward button is disabled\n\t\t\t\n\t\t\tnodeHistory.length = nodeHistoryPosition;\n\t\t}\n\t\t\n\t\tsetSelectedNode(newNode);\n\t\tupdateDatasetButtons();\n\t\tupdateView();\n\t}\n}\nfunction nextDataset()\n{\n\tvar newDataset = currentDataset;\n\t\n\tdo\n\t{\n\t\tif ( newDataset == datasets - 1 )\n\t\t{\n\t\t\tnewDataset = 0;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tnewDataset++;\n\t\t}\n\t}\n\twhile ( datasetDropDown.options[newDataset].disabled )\n\t\n\tselectDataset(newDataset);\n}\nfunction onDatasetChange()\n{\n\tselectDataset(datasetDropDown.selectedIndex);\n}\nfunction onKeyDown(event)\n{\n\tif\n\t(\n\t\tevent.keyCode == 37 &&\n\t\tdocument.activeElement.id != \'search\' &&\n\t\tdocument.activeElement.id != \'linkText\'\n\t)\n\t{\n\t\tnavigateBack();\n\t\tevent.preventDefault();\n\t}\n\telse if\n\t(\n\t\tevent.keyCode == 39 &&\n\t\tdocument.activeElement.id != \'search\' &&\n\t\tdocument.activeElement.id != \'linkText\'\n\t)\n\t{\n\t\tnavigateForward();\n\t\tevent.preventDefault();\n\t}\n\telse if ( event.keyCode == 38 && datasets > 1 )\n\t{\n\t\tprevDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n\telse if ( event.keyCode == 40 && datasets > 1 )\n\t{\n\t\tnextDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n\telse if ( event.keyCode == 9 && datasets > 1 )\n\t{\n\t\tselectLastDataset();\n\t\tevent.preventDefault();\n\t}\n\telse if ( event.keyCode == 83 )\n\t{\n\t\tprogress += .2;\n\t}\n\telse if ( event.keyCode == 66 )\n\t{\n\t\tprogress -= .2;\n\t}\n\telse if ( event.keyCode == 70 )\n\t{\n\t\tprogress = 1;\n\t}\n}\nfunction onKeyPress(event)\n{\n\tif ( event.keyCode == 38 && datasets > 1 )\n\t{\n//\t\tprevDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n\telse if ( event.keyCode == 40 && datasets > 1 )\n\t{\n//\t\tnextDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n}\nfunction onKeyUp(event)\n{\n\tif ( event.keyCode == 27 && document.activeElement.id == \'search\' )\n\t{\n\t\tsearch.value = \'\';\n\t\tonSearchChange();\n\t}\n\telse if ( event.keyCode == 38 && datasets > 1 )\n\t{\n//\t\tprevDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n\telse if ( event.keyCode == 40 && datasets > 1 )\n\t{\n//\t\tnextDataset();\n\t\t\n\t\t//if ( document.activeElement.id == \'datasets\' )\n\t\t{\n\t\t\tevent.preventDefault();\n\t\t}\n\t}\n}\nfunction onSearchChange()\n{\n\tnSearchResults = 0;\n\thead.search();\n\t\n\tif ( search.value == \'\' )\n\t{\n\t\tsearchResults.innerHTML = \'\';\n\t}\n\telse\n\t{\n\t\tsearchResults.innerHTML = nSearchResults + \' results\';\n\t}\n\t\n\tsetFocus(selectedNode);\n\tdraw();\n}\nfunction post(url, variable, value, postWindow)\n{\n\tvar form = document.createElement(\'form\');\n\tvar input = document.createElement(\'input\');\n\tvar inputDataset = document.createElement(\'input\');\n\t\n\tform.appendChild(input);\n\tform.appendChild(inputDataset);\n\t\n\tform.method = "POST";\n\tform.action = url;\n\t\n\tif ( postWindow == undefined )\n\t{\n\t\tform.target = \'_blank\';\n\t\tpostWindow = window;\n\t}\n\t\n\tinput.type = \'hidden\';\n\tinput.name = variable;\n\tinput.value = value;\n\t\n\tinputDataset.type = \'hidden\';\n\tinputDataset.name = \'dataset\';\n\tinputDataset.value = currentDataset;\n\t\n\tpostWindow.document.body.appendChild(form);\n\tform.submit();\n}\nfunction prevDataset()\n{\n\tvar newDataset = currentDataset;\n\t\n\tdo\n\t{\n\t\tif ( newDataset == 0 )\n\t\t{\n\t\t\tnewDataset = datasets - 1;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tnewDataset--;\n\t\t}\n\t}\n\twhile ( datasetDropDown.options[newDataset].disabled );\n\t\n\tselectDataset(newDataset);\n}\nfunction radiusDecrease()\n{\n\tif ( bufferFactor < .309 )\n\t{\n\t\tbufferFactor += .03;\n\t\tupdateViewNeeded = true;\n\t}\n}\nfunction radiusIncrease()\n{\n\tif ( bufferFactor > .041 )\n\t{\n\t\tbufferFactor -= .03;\n\t\tupdateViewNeeded = true;\n\t}\n}\nfunction resetKeyOffset()\n{\n\tcurrentKey = 1;\n\tkeyMinTextLeft = centerX + gRadius + buffer - buffer / (keys + 1) / 2 + fontSize / 2;\n\tkeyMinAngle = 0;\n}\nfunction rgbText(r, g, b)\n{\n\tvar rgbArray =\n\t[\n\t\t"rgb(",\n\t\tMath.floor(r),\n\t\t",",\n\t\tMath.floor(g),\n\t\t",",\n\t\tMath.floor(b),\n\t\t")"\n\t];\n\t\n\treturn rgbArray.join(\'\');\n}\nfunction round(number)\n{\n\tif ( number >= 1 || number <= -1 )\n\t{\n\t\treturn number.toFixed(0);\n\t}\n\telse\n\t{\n\t\treturn number.toPrecision(1);\n\t}\n}\nfunction roundedRectangle(x, y, width, height, radius)\n{\n\tif ( radius * 2 > width )\n\t{\n\t\tradius = width / 2;\n\t}\n\t\n\tif ( radius * 2 > height )\n\t{\n\t\tradius = height / 2;\n\t}\n\t\n\tcontext.beginPath();\n\tcontext.arc(x + radius, y + radius, radius, Math.PI, Math.PI * 3 / 2, false);\n\tcontext.lineTo(x + width - radius, y);\n\tcontext.arc(x + width - radius, y + radius, radius, Math.PI * 3 / 2, Math.PI * 2, false);\n\tcontext.lineTo(x + width, y + height - radius);\n\tcontext.arc(x + width - radius, y + height - radius, radius, 0, Math.PI / 2, false);\n\tcontext.lineTo(x + radius, y + height);\n\tcontext.arc(x + radius, y + height - radius, radius, Math.PI / 2, Math.PI, false);\n\tcontext.lineTo(x, y + radius);\n}\nfunction passClick(e)\n{\n\tmouseClick(e);\n}\nfunction searchResultString(results)\n{\n\tvar searchResults = this.searchResults;\n\t\n\tif ( this.isSearchResult )\n\t{\n\t\t// don\'t count ourselves\n\t\tsearchResults--;\n\t}\n\t\n\treturn \' - \' + results + (results > 1 ? \' results\' : \' result\');\n}\nfunction setCallBacks()\n{\n\tcanvas.onselectstart = function(){return false;} // prevent unwanted highlighting\n\toptions.onselectstart = function(){return false;} // prevent unwanted highlighting\n\tdocument.onmousemove = mouseMove;\n\twindow.onblur = focusLost;\n\twindow.onmouseout = focusLost;\n\tdocument.onkeyup = onKeyUp;\n\tdocument.onkeydown = onKeyDown;\n\tcanvas.onmousedown = mouseClick;\n\tdocument.onmouseup = mouseUp;\n\tkeyControl.onclick = toggleKeys;\n\tcollapseCheckBox = document.getElementById(\'collapse\');\n\tcollapseCheckBox.checked = collapse;\n\tcollapseCheckBox.onclick = handleResize;\n\tcollapseCheckBox.onmousedown = suppressEvent;\n\tmaxAbsoluteDepthText = document.getElementById(\'maxAbsoluteDepth\');\n\tmaxAbsoluteDepthButtonDecrease = document.getElementById(\'maxAbsoluteDepthDecrease\');\n\tmaxAbsoluteDepthButtonIncrease = document.getElementById(\'maxAbsoluteDepthIncrease\');\n\tmaxAbsoluteDepthButtonDecrease.onclick = maxAbsoluteDepthDecrease;\n\tmaxAbsoluteDepthButtonIncrease.onclick = maxAbsoluteDepthIncrease;\n\tmaxAbsoluteDepthButtonDecrease.onmousedown = suppressEvent;\n\tmaxAbsoluteDepthButtonIncrease.onmousedown = suppressEvent;\n\tfontSizeText = document.getElementById(\'fontSize\');\n\tfontSizeButtonDecrease = document.getElementById(\'fontSizeDecrease\');\n\tfontSizeButtonIncrease = document.getElementById(\'fontSizeIncrease\');\n\tfontSizeButtonDecrease.onclick = fontSizeDecrease;\n\tfontSizeButtonIncrease.onclick = fontSizeIncrease;\n\tfontSizeButtonDecrease.onmousedown = suppressEvent;\n\tfontSizeButtonIncrease.onmousedown = suppressEvent;\n\tradiusButtonDecrease = document.getElementById(\'radiusDecrease\');\n\tradiusButtonIncrease = document.getElementById(\'radiusIncrease\');\n\tradiusButtonDecrease.onclick = radiusDecrease;\n\tradiusButtonIncrease.onclick = radiusIncrease;\n\tradiusButtonDecrease.onmousedown = suppressEvent;\n\tradiusButtonIncrease.onmousedown = suppressEvent;\n\tmaxAbsoluteDepth = 0;\n\tbackButton = document.getElementById(\'back\');\n\tbackButton.onclick = navigateBack;\n\tbackButton.onmousedown = suppressEvent;\n\tforwardButton = document.getElementById(\'forward\');\n\tforwardButton.onclick = navigateForward;\n\tforwardButton.onmousedown = suppressEvent;\n\tsnapshotButton = document.getElementById(\'snapshot\');\n\tsnapshotButton.onclick = snapshot;\n\tsnapshotButton.onmousedown = suppressEvent;\n\tdetailsName = document.getElementById(\'detailsName\');\n\tdetailsExpand = document.getElementById(\'detailsExpand\');\n\tdetailsInfo = document.getElementById(\'detailsInfo\');\n\tsearch = document.getElementById(\'search\');\n\tsearch.onkeyup = onSearchChange;\n\tsearch.onmousedown = suppressEvent;\n\tsearchResults = document.getElementById(\'searchResults\');\n\tuseHueDiv = document.getElementById(\'useHueDiv\');\n\tlinkButton = document.getElementById(\'linkButton\');\n\tlinkButton.onclick = showLink;\n\tlinkButton.onmousedown = suppressEvent;\n\tlinkText = document.getElementById(\'linkText\');\n\tlinkText.onblur = hideLink;\n\tlinkText.onmousedown = suppressEvent;\n\thide(linkText);\n\tvar helpButton = document.getElementById(\'help\');\n\thelpButton.onmousedown = suppressEvent;\n\tvar searchClear = document.getElementById(\'searchClear\');\n\tsearchClear.onmousedown = suppressEvent;\n\tif ( datasets > 1 )\n\t{\n\t\tdatasetDropDown.onmousedown = suppressEvent;\n\t\tvar prevDatasetButton = document.getElementById(\'prevDataset\');\n\t\tprevDatasetButton.onmousedown = suppressEvent;\n\t\tvar nextDatasetButton = document.getElementById(\'nextDataset\');\n\t\tnextDatasetButton.onmousedown = suppressEvent;\n\t\tvar lastDatasetButton = document.getElementById(\'lastDataset\');\n\t\tlastDatasetButton.onmousedown = suppressEvent;\n\t}\n\t\n\timage = document.getElementById(\'hiddenImage\');\n\t\n\tif ( image.complete )\n\t{\n\t\thiddenPattern = context.createPattern(image, \'repeat\');\n\t}\n\telse\n\t{\n\t\timage.onload = function()\n\t\t{\n\t\t\thiddenPattern = context.createPattern(image, \'repeat\');\n\t\t}\n\t}\n\t\n\tvar loadingImageElement = document.getElementById(\'loadingImage\');\n\t\n\tif ( loadingImageElement )\n\t{\n\t\tloadingImage = loadingImageElement.src;\n\t}\n}\nfunction selectDataset(newDataset)\n{\n\tlastDataset = currentDataset;\n\tcurrentDataset = newDataset\n\tif ( datasets > 1 )\n\t{\n\t\tdatasetDropDown.selectedIndex = currentDataset;\n\t\tupdateDatasetButtons();\n\t\tdatasetAlpha.start = 1.5;\n\t\tdatasetChanged = true;\n\t}\n\thead.setMagnitudes(0);\n\thead.setDepth(1, 1);\n\thead.setMaxDepths();\n\thandleResize();\n}\nfunction selectLastDataset()\n{\n\tselectDataset(lastDataset);\n\thandleResize();\n}\nfunction selectNode(newNode)\n{\n\tif ( selectedNode != newNode )\n\t{\n\t\t// truncate history at current location to create a new branch\n\t\t//\n\t\tnodeHistory.length = nodeHistoryPosition;\n\t\t\n\t\tif ( selectedNode != 0 )\n\t\t{\n\t\t\tnodeHistory.push(selectedNode);\n\t\t\tnodeHistoryPosition++;\n\t\t}\n\t\t\n\t\tsetSelectedNode(newNode);\n\t\t//updateView();\n\t}\n\t\n\tupdateDatasetButtons();\n}\nfunction setFocus(node)\n{\n\tif ( node == focusNode )\n\t{\n//\t\treturn;\n\t}\n\t\n\tfocusNode = node;\n\t\n\tif ( node.href )\n\t{\n\t\tdetailsName.innerHTML =\n\t\t\t\'<a target="_blank" href="\' + node.href + \'">\' + node.name + \'</a>\';\n\t}\n\telse\n\t{\n\t\tdetailsName.innerHTML = node.name;\n\t}\n\t\n\tvar table = \'<table>\';\n\t\n\ttable += \'<tr><td></td></tr>\';\n\t\n\tfor ( var i = 0; i < node.attributes.length; i++ )\n\t{\n\t\tif ( attributes[i].displayName && node.attributes[i] != undefined )\n\t\t{\n\t\t\tvar index = node.attributes[i].length == 1 && attributes[i].mono ? 0 : currentDataset;\n\t\t\t\n\t\t\tif ( typeof node.attributes[i][currentDataset] == \'number\' || node.attributes[i][index] != undefined && node.attributes[i][currentDataset] != \'\' )\n\t\t\t{\n\t\t\t\tvar value = node.attributes[i][index];\n\t\t\t\t\n\t\t\t\tif ( attributes[i].listNode != undefined )\n\t\t\t\t{\n\t\t\t\t\tvalue =\n\t\t\t\t\t\t\'<a href="" onclick="showList(\' +\n\t\t\t\t\t\tattributeIndex(attributes[i].listNode) + \',\' + i +\n\t\t\t\t\t\t\',false);return false;" title="Show list">\' +\n\t\t\t\t\t\tvalue + \'</a>\';\n\t\t\t\t}\n\t\t\t\telse if ( attributes[i].listAll != undefined )\n\t\t\t\t{\n\t\t\t\t\tvalue =\n\t\t\t\t\t\t\'<a href="" onclick="showList(\' +\n\t\t\t\t\t\tattributeIndex(attributes[i].listAll) + \',\' + i +\n\t\t\t\t\t\t\',true);return false;" title="Show list">\' +\n\t\t\t\t\t\tvalue + \'</a>\';\n\t\t\t\t}\n\t\t\t\telse if ( attributes[i].dataNode != undefined && dataEnabled )\n\t\t\t\t{\n\t\t\t\t\tvalue =\n\t\t\t\t\t\t\'<a href="" onclick="showData(\' +\n\t\t\t\t\t\tattributeIndex(attributes[i].dataNode) + \',\' + i +\n\t\t\t\t\t\t\',false);return false;" title="Show data">\' +\n\t\t\t\t\t\tvalue + \'</a>\';\n\t\t\t\t}\n\t\t\t\telse if ( attributes[i].dataAll != undefined && dataEnabled )\n\t\t\t\t{\n\t\t\t\t\tvalue =\n\t\t\t\t\t\t\'<a href="" onclick="showData(\' +\n\t\t\t\t\t\tattributeIndex(attributes[i].dataAll) + \',\' + i +\n\t\t\t\t\t\t\',true);return false;" title="Show data">\' +\n\t\t\t\t\t\tvalue + \'</a>\';\n\t\t\t\t}\n\t\t\t\t\n\t\t\t\ttable +=\n\t\t\t\t\t\'<tr><td><strong>\' + attributes[i].displayName + \':</strong></td><td>\' +\n\t\t\t\t\tvalue + \'</td></tr>\';\n\t\t\t}\n\t\t}\n\t}\n\t\n\ttable += \'</table>\';\n\tdetailsInfo.innerHTML = table;\n\t\n\tdetailsExpand.disabled = !focusNode.hasChildren() || focusNode == selectedNode;\n}\nfunction setSelectedNode(newNode)\n{\n\tif ( selectedNode && selectedNode.hasParent(newNode) )\n\t{\n\t\tzoomOut = true;\n\t}\n\telse\n\t{\n\t\tzoomOut = false;\n\t}\n\t\n\tselectedNodeLast = selectedNode;\n\tselectedNode = newNode;\n\t\n\t//if ( focusNode != selectedNode )\n\t{\n\t\tsetFocus(selectedNode);\n\t}\n}\nfunction waitForData(dataWindow, target, title, time, postUrl, postVar)\n{\n\tif ( nodeData.length == target )\n\t{\n\t\tif ( postUrl != undefined )\n\t\t{\n\t\t\tfor ( var i = 0; i < nodeData.length; i++ )\n\t\t\t{\n\t\t\t\tnodeData[i] = nodeData[i].replace(/\\n/g, \',\');\n\t\t\t}\n\t\t\t\n\t\t\tvar postString = nodeData.join(\'\');\n\t\t\tpostString = postString.slice(0, -1);\n\t\t\t\n\t\t\tdataWindow.document.body.removeChild(dataWindow.document.getElementById(\'loading\'));\n\t\t\tdocument.body.removeChild(document.getElementById(\'data\'));\n\t\t\t\n\t\t\tpost(postUrl, postVar, postString, dataWindow);\n\t\t}\n\t\telse\n\t\t{\n\t\t\t//dataWindow.document.body.removeChild(dataWindow.document.getElementById(\'loading\'));\n\t\t\t//document.body.removeChild(document.getElementById(\'data\'));\n\t\t\t\n\t\t\tdataWindow.document.open();\n\t\t\tdataWindow.document.write(\'<pre>\' + nodeData.join(\'\') + \'</pre>\');\n\t\t\tdataWindow.document.close();\n\t\t}\n\t\t\n\t\tdataWindow.document.title = title; // replace after document.write()\n\t}\n\telse\n\t{\n\t\tvar date = new Date();\n\t\t\n\t\tif ( date.getTime() - time > 10000 )\n\t\t{\n\t\t\tdataWindow.document.body.removeChild(dataWindow.document.getElementById(\'loading\'));\n\t\t\tdocument.body.removeChild(document.getElementById(\'data\'));\n\t\t\tdataWindow.document.body.innerHTML =\n\t\t\t\t\'Timed out loading supplemental files for:<br/>\' + document.location;\n\t\t}\n\t\telse\n\t\t{\n\t\t\tsetTimeout(function() {waitForData(dataWindow, target, title, time, postUrl, postVar);}, 100);\n\t\t}\n\t}\n}\nfunction data(newData)\n{\n\tnodeData.push(newData);\n}\nfunction enableData()\n{\n\tdataEnabled = true;\n}\nfunction showData(indexData, indexAttribute, summary)\n{\n\tvar dataWindow = window.open(\'\', \'_blank\');\n\tvar title = \'Krona - \' + attributes[indexAttribute].displayName + \' - \' + focusNode.name;\n\tdataWindow.document.title = title;\n\t\n\tnodeData = new Array();\n\t\n\tif ( dataWindow && dataWindow.document && dataWindow.document.body != null )\n\t{\n\t\t//var loadImage = document.createElement(\'img\');\n\t\t//loadImage.src = "file://localhost/Users/ondovb/Krona/KronaTools/img/loading.gif";\n\t\t//loadImage.id = "loading";\n\t\t//loadImage.alt = "Loading...";\n\t\t//dataWindow.document.body.appendChild(loadImage);\n\t\tdataWindow.document.body.innerHTML =\n\t\t\t\'<img id="loading" src="\' + loadingImage + \'" alt="Loading..."></img>\';\n\t}\n\t\n\tvar scripts = document.createElement(\'div\');\n\tscripts.id = \'data\';\n\tdocument.body.appendChild(scripts);\n\t\n\tvar files = focusNode.getData(indexData, summary);\n\t\n\tvar date = new Date();\n\tvar time = date.getTime();\n\t\n\tfor ( var i = 0; i < files.length; i++ )\n\t{\n\t\tvar script = document.createElement(\'script\');\n\t\tscript.src = files[i] + \'?\' + time;\n\t\tscripts.appendChild(script);\n\t}\n\t\n\twaitForData(dataWindow, files.length, title, time, attributes[indexAttribute].postUrl, attributes[indexAttribute].postVar);\n\t\n\treturn false;\n}\nfunction showList(indexList, indexAttribute, summary)\n{\n\tvar list = focusNode.getList(indexList, summary);\n\t\n\tif ( attributes[indexAttribute].postUrl != undefined )\n\t{\n\t\tpost(attributes[indexAttribute].postUrl, attributes[indexAttribute].postVar, list.join(\',\'));\n\t}\n\telse\n\t{\n\t\tvar dataWindow = window.open(\'\', \'_blank\');\n\t\t\n\t\tif ( true || navigator.appName == \'Microsoft Internet Explorer\' ) // :(\n\t\t{\n\t\t\tdataWindow.document.open();\n\t\t\tdataWindow.document.write(\'<pre>\' + list.join(\'\\n\') + \'</pre>\');\n\t\t\tdataWindow.document.close();\n\t\t}\n\t\telse\n\t\t{\n\t\t\tvar pre = document.createElement(\'pre\');\n\t\t\tdataWindow.document.body.appendChild(pre);\n\t\t\tpre.innerHTML = list;\n\t\t}\n\t\t\n\t\tdataWindow.document.title = \'Krona - \' + attributes[indexAttribute].displayName + \' - \' + focusNode.name;\n\t}\n}\nfunction snapshot()\n{\n\tsvg = svgHeader();\n\t\n\tresetKeyOffset();\n\t\n\tsnapshotMode = true;\n\t\n\tselectedNode.draw(false, true);\n\tselectedNode.draw(true, true);\n\t\n\tif ( focusNode != 0 && focusNode != selectedNode )\n\t{\n\t\tcontext.globalAlpha = 1;\n\t\tfocusNode.drawHighlight(true);\n\t}\n\t\n\tif ( hueDisplayName && useHue() )\n\t{\n\t\tdrawLegendSVG();\n\t}\n\t\n\tsnapshotMode = false;\n\t\n\tsvg += svgFooter();\n\t\n\tsnapshotWindow = window.open\n\t(\n\t\t\'data:image/svg+xml;charset=utf-8,\' + encodeURIComponent(svg),\n\t\t\'_blank\'\n\t);\n/*\tvar data = window.open(\'data:text/plain;charset=utf-8,hello\', \'_blank\');\n\tvar data = window.open(\'\', \'_blank\');\n\tdata.document.open(\'text/plain\');\n\tdata.document.write(\'hello\');\n\tdata.document.close();\n\tvar button = document.createElement(\'input\');\n\tbutton.type = \'button\';\n\tbutton.value = \'save\';\n\tbutton.onclick = save;\n\tdata.document.body.appendChild(button);\n//\tsnapshotWindow.document.write(svg);\n//\tsnapshotWindow.document.close();\n*/\t\n}\nfunction save()\n{\n\talert(document.body.innerHTML);\n}\nfunction spacer()\n{\n\tif ( snapshotMode )\n\t{\n\t\treturn \'&#160;&#160;&#160;\';\n\t}\n\telse\n\t{\n\t\treturn \'   \';\n\t}\n}\nfunction suppressEvent(e)\n{\n\te.cancelBubble = true;\n\tif (e.stopPropagation) e.stopPropagation();\n}\nfunction svgFooter()\n{\n\treturn \'</svg>\';\n}\nfunction svgHeader()\n{\n\tvar patternWidth = fontSize * .6;//radius / 50;\n\t\n\treturn \'\\\n<?xml version="1.0" standalone="no"?>\\\n<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" \\\n\t"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\\\n<svg width="\' + imageWidth + \'" height="\' + imageHeight + \'" version="1.1"\\\n\txmlns="http://www.w3.org/2000/svg">\\\n<title>Krona (snapshot) - \' +\n(datasets > 1 ? datasetNames[currentDataset] + \' - \' : \'\') + selectedNode.name +\n\'</title>\\\n<defs>\\\n\t<style type="text/css">\\\n\ttext {font-size: \' + fontSize + \'px; font-family: \' + fontFamily + \'; dominant-baseline:central}\\\n\tpath {stroke-width:\' + thinLineWidth * fontSize / 12 + \';}\\\n\tpath.wedge {stroke:none}\\\n\tpath.line {fill:none;stroke:black;}\\\n\tline {stroke:black;stroke-width:\' + thinLineWidth * fontSize / 12 + \';}\\\n\tline.tick {stroke-width:\' + thinLineWidth * fontSize / 6 + \';}\\\n\tline.pattern {stroke-width:\' + thinLineWidth * fontSize / 18 + \';}\\\n\tcircle {fill:none;stroke:black;stroke-width:\' + thinLineWidth * fontSize / 12 + \';}\\\n\trect {stroke:black;stroke-width:\' + thinLineWidth * fontSize / 12 + \';}\\\n\t.highlight {stroke:black;stroke-width:\'+ highlightLineWidth * fontSize / 12 + \';}\\\n\t.searchHighlight {fill:rgb(255, 255, 100);stroke:none;}\\\n\t</style>\\\n<pattern id="hiddenPattern" patternUnits="userSpaceOnUse" \\\nx="0" y="0" width="\' + patternWidth + \'" height="\' + patternWidth + \'">\\\n<line class="pattern" x1="0" y1="0" x2="\' + patternWidth / 2 + \'" y2="\' + patternWidth / 2 + \'"/>\\\n<line class="pattern" x1="\' + patternWidth / 2 + \'" y1="\' + patternWidth +\n\'" x2="\' + patternWidth + \'" y2="\' + patternWidth / 2 + \'"/>\\\n</pattern>\\\n</defs>\\\n\';\n}\nfunction svgText(text, x, y, anchor, bold, color)\n{\n\tif ( typeof(anchor) == \'undefined\' )\n\t{\n\t\tanchor = \'start\';\n\t}\n\t\n\tif ( color == undefined )\n\t{\n\t\tcolor = \'black\';\n\t}\n\t\n\treturn \'<text x="\' + x + \'" y="\' + y +\n\t\t\'" style="font-color:\' + color + \';font-weight:\' + (bold ? \'bold\' : \'normal\') +\n\t\t\'" text-anchor="\' + anchor + \'">\' + text + \'</text>\';\n}\nfunction toggleKeys()\n{\n\tif ( showKeys )\n\t{\n\t\tkeyControl.value = \'...\';\n\t\tshowKeys = false;\n\t}\n\telse\n\t{\n\t\tkeyControl.value = \'x\';\n\t\tshowKeys = true;\n\t}\n\t\n\tupdateKeyControl();\n\t\n\tif ( progress == 1 )\n\t{\n\t\tdraw();\n\t}\n}\nfunction update()\n{\n\tif ( ! head )\n\t{\n\t\treturn;\n\t}\n\t\n\tif ( mouseDown && focusNode != selectedNode )\n\t{\n\t\tvar date = new Date();\n\t\t\n\t\tif ( date.getTime() - mouseDownTime > quickLookHoldLength )\n\t\t{\n\t\t\tif ( focusNode.hasChildren() )\n\t\t\t{\n\t\t\t\texpand(focusNode);\n\t\t\t\tquickLook = true;\n\t\t\t}\n\t\t}\n\t}\n\t\n\tif ( updateViewNeeded )\n\t{\n\t\tresize();\n\t\tmouseX = -1;\n\t\tmouseY = -1;\n\t\t\n\t\tcollapse = collapseCheckBox.checked;\n\t\tcompress = true;//compressCheckBox.checked;\n\t\tshorten = true;//shortenCheckBox.checked;\n\t\t\n\t\tcheckSelectedCollapse();\n\t\tupdateMaxAbsoluteDepth();\n\t\t\n\t\tif ( focusNode.getCollapse() || focusNode.depth > maxAbsoluteDepth )\n\t\t{\n\t\t\tsetFocus(selectedNode);\n\t\t}\n\t\telse\n\t\t{\n\t\t\tsetFocus(focusNode);\n\t\t}\n\t\t\n\t\tupdateView();\n\t\t\n\t\tupdateViewNeeded = false;\n\t}\n\t\n\tvar date = new Date();\n\tprogress = (date.getTime() - tweenStartTime) / tweenLength;\n//\tprogress += .01;\n\t\n\tif ( progress >= 1 )\n\t{\n\t\tprogress = 1;\n\t}\n\t\n\tif ( progress != progressLast )\n\t{\n\t\ttweenFactor =// progress;\n\t\t\t(1 / (1 + Math.exp(-tweenCurvature * (progress - .5))) - .5) /\n\t\t\t(tweenMax - .5) / 2 + .5;\n\t\t\n\t\tif ( progress == 1 )\n\t\t{\n\t\t\tsnapshotButton.disabled = false;\n\t\t\tzoomOut = false;\n\t\t\t\n\t\t\t//updateKeyControl();\n\t\t\t\n\t\t\tif ( ! quickLook )\n\t\t\t{\n\t\t\t\t//checkHighlight();\n\t\t\t}\n\t\t\t\n\t\t\t\n\t\t\tif ( fpsDisplay )\n\t\t\t{\n\t\t\t\tfpsDisplay.innerHTML = \'fps: \' + Math.round(tweenFrames * 1000 / tweenLength);\n\t\t\t}\n\t\t}\n\t\t\n\t\tdraw();\n\t}\n\t\n\tprogressLast = progress;\n}\nfunction updateDatasetButtons()\n{\n\tif ( datasets == 1 )\n\t{\n\t\treturn;\n\t}\n\t\n\tvar node = selectedNode ? selectedNode : head;\n\t\n\tdatasetButtonLast.disabled =\n\t\tnode.attributes[magnitudeIndex][lastDataset] == 0;\n\t\n\tdatasetButtonPrev.disabled = true;\n\tdatasetButtonNext.disabled = true;\n\t\n\tfor ( var i = 0; i < datasets; i++ )\n\t{\n\t\tvar disable = node.attributes[magnitudeIndex][i] == 0;\n\t\t\n\t\tdatasetDropDown.options[i].disabled = disable;\n\t\t\n\t\tif ( ! disable )\n\t\t{\n\t\t\tif ( i != currentDataset )\n\t\t\t{\n\t\t\t\tdatasetButtonPrev.disabled = false;\n\t\t\t\tdatasetButtonNext.disabled = false;\n\t\t\t}\n\t\t}\n\t}\n}\nfunction updateDatasetWidths()\n{\n\tif ( datasets > 1 )\n\t{\n\t\tfor ( var i = 0; i < datasets; i++ )\n\t\t{\n\t\t\tcontext.font = fontBold;\n\t\t\tvar dim = context.measureText(datasetNames[i]);\n\t\t\tdatasetWidths[i] = dim.width;\n\t\t}\n\t}\n}\nfunction updateKeyControl()\n{\n\tif ( keys == 0 )//|| progress != 1 )\n\t{\n\t\tkeyControl.style.visibility = \'hidden\';\n\t}\n\telse\n\t{\n\t\tkeyControl.style.visibility = \'visible\';\n\t\tkeyControl.style.right = margin + \'px\';\n\t\t\n\t\tif ( showKeys )\n\t\t{\n\t\t\tkeyControl.style.top =\n\t\t\t\timageHeight -\n\t\t\t\t(\n\t\t\t\t\tkeys * (keySize + keyBuffer) -\n\t\t\t\t\tkeyBuffer +\n\t\t\t\t\tmargin +\n\t\t\t\t\tkeyControl.clientHeight * 1.5\n\t\t\t\t) + \'px\';\n\t\t}\n\t\telse\n\t\t{\n\t\t\tkeyControl.style.top =\n\t\t\t\t(imageHeight - margin - keyControl.clientHeight) + \'px\';\n\t\t}\n\t}\n}\nfunction updateView()\n{\n\tif ( selectedNode.depth > maxAbsoluteDepth - 1 )\n\t{\n\t\tmaxAbsoluteDepth = selectedNode.depth + 1;\n\t}\n\t\n\thighlightedNode = selectedNode;\n\t\n\tangleFactor = 2 * Math.PI / (selectedNode.magnitude);\n\t\n\tmaxPossibleDepth = Math.floor(gRadius / (fontSize * minRingWidthFactor));\n\t\n\tif ( maxPossibleDepth < 4 )\n\t{\n\t\tmaxPossibleDepth = 4;\n\t}\n\t\n\tvar minRadiusInner = fontSize * 8 / gRadius;\n\tvar minRadiusFirst = fontSize * 6 / gRadius;\n\tvar minRadiusOuter = fontSize * 5 / gRadius;\n\t\n\tif ( .25 < minRadiusInner )\n\t{\n\t\tminRadiusInner = .25;\n\t}\n\t\n\tif ( .15 < minRadiusFirst )\n\t{\n\t\tminRadiusFirst = .15;\n\t}\n\t\n\tif ( .15 < minRadiusOuter )\n\t{\n\t\tminRadiusOuter = .15;\n\t}\n\t\n\t// visibility of nodes depends on the depth they are displayed at,\n\t// so we need to set the max depth assuming they can all be displayed\n\t// and iterate it down based on the deepest child node we can display\n\t//\n\tvar maxDepth;\n\tvar newMaxDepth = selectedNode.getMaxDepth() - selectedNode.getDepth() + 1;\n\t//\n\tdo\n\t{\n\t\tmaxDepth = newMaxDepth;\n\t\t\n\t\tif ( ! compress && maxDepth > maxPossibleDepth )\n\t\t{\n\t\t\tmaxDepth = maxPossibleDepth;\n\t\t}\n\t\t\n\t\tif ( compress )\n\t\t{\n\t\t\tcompressedRadii = new Array(maxDepth);\n\t\t\t\n\t\t\tcompressedRadii[0] = minRadiusInner;\n\t\t\t\n\t\t\tvar offset = 0;\n\t\t\t\n\t\t\twhile\n\t\t\t(\n\t\t\t\tlerp\n\t\t\t\t(\n\t\t\t\t\tMath.atan(offset + 2),\n\t\t\t\t\tMath.atan(offset + 1),\n\t\t\t\t\tMath.atan(maxDepth + offset - 1),\n\t\t\t\t\tminRadiusInner,\n\t\t\t\t\t1 - minRadiusOuter\n\t\t\t\t) - minRadiusInner > minRadiusFirst &&\n\t\t\t\toffset < 10\n\t\t\t)\n\t\t\t{\n\t\t\t\toffset++;\n\t\t\t}\n\t\t\t\n\t\t\toffset--;\n\t\t\t\n\t\t\tfor ( var i = 1; i < maxDepth; i++ )\n\t\t\t{\n\t\t\t\tcompressedRadii[i] = lerp\n\t\t\t\t(\n\t\t\t\t\tMath.atan(i + offset),\n\t\t\t\t\tMath.atan(offset),\n\t\t\t\t\tMath.atan(maxDepth + offset - 1),\n\t\t\t\t\tminRadiusInner,\n\t\t\t\t\t1 - minRadiusOuter\n\t\t\t\t)\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tnodeRadius = 1 / maxDepth;\n\t\t}\n\t\t\n\t\tnewMaxDepth = selectedNode.maxVisibleDepth(maxDepth);\n\t\t\n\t\tif ( compress )\n\t\t{\n\t\t\tif ( newMaxDepth <= maxPossibleDepth )\n\t\t\t{\n//\t\t\t\tcompress\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tif ( newMaxDepth > maxPossibleDepth )\n\t\t\t{\n\t\t\t\tnewMaxDepth = maxPossibleDepth;\n\t\t\t}\n\t\t}\n\t}\n\twhile ( newMaxDepth < maxDepth );\n\t\n\tmaxDisplayDepth = maxDepth;\n\t\n\tlightnessFactor = (lightnessMax - lightnessBase) / (maxDepth > 8 ? 8 : maxDepth);\n\tkeys = 0;\n\t\n\tnLabelOffsets = new Array(maxDisplayDepth - 1);\n\tlabelOffsets = new Array(maxDisplayDepth - 1);\n\tlabelLastNodes = new Array(maxDisplayDepth - 1);\n\tlabelFirstNodes = new Array(maxDisplayDepth - 1);\n\t\n\tfor ( var i = 0; i < maxDisplayDepth - 1; i++ )\n\t{\n\t\tif ( compress )\n\t\t{\n\t\t\tif ( i == maxDisplayDepth - 1 )\n\t\t\t{\n\t\t\t\tnLabelOffsets[i] = 0;\n\t\t\t}\n\t\t\telse\n\t\t\t{\n\t\t\t\tvar width =\n\t\t\t\t\t(compressedRadii[i + 1] - compressedRadii[i]) *\n\t\t\t\t\tgRadius;\n\t\t\t\t\n\t\t\t\tnLabelOffsets[i] = Math.floor(width / fontSize / 1.2);\n\t\t\t\t\n\t\t\t\tif ( nLabelOffsets[i] > 2 )\n\t\t\t\t{\n\t\t\t\t\tnLabelOffsets[i] = min\n\t\t\t\t\t(\n\t\t\t\t\t\tMath.floor(width / fontSize / 1.75),\n\t\t\t\t\t\t5\n\t\t\t\t\t);\n\t\t\t\t}\n\t\t\t}\n\t\t}\n\t\telse\n\t\t{\n\t\t\tnLabelOffsets[i] = Math.max\n\t\t\t(\n\t\t\t\tMath.floor(Math.sqrt((nodeRadius * gRadius / fontSize)) * 1.5),\n\t\t\t\t3\n\t\t\t);\n\t\t}\n\t\t\n\t\tlabelOffsets[i] = Math.floor((nLabelOffsets[i] - 1) / 2);\n\t\tlabelLastNodes[i] = new Array(nLabelOffsets[i] + 1);\n\t\tlabelFirstNodes[i] = new Array(nLabelOffsets[i] + 1);\n\t\t\n\t\tfor ( var j = 0; j <= nLabelOffsets[i]; j++ )\n\t\t{\n\t\t\t// these arrays will allow nodes with neighboring labels to link to\n\t\t\t// each other to determine max label length\n\t\t\t\n\t\t\tlabelLastNodes[i][j] = 0;\n\t\t\tlabelFirstNodes[i][j] = 0;\n\t\t}\n\t}\n\t\n\tfontSizeText.innerHTML = fontSize;\n\tfontNormal = fontSize + \'px \' + fontFamily;\n\tcontext.font = fontNormal;\n\tfontBold = \'bold \' + fontSize + \'px \' + fontFamily;\n\ttickLength = fontSize * .7;\n\t\n\thead.setTargets(0);\n\t\n\tkeySize = ((imageHeight - margin * 3) * 1 / 2) / keys * 3 / 4;\n\t\n\tif ( keySize > fontSize * maxKeySizeFactor )\n\t{\n\t\tkeySize = fontSize * maxKeySizeFactor;\n\t}\n\t\n\tkeyBuffer = keySize / 3;\n\t\n\tfontSizeLast = fontSize;\n\t\n\tif ( datasetChanged )\n\t{\n\t\tdatasetChanged = false;\n\t}\n\telse\n\t{\n\t\tdatasetAlpha.start = 0;\n\t}\n\t\n\tvar date = new Date();\n\ttweenStartTime = date.getTime();\n\tprogress = 0;\n\ttweenFrames = 0;\n\t\n\tupdateKeyControl();\n\tupdateDatasetWidths();\n\t\n\tdocument.title = \'Krona - \' + selectedNode.name;\n\tupdateNavigationButtons();\n\tsnapshotButton.disabled = true;\n\t\n\tmaxAbsoluteDepthText.innerHTML = maxAbsoluteDepth - 1;\n\t\n\tmaxAbsoluteDepthButtonDecrease.disabled = (maxAbsoluteDepth == 2);\n\tmaxAbsoluteDepthButtonIncrease.disabled = (maxAbsoluteDepth == head.maxDepth);\n\t\n\tif ( collapse != collapseLast && search.value != \'\' )\n\t{\n\t\tonSearchChange();\n\t\tcollapseLast = collapse;\n\t}\n}\nfunction updateMaxAbsoluteDepth()\n{\n\twhile ( maxAbsoluteDepth > 1 && selectedNode.depth > maxAbsoluteDepth - 1 )\n\t{\n\t\tselectedNode = selectedNode.getParent();\n\t}\n}\nfunction updateNavigationButtons()\n{\n\tbackButton.disabled = (nodeHistoryPosition == 0);\n//\tupButton.disabled = (selectedNode.getParent() == 0);\n\tforwardButton.disabled = (nodeHistoryPosition == nodeHistory.length);\n}\nfunction useHue()\n{\n\treturn useHueCheckBox && useHueCheckBox.checked;\n}\n/*\nfunction zoomOut()\n{\n\treturn (\n\t\tselectedNodeLast != 0 &&\n\t\tselectedNodeLast.getDepth() < selectedNode.getDepth());\n}\n*/'
_KRONA_IMAGE_SHORTCUT = 'data:image/x-icon;base64,AAABAAMAEBAAAAEAIABoBAAANgAAABgYAAABACAAiAkAAJ4EAAAgIAAAAQAgAKgQAAAmDgAAKAAAABAAAAAgAAAAAQAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP///wCAgIAC////AP///wC0tKJHlZWSqI6OmuRxcXn9koaK9J2Ym8uNm5V73PPoFv///wAAAAAD////AP///wCAgIAC////AP//7xCHh4OodHSi/2xqw/9nZ9H/XVuU/9eFqv/njb7/yYiq/5GAiOifqKJV////AAAAAAL///8A////AO3t2w56eoHPdHLH/3Vz3vx4dtz7dXTb/mVimv/XhKn/8ZPG/fSVx/vvlcT/kG+E/32Himb///8AAAAAA////wB0dHageHbO/3585fl5d9r/dnTX/3V02/5mYpz80X2i/eaIuv/ukcD/wXak/I6Maf+OkXD9iIKjL////wCEhFo+c3Kv/4F+6P18etz/e3nd/3l33P57eeT/aWag/9F5oP/sib3+tmuZ/3uEXf/P1YX7zdKG/3x9dbX///8Ac3J9nYKA4f+CgOb7fnzf/3174P58et//bm2v8Vtgc7SdbIHNqWmM/4aLZP/L0YT/0dWH/tfdiP+go3P6YlyWJ3BvmduGg+n/g4Hk/oF+4/6AfuX/bm2k7pCUUUX///8Ajv/jCVBmXqPBxoH/1NmI/c7Uhf/b4Yn9wMWB/2BdcGBraqX6iYbs/4SB5P+Eguj7hYPm/2Vlcov///8AACRtB////wAAAIAWrLB29uHojv/X3on+3eSK+9PYh/9XWVN7bWum+ouI7/+Gg+f/hoTq+4iF6f9lZXKL////ACQAbQf///8AAABoFoF+Xva0snL/y818/uDmiPvZ34n/WVlTe3Jxm9uMiu//iofr/oiG6f6Jhu3/cnGo7oyQQ0X///8A/znjCXF8dKNotpL/Za2O/WuVe/+FkGv9mJ5w/2pqaGBzc3+djYrr/42K8fuKh+v/i4js/oqH6/93dbfxYmJztGaTec1yxJr/gOGu/4nptf+K6Lf+i+q8/2eXhfqDNG8ngIBKPnt6uP+Rjvf9jInr/4yJ7f+Miu3+jor0/25xoP920pz/gOGw/oHerv+H4rL/kPC9+5Ttvf98ioC1////AP///wB0dHSgioff/5SR+/mOi+7/jYvt/4+L8P5scJ78etKe/YXis/+H47P/kO+8/JLxv/92m4f9nXKNL////wD///8A7e3bDnt7hM+Ihdv/ko/4/JSS+PuUkPb+bnKg/4Lapv+Q7739kfO/+43puf9/qJH/mYeUZv///wAAAAADgICAAv///wD//+8QhoaBqIB/sf+Iht7/jYjw/2tvnf+D3Kf/ieS1/4bIpP+AkYfoqJOfVf///wAAAAAC////AP///wCAgIAC////AP///wCwsJdHlZWQqJORn+RzdHv9hZOK9Jecmcuej5l7/+jzFv///wAAAAAD////AP///wD4PwAA4A8AAMAHAACAAwAAgAEAAAABAAADgQAAA8EAAAPBAAADgQAAAAEAAIABAACAAwAAwAcAAOAPAAD4PwAAKAAAABgAAAAwAAAAAQAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP///wD///8A////AP///wEAAAAB////AP///wi6urdDubmxj6ennsuQkIzvc3Nz/YSIh/WWoZzWrrezoLe6t1X///8T////AP///wCAgIAC////AP///wD///8A////AP///wD///8A////Af///wD///8Ax8fDQJeYkLV4eH39WlmM/19dpP9dX6b/WFRr/611j/+/fJ//om2I/4N3ff+Jko7Ntri2Xv///wD///8AgICAAv///wD///8A////AP///wD///8B////AP///wOfn5t9c3N29GRjn/9oZ8P/c3LU/XVy3P9rbdH/YVt8/9eIrv/5mcz/75bE/t+Mt/67epz/fnN5/4uTkKP///8W////AAAAAAL///8A////AP///wH///8A////A5iYkZJoaH3/bWvB/3h23Px4dt3+dHPX/3Rx1v9sbcv/YFp6/82ApP/rj8D/65DA//KVxv/1mcn83424/5t1iP98gX7A////F////wCAgIAC////AAAAAAL///8AkJCHdWZlf/9zcs7/fHrj/Hh32/92ddn/dXTX/3Zz2P9vcM//YFt7/8t9of/pi73/54u7/+eNvf/vk8P/7ZbE/pxsh/xLVU7/iomJqP///wD///8A////Af///wCfn5Y4ZmV07XZ0zv9+fOX8enjc/3p43P94dtv/dnXZ/3d12f9xcdD/YFp6/8h5nf/mh7r/44e3/+qMvf/ljrv/jV59/3F+V/68wH3/gYNp/4uIlWf///8AAAAAAv///wJ3d22pc3G6/4KA5/t9e9//fHre/3t53f96eNz/eHbZ/3l23P50dNb9Ylx+/ch2m/3jg7b+44W2/+CItv+NXXz/bnpW/8XIgv/Y3Yv8ur57/3Fxbtni4uIa////AHl5XTdoZ4b4f3zb/4F+5P5+fN//fXvf/3x63v96edz/ennf/nt54P9xccb/Xlxy/7Rvj//gf7P/3IKw/olZeP9te1b/xsqD/9DVh//P1IX/09iG/pWYcf9nZ3Nt////AISEb35ycLb/hYPr/YB+4f9/feH/fnzh/3173/99e+L+enjX/2dmmfhlZXiwYGBiknJkaqGVZ37pg11x/3eBXP7GyoP/0NWH/87Uhf/R1ob/2eCK+8HFff9tbXWv////BG9vdb16eMn/hYPq/4F/4/+BfuP/gH7i/4B94v59e93/Y2KO+ISEbX7y8uQT////AP///weAjo5aT1tP5L7Cf//T2Ij+ztOF/9HWhv/S14b/2d+K/s3Sgv+BgnLcg4OSI11deeaBf9r/hYPo/4OB5P+DgeX/gX/i/4OB6v10c77/cXFjmv///wL///8A////AP///wD///8AeHN9apuecP/U2Yf/0teG/tPZhv/V24f/2uCJ/tTahv+KjGj6AABINV9egfuGg+P/hoTp/4WC5v+EgeX/goDk/4WD6/1wbq//WFgAQ////wAzMzMFAAAAAVVVVQP///8AAABjEoiLZ/Df5Y7/2d+J/dfch//W3If/2d+H/9vhif6ZnW3/AAA5P19egvuHheX/iIXq/4eE6P+Gg+f/hILl/4eE7P1xb7D/WFgAQ////wAzMzMFAAAAAVVVVQP///8AAABVEnd1WvDEx37/y9GD/dnfif/g54z/4OeL/97kiv6anWz/AAA1P15eeuaFg93/ioft/4iF6f+Hhen/hoPn/4mG7/14d8L/cXFimv///wP///8A////AP///wD///8AenV3a1dyZf9QY1//bGVW/pORY/+trnD/xch7/tTZhP+Pkmn6AABINXBwdb1/fc7/jYry/4mH6v+Jhur/iIbp/4iG6/6Gg+X/Z2aT+IKEaX7y8uQT////AP///wiXgY9bXndp5HXLnv95z6P+dsif/2Wrjv9RcWv/Z3Rk/nyBYv9xcWbbfHyDI4SEa354d73/kI31/YuI6/+LiOz/iofr/4mG6v+Lh+3+hoPi/21soPhnZ3qwYGBikmNza6JhinXpcsWa/3vdq/6C4a//iOi1/4zpt/+O6Lf/kOu8+3O1l/9jVmKv////BHl5WDdraor4iofn/4+M8v6Miez/i4nt/4uI7P+LiOv/i4jt/ouJ8P9+e83/XGBt/223jv952aj/e9uq/n3aqv+D36//huKy/4rntv+Q7r3/lvTC/nupj/9zYmxt////AP///wJ3d2qpfXvF/5KP9vyNiu7/jYru/4yJ7f+Miu7/i4jr/46L8P6FgN79XmZ5/XTMm/1+3q7+f9ys/4Lgr/+G5LP/iue2/4zot/+V9cL8iNWs/3BxcNni4uIa////AP///wCfn5I4aGd17YeE3v+Ukfj8jovu/46L7/+Oi+//jYru/46M8f+FgNv/XmV3/3jMnf+D4rH/g+Cw/4fks/+K57b/jem4/5TywP6Q57n/bYh5/5KIjWf///8AAAAAAgAAAAL///8AkJCFdWtqhf+IheL/lZL5/JCN8P+PjO//j4zv/5GO8/+Hgt3/X2Z4/33Rof+I5rb/h+Sz/4nmtf+N6rn/lPPB/pHruvxxm4T/h3yDqP///wD///8A////Af///wH///8A////A5aWj5JtbIP/g4DW/5SR9vyVkvj+kY7y/5KP9P+Ig97/X2d4/4DUpf+N6rn/jeq4/5Lwvv+U88D8h9qu/3KVgv+Lgoe/////F////wCAgIAC////AP///wD///8B////AP///wOfn5l9dHR39HVzsf+Fg97/k4/x/ZiW/v+Nh+f/YWl7/4ngrv+W+MX/kuy8/onesP52uZb/dH54/5aOk6P///8W////AAAAAAL///8A////AP///wD///8A////Af///wD///8Ax8fDQJWVjbV5eX79amic/3h2vv93crf/WV9r/3W1kP95vZn/aZ2C/3eCfP+Ui5DNu7i7Xv///wD///8AgICAAv///wD///8A////AP///wD///8A////AP///wEAAAAB////AP///wi6urdDuLitj6WmmcuPj4rvc3Nz/YiEh/Whlp3Wt66zoLq3ulX///8T////AP///wD///8C////AP///wD///8A////AP8B/wD8AH8A+AAfAOAADwDgAAcAwAAHAIAAAwCAAAMAgAABAAB8AQAAfgEAAP4BAAD+AQAAfgEAAHwBAIAAAQCAAAMAgAADAMAABwDgAAcA4AAPAPgAHwD8AH8A/wH/ACgAAAAgAAAAQAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD///8A////AP///wD///8A////AP///wD///8A////AP///wD///8A////KO3t7XXFxcWwqqqq25eXl/RwcHD9mJiY9Kurq9vHx8ew8PDwdf///yj///8A////AP///wD///8A////AP///wD///8A////AP///wD///8AAAAAAP///wD///8A////AP///wD///8A////AP///wD///8A9/f3TLS0tLpvb3H/UVFR/1VVWf9bW27/WVlz/1FRUf96ZG//c2Nr/1lWV/9RUVH/d3V2/7m5ubr5+flM////AP///wD///8A////AP///wD///8A////AP///wAAAAAA////AP///wD///8A////AP///wD///8A/f39HLW1tatbW13/V1de/2FgkP9rabv/cG7S/29u0/9paL7/UVFR/9GIrv/wlcX/75bE/8+Irf+ZcYb/XVhb/2FfYP+8vLyr/v7+HP///wD///8A////AP///wD///8A////AAAAAAD///8A////AP///wD///8A////AOvr60l2dnnpVVVZ/2ZlnP9ycdL/c3LW/3Jx1f9xcNT/cG/T/2ppvv9RUVH/zoSs/+yRwf/uk8L/7pPC/++UxP/pk8D/pXaP/1hVV/9/fn/p7+/vSf///wD///8A////AP///wD///8AAAAAAP///wD///8A////AP///wDg4OBWX19h/11dcv9zccr/d3XZ/3Z02P91c9f/dHLX/3Nx1v9xcNT/a2q//1FRUf/Ngqn/6o6+/+uQwP/tksH/7pPD/++Vxf/xl8b/3I64/3Jiav9lY2T/2dnZVv///wD///8A////AP///wAAAAAA////AP///wD///8A5OTkQ1tbXf9hYIH/eHfY/3l32/94dtr/d3ba/3Z12f90c9f/c3LW/3Nx1f9tbMH/UVFR/8l/pf/mirr/6Iy8/+mOvv/qj7//7JHB/++VxP/xl8b/04qx/1RUUv9ZWVj/6+vrQ////wD///8A////AAAAAAD///8A////APX19Q9mZmjfYF98/3t53P97ed3/enjc/3l32/93dtr/dnXZ/3Z02P91c9f/dXTY/25tw/9RUVH/yHyj/+WIuf/libn/54q7/+mNvf/rj7//7JHB/9WJsf9XVlX/Z2hb/3BxX/9vb23f+fn5D////wD///8AAAAAAP///wD///8AkZCSnFpaZ/97edn/fHre/3t53v97ed3/enjc/3l33P94dtv/dnXZ/3V02P92dNj/cG7D/1FRUf/FeaH/4oW2/+OGt//liLj/54q6/+mOvv/Sha3/V1VV/2doW//FyYL/xcqC/19fWf+fn5+c////AP///wAAAAAA////ANjY2DhRUVL/dXPA/3994f9+fOD/fXvf/3x63v97ed3/enjc/3l32/94dtr/eHba/3d12f9xb8T/UVFR/8J1nf/fgLH/4IOz/+KEtf/liLn/zoCq/1dVVf9naFr/w8eA/83Shf/P1Yb/q652/1JSUf/h4eE4////AAAAAAD///8Afn6ApWNigv+BfuL/gH3h/3584P99e9//fHre/3x63v97ed3/eXjc/3l32/95d9v/d3ba/3Fwxf9RUVH/wXSc/9x9r//fgbL/4YO0/8t9pv9XVVT/Z2ha/8LHgf/N0oX/ztOF/8/Uhv/S14b/c3Vh/42Ni6X///8AAAAAAOfn5xJRUVH/eXbG/4F/4/+AfuL/f33h/3584P99e+D/fHrf/3x63v97ed3/e3ne/3p43f9ycMH/Y2OR/1FRUf+NY3j/w3Sd/9x9rv/IeaL/V1VU/2doWv/DyIH/zNGF/87Thf/Q1Yb/0daG/9TZh/+ytnj/UlJS/+/v7xIAAAAAq6urYF1cav+DgeX/goDk/4F/4/+AfuL/f33h/3994f9+fOD/fXvf/3x63v95d9X/YF99/1FRUb9sbG+AZmVmgG1qa4BRUVG/fWBv/1dVVP9naFr/xMiB/83Shv/O04b/z9WG/9HWhv/S14b/1NqH/9TZhv9fYFn/t7e3YAAAAAB6enyeaGeV/4SC5v+DgeX/goDk/4F+4/+BfuP/f33h/3584P9+fOD/fHrb/1xcbv9sbG6g39/fQP///wD///8A////ANra2kBbWlugVVVT/8LHgf/O04X/ztOF/8/Uhf/Q1Yb/0teG/9PZh//V24f/2N6I/4aIZv+IiIeeAAAAAFpaXMx0crX/hYLm/4SC5v+DgeX/goDk/4KA5P+Bf+P/gH7i/3994f9oZ5v/X19hv/Hx8SD///8A////AP///wD///8A////AOfn5yBcXFq/m55w/9DVhv/R1ob/0teG/9LXhv/T2Yb/1dqH/9jdiP/Z34j/o6Zy/2VlY8wAAAAAUVFR7Hx7zP+GhOj/hYLm/4SB5f+Egub/g4Hl/4KA5P+Bf+P/gX/j/1paZP+xsbFA////AP///wD///8A////AP///wD///8A////AKampkBqa17/0teG/9LXhv/T2Yf/1dqH/9bciP/Y3oj/2N6H/9rgiP+5vXr/UlJS7AAAAABRUVH7gX7W/4eE6P+GhOj/hYPn/4WC5v+EgeX/g4Hl/4OB5f9/fNv/UVFR/9zc3AD///8A////AP///wD///8A////AP///wD///8A0NDQAFRUU//R14b/1dqH/9Xbh//W3If/192H/9jeh//a4Ij/3OKI/8PJfv9RUVH7AAAAAFFRUfuCgNj/iIbq/4eE6P+HhOj/hoPn/4WC5v+Egub/hILm/4B+3f9RUVH/3NzcAP///wD///8A////AP///wD///8A////AP///wDLy8sAUlJR/87ThP/W3Ij/192H/9nfiP/Z34j/2uCI/9vhiP/d44j/xMh+/1FRUfsAAAAAUVFR7H99zv+Jh+v/iIXq/4iF6v+Hhen/hoTo/4aD5/+Fgub/hYLm/1taZP+xsbFA////AP///wD///8A////AP///wD///8A////AJWVlUBSUlH/UlJR/25vXf+Ul2z/ur97/9jeh//c44j/3uWJ/97lif+8wXr/UVFR7AAAAABaWlzMd3W4/4qH6/+Kh+v/iYbq/4iF6f+HhOj/iIXp/4aE6P+GhOj/a2qe/19fYr/x8fEg////AP///wD///8A////AP///wDp6ekgW1xcv2SWe/9nm4D/Wm9k/1daVv9RUVH/VlZT/3p8Yf+ipnH/ys+A/6ircv9lZmPMAAAAAHp6fJ5rapf/i4js/4uJ7f+KiOz/iYbq/4mG6v+Ihur/iIbq/4eE6P+FguP/Xl5w/2xsb6Df399A////AP///wD///8A2traQGVnZqBccGX/d9Sk/3vZqf+B367/g+Gw/33No/9qk33/W2lh/1ZXVP9RUVH/U1NR/3d3dZ4AAAAAq6urYF5ea/+Ni+//jInt/4uI7P+LiOz/iofr/4mG6v+Jhur/iIXp/4iF6f+EguD/ZGSB/1FRUb9ra2+AZWdmgGpubIBRUVG/Xn1u/3TOn/9516f/ftyr/4Herf+G47L/iOa1/4vot/+Q7bz/k+++/4PHo/9SVFP/oaGhYAAAAADn5+cSUVFR/4F/z/+Niu7/jInt/4yJ7f+LiOz/i4nt/4qI7P+LiOz/iofr/4mG6v+Jhur/fnzN/2ZljP9RUVH/YpB5/2++lf921KT/edem/3vZqf+A3q3/hOGw/4bjs/+J5rX/juu6/5Dsu/+U8cD/g8Wi/1JTUv/t7e0SAAAAAP///wB+foClZ2aG/46L8P+Oi/D/jYru/4yJ7f+Mie3/i4js/4uI7P+LiOz/iofr/4qH6/+KiOz/eXjA/1FRUf9zyJz/eden/3zaqf982qn/ftys/4Hfrv+G47L/iOW0/4rntv+N6rn/ku++/5XxwP9le2//iYqKpf///wAAAAAA////ANjY2DhRUVL/gH7L/4+M8P+Oi+//jovv/42K7v+Niu7/jIru/4yK7v+Mie3/jInt/4uI7P96eMD/UVFR/3bLoP972an/ftyr/4Herv+D4K//heOy/4jmtf+L6Lf/juu5/5DtvP+V8cD/gcGf/1NUU//l5eU4////AAAAAAD///8A////AJGRkpxdXGn/jYrp/5CN8f+PjPD/j4zw/46L7/+Oi+//jYru/42K7v+Niu7/jYru/3x5wv9RUVH/ec6i/4Herv+B367/hOGw/4bjsv+I5bT/i+i3/43quf+Q7bz/k+++/4/jtv9bY17/nZ6enP///wD///8AAAAAAP///wD///8A9vb2D2Zmad9mZYL/kI3w/5CN8f+QjfH/j4zw/4+M8P+PjPD/j4zw/46L7/+Oi+//fHrD/1FRUf9/0qb/hOGw/4Xjsv+F47L/iOa1/4vot/+O67r/kO28/5Lvvf+S7b3/ZHht/3N1dN/6+voP////AP///wAAAAAA////AP///wD///8A5OTkQ1xcXv9oZ4j/kI3v/5GN8f+RjfH/kI3x/5CN8f+QjfH/kI3x/5CN8f99e8P/UVFR/3/UqP+I5bT/iOW0/4vot/+L6Lf/juu5/4/su/+S773/kem6/2Z9cf9gYmH/6+vrQ////wD///8A////AAAAAAD///8A////AP///wD///8A4ODgVl9fYf9iYnf/iojh/5KP8/+Sj/P/kY7y/5GO8v+RjvL/kY7y/358xf9RUVH/hdmt/4vot/+O67n/juu5/5Dtu/+Q7bv/ku++/4rXrv9hcGj/aGpp/+np6Vb///8A////AP///wD///8AAAAAAP///wD///8A////AP///wD///8A7OzsSXd3eelWVlr/dXOr/5CO7v+TkPT/k5D0/5KP8/+Sj/P/f33F/1FRUf+I267/kO28/5Dsu/+Q7Lv/ku+9/5Douv9zoon/U1VU/4GDgunx8fFJ////AP///wD///8A////AP///wAAAAAA////AP///wD///8A////AP///wD///8A/f39HLW1tatbW13/WVlg/3Bunv+Gg9T/k4/y/5SR9f+Bfsf/UVFR/4vfsv+S773/kOu7/4TKpf9ul4H/VlpY/19hYP/AwMCr////HP///wD///8A////AP///wD///8A////AAAAAAD///8A////AP///wD///8A////AP///wD///8A////APf390y0tLS6cG9y/1FRUf9WVlr/Y2J2/2Fgdv9RUVH/ZYBx/2N1a/9VV1b/UVFR/3R2df+/v7+6+vr6TP///wD///8A////AP///wD///8A////AP///wD///8AAAAAAP///wD///8A////AP///wD///8A////AP///wD///8A////AP///wD///8o7e3tdcXFxbCqqqrbl5eX9HBwcP2YmJj0q6ur28fHx7Dv7+91////KP///wD///8A////AP///wD///8A////AP///wD///8A////AP///wAAAAAA//Af//+AA//+AAD//AAAf/gAAD/wAAAf4AAAD8AAAAfAAAAHgAAAA4AAAAOAAAADAAfAAQAP4AEAH/ABAB/wAQAf8AEAH/ABAA/gAQAHwAGAAAADgAAAA4AAAAPAAAAHwAAAB+AAAA/wAAAf+AAAP/wAAH/+AAD//4AD///wH/8='
_KRONA_IMAGE_HIDDEN = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAAICAYAAADED76LAAAAAXNSR0IArs4c6QAAAAZiS0dEAP8A/wD/oL2nkwAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB9oLCBQhNQwWVnsAAAAidEVYdENvbW1lbnQAQ3JlYXRlZCB3aXRoIEdJTVAgb24gYSBNYWOHqHdDAAABE0lEQVQYGQEIAff+AwAAABkAAAAAAAAA+gAAAAAAAAAAAAAAAAAAAAAAAAAMAwAAAAAAAAANAAAAAAAAAPoAAAAAAAAADAAAAAYAAAD0AwAAAPoAAAAAAAAAAAAAAPoAAAAMAAAADQAAAPoAAAD6AAAAAAAAAAAAAAAAAAAAAAwAAAAZAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAABkAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAGQAAAAwAAAAAAAAADAAAAAwAAAAABAAAAAAAAAAAAAAA8wAAAPQAAAAAAAAAAAAAAA0AAAAAAAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAZRssKC5OpXwYAAAAASUVORK5CYII='
_KRONA_IMAGE_LOADING = 'data:image/gif;base64,R0lGODlhEAAQAPIAAP///wAAAMLCwkJCQgAAAGJiYoKCgpKSkiH/C05FVFNDQVBFMi4wAwEAAAAh/hpDcmVhdGVkIHdpdGggYWpheGxvYWQuaW5mbwAh+QQJCgAAACwAAAAAEAAQAAADMwi63P4wyklrE2MIOggZnAdOmGYJRbExwroUmcG2LmDEwnHQLVsYOd2mBzkYDAdKa+dIAAAh+QQJCgAAACwAAAAAEAAQAAADNAi63P5OjCEgG4QMu7DmikRxQlFUYDEZIGBMRVsaqHwctXXf7WEYB4Ag1xjihkMZsiUkKhIAIfkECQoAAAAsAAAAABAAEAAAAzYIujIjK8pByJDMlFYvBoVjHA70GU7xSUJhmKtwHPAKzLO9HMaoKwJZ7Rf8AYPDDzKpZBqfvwQAIfkECQoAAAAsAAAAABAAEAAAAzMIumIlK8oyhpHsnFZfhYumCYUhDAQxRIdhHBGqRoKw0R8DYlJd8z0fMDgsGo/IpHI5TAAAIfkECQoAAAAsAAAAABAAEAAAAzIIunInK0rnZBTwGPNMgQwmdsNgXGJUlIWEuR5oWUIpz8pAEAMe6TwfwyYsGo/IpFKSAAAh+QQJCgAAACwAAAAAEAAQAAADMwi6IMKQORfjdOe82p4wGccc4CEuQradylesojEMBgsUc2G7sDX3lQGBMLAJibufbSlKAAAh+QQJCgAAACwAAAAAEAAQAAADMgi63P7wCRHZnFVdmgHu2nFwlWCI3WGc3TSWhUFGxTAUkGCbtgENBMJAEJsxgMLWzpEAACH5BAkKAAAALAAAAAAQABAAAAMyCLrc/jDKSatlQtScKdceCAjDII7HcQ4EMTCpyrCuUBjCYRgHVtqlAiB1YhiCnlsRkAAAOwAAAAAAAAAAAA=='
_KRONA_IMAGE_LOGO = 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAANgAAAA8CAYAAAAE9XR5AAAC0WlDQ1BJQ0MgUHJvZmlsZQAAOI2NlM9LFGEYx7+zjRgoQWBme4ihQ0ioTBZlROWuv9i0bVl/lBLE7Oy7u5Ozs9PM7JoiEV46ZtE9Kg8e+gM8eOiUl8LALALpblFEgpeS7Xlnxt0R7ccLM/N5nx/f53nf4X2BGlkxTT0kAXnDsZJ9Uen66JhU+xEhHEEdwqhTVNuMJBIDoMFjsWtsvofAvyute/v/OurStpoHhP1A6Eea2Sqw7xfZC1lqBBC5XsOEYzrE9zhbnv0x55TH8659KNlFvEh8QDUtHv+auEPNKWmgRiRuyQZiUgHO60XV7+cgPfXMGB6k73Hq6S6ze3wWZtJKdz9xG/HnNOvu4ZrE8xmtN0bcTM9axuod9lg4oTmxIY9DI4YeH/C5yUjFr/qaoulEk9v6dmmwZ9t+S7mcIA4TJ8cL/TymkXI7p3JD1zwW9KlcV9znd1Yxyeseo5g5U3f/F/UWeoVR6GDQYNDbgIQk+hBFK0xYKCBDHo0iNLIyN8YitjG+Z6SORIAl8q9TzrqbcxtFyuZZI4jGMdNSUZDkD/JXeVV+Ks/JX2bDxeaqZ8a6qanLD76TLq+8ret7/Z48fZXqRsirI0vWfGVNdqDTQHcZYzZcVeI12P34ZmCVLFCpFSlXadytVHJ9Nr0jgWp/2j2KXZpebKrWWhUXbqzUL03v2KvCrlWxyqp2zqtxwXwmHhVPijGxQzwHSbwkdooXxW6anRcHKhnDpKJhwlWyoVCWgUnymjv+mRcL76y5o6GPGczSVImf/4RVyGg6CxzRf7j/c/B7xaOxIvDCBg6frto2ku4dIjQuV23OFeDCN7oP3lZtzXQeDj0BFs6oRavkSwvCG4pmdxw+6SqYk5aWzTlSuyyflSJ0JTEpZqhtLZKi65LrsiWL2cwqsXQb7Mypdk+lnnal5lO5vEHnr/YRsPWwXP75rFzeek49rAEv9d/AvP1FSOihagAAAAZiS0dEAP8A/wD/oL2nkwAAAAlwSFlzAAALEwAACxMBAJqcGAAAAAd0SU1FB+AEEhYRJKJxkPUAACAASURBVHgBAGaAmX8B/////wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAvUj3EAAAIABJREFUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA////AAAAAAD+/v4AAAAAAP///wD///8AAAAAAP///wAAAAAAAQEBAAEBAQAAAAAAAQEBAAICAgAAAAAAAQEBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP///wAAAAAA/v7+AP7+/gD9/f0A/v7+AP39/QD+/v4A/f39AP39/QD9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'

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
    Create a compact Krona working table containing taxonomy + total abundance.
    Never copy the full feature x sample abundance matrix.
    """
    tax_cols = [c for c in TAX_RANKS if c in final_df.columns]
    df = final_df.loc[:, tax_cols].copy()

    if "Total_Abundance" in final_df.columns:
        abundance = pd.to_numeric(
            final_df["Total_Abundance"],
            errors="coerce"
        ).fillna(0)
    else:
        totals = np.zeros(len(final_df), dtype=RAM_SAFE_SUM_DTYPE)
        valid_sites = [c for c in site_cols if c in final_df.columns]
        for start in range(0, len(valid_sites), RAM_SAFE_SITE_CHUNK):
            cols = valid_sites[start:start + RAM_SAFE_SITE_CHUNK]
            block = final_df.loc[:, cols].to_numpy(dtype=RAM_SAFE_SUM_DTYPE, copy=True)
            totals += block.sum(axis=1, dtype=RAM_SAFE_SUM_DTYPE)
            del block
        abundance = pd.Series(totals, index=final_df.index)
        del totals

    for rank in TAX_RANKS:
        if rank not in df.columns:
            df[rank] = ""
        df[rank] = df[rank].apply(lambda x: _clean_krona_taxon(x, rank))

    df["_krona_abundance"] = abundance.to_numpy(copy=False)
    if KRONA_MIN_ABUNDANCE is not None:
        df = df.loc[df["_krona_abundance"] >= KRONA_MIN_ABUNDANCE].copy()

    gc.collect()
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
        ktImportText Krona_Total_Abundance_with_header_{collapse_label}.tsv -o Krona_Total_Abundance_{collapse_label}.html
    """
    krona_dir.mkdir(parents=True, exist_ok=True)

    # Include collapse strategy in Krona output filenames to avoid overwriting
    # files generated with species_only, genus, or all strategies.
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()

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

    # V1.7.7 change:
    # Export only one Krona total table, with headers, to avoid redundant files.
    # The collapse strategy suffix prevents overwriting between species_only, genus, and all.
    krona_path = krona_dir / f"Krona_Total_Abundance_with_header_{collapse_label}.tsv"
    krona_table[["abundance"] + used_ranks].to_csv(
        krona_path,
        sep="\t",
        index=False,
        header=True,
        encoding="utf-8"
    )

    print(f"   ✅ Krona total table with header: {krona_path}")
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

    # Include collapse strategy in per-sample Krona filenames.
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()

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
        out_path = sample_dir / f"{safe_sample}_Krona_{collapse_label}.tsv"
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


def _krona_dataset_name() -> str:
    """Return the single dataset label used by the total Krona visualization."""
    return f"{MODE}_{SUBSET_MODE}_{str(COLLAPSE_STRATEGY).strip().lower()}"


def _build_krona_tree(krona_table: pd.DataFrame, used_ranks: list) -> dict:
    """Build a nested Krona tree from the already-aggregated MetaDiv table."""
    root = {"name": _krona_dataset_name(), "magnitude": 0, "children": {}}
    columns = list(krona_table.columns)
    for row in krona_table.itertuples(index=False, name=None):
        row_map = dict(zip(columns, row))
        try:
            amount = int(row_map.get("abundance", 0) or 0)
        except Exception:
            amount = 0
        if amount <= 0:
            continue
        root["magnitude"] += amount
        current = root
        for rank in used_ranks:
            taxon = str(row_map.get(rank, "") or "").strip()
            if not taxon:
                taxon = f"Unclassified_{rank}"
            child = current["children"].get(taxon)
            if child is None:
                child = {"name": taxon, "magnitude": 0, "children": {}}
                current["children"][taxon] = child
            child["magnitude"] += amount
            current = child
    return root


def _krona_node_xml(node: dict, indent: int = 12) -> str:
    """Serialize one Krona node recursively using the legacy Excel schema."""
    pad = " " * indent
    name = html_lib.escape(str(node.get("name", "")), quote=True)
    magnitude = int(node.get("magnitude", 0) or 0)
    lines = [
        f'{pad}<node name="{name}">',
        f'{pad}  <magnitude><val>{magnitude}</val></magnitude>',
        f'{pad}  <score><val></val></score>',
    ]
    children = node.get("children", {})
    ordered_children = sorted(
        children.values(),
        key=lambda x: (-int(x.get("magnitude", 0) or 0), str(x.get("name", "")))
    )
    for child in ordered_children:
        lines.append(_krona_node_xml(child, indent + 2))
    lines.append(f'{pad}</node>')
    return "\n".join(lines)


def build_krona_html(krona_table: pd.DataFrame, krona_dir: Path) -> Path:
    """
    Render a self-contained Krona HTML from MetaDiv's existing Krona table.

    Mapping:
      Dataset = one MetaDiv subset label for all rows
      Amount  = abundance
      Score   = blank
      Category 1..N = domain..selected terminal taxonomic rank

    No abundance recalculation, Excel, KronaTools, network access, or extra
    Python package is required.
    """
    used_ranks = krona_table.attrs.get("used_ranks", get_krona_ranks())
    if krona_table.empty:
        raise ValueError("Cannot build Krona HTML from an empty Krona table.")
    tree = _build_krona_tree(krona_table, used_ranks)
    dataset_name = html_lib.escape(_krona_dataset_name(), quote=False)
    node_xml = _krona_node_xml(tree, indent=12)
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()
    html_path = krona_dir / f"Krona_Total_Abundance_{collapse_label}.html"
    document = (
        '<!DOCTYPE html>\n'
        '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
        '<head>\n<meta charset="UTF-8"/>\n'
        f'<title>MetaDiv Builder Krona - {dataset_name}</title>\n'
        f'<link rel="shortcut icon" href="{_KRONA_IMAGE_SHORTCUT}"/>\n'
        '<script>\n' + _KRONA_LEGACY_JS + '\n</script>\n</head>\n<body>\n'
        f'<img id="hiddenImage" src="{_KRONA_IMAGE_HIDDEN}" style="display:none"/>\n'
        f'<img id="loadingImage" src="{_KRONA_IMAGE_LOADING}" style="display:none"/>\n'
        f'<img id="logo" src="{_KRONA_IMAGE_LOGO}" style="display:none"/>\n'
        '<noscript>Javascript must be enabled to view this page.</noscript>\n'
        '<div style="display:none">\n  <krona>\n'
        '    <attributes magnitude="magnitude">\n'
        '      <attribute display="Amount">magnitude</attribute>\n'
        '      <attribute display="Score">score</attribute>\n'
        '    </attributes>\n'
        '    <datasets>\n'
        f'      <dataset>{dataset_name}</dataset>\n'
        '    </datasets>\n'
        + node_xml + '\n'
        '  </krona>\n</div>\n</body>\n</html>\n'
    )
    html_path.write_text(document, encoding="utf-8")
    print(f"   ✅ Self-contained Krona HTML: {html_path}")
    print(f"   • Dataset: {_krona_dataset_name()}")
    print("   • Amount: abundance")
    print("   • Score: blank")
    print(f"   • Categories: {', '.join(used_ranks)}")
    return html_path


def build_krona_info_file(krona_dir: Path, krona_table: pd.DataFrame, per_sample_files: int, html_path: Path = None):
    """
    Write a short usage guide for KronaTools and Excel-template workflows.
    """
    used_ranks = krona_table.attrs.get("used_ranks", get_krona_ranks())
    collapse_label = str(COLLAPSE_STRATEGY).strip().lower()
    original_features = krona_table.attrs.get("original_features", "NA")
    grouped_rows_before_limit = krona_table.attrs.get("grouped_rows_before_limit", "NA")
    limit_note = krona_table.attrs.get("limit_note", "NA")

    info = f"""# MetaDiv Builder Krona Export

Collapse strategy suffix: {collapse_label}

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

Main Krona outputs:
- Krona_Total_Abundance_with_header_{collapse_label}.tsv
- {html_path.name if html_path is not None else "HTML not generated"}

Note:
Only one total Krona table is exported in v1.7.7 to avoid redundant files.
The exported table includes headers and the collapse-strategy suffix.

Per-sample Krona files:
- {per_sample_files}

Render with KronaTools, if desired:
-------------------------------
ktImportText Krona_Total_Abundance_with_header_{collapse_label}.tsv -o Krona_Total_Abundance_{collapse_label}.html

Recommended settings for very large meta-analyses:
--------------------------------------------------
KRONA_MODE = "GENUS"
KRONA_MAX_ROWS = 30000
KRONA_LIMIT_ACTION = "TOP"

Notes:
------
This module renders a self-contained Krona HTML directly from Python using
the same aggregated Krona table. The legacy Krona JavaScript and images are embedded
in MetaDiv Builder, so no KronaTools/Excel installation or extra Python dependency
is required. The TSV is retained as a transparent tabular record of the plotted data.
"""
    with open(krona_dir / f"KRONA_INFO_{collapse_label}.txt", "w", encoding="utf-8") as f:
        f.write(info)


def export_krona_package(final_df: pd.DataFrame, site_cols: list, output_dir: Path, report_lines: list):
    """
    Export Krona-compatible visualization tables.

    This module writes the standard aggregated table and renders a self-contained
    interactive Krona HTML directly from the same table. KronaTools and Excel are
    not required.
    """
    if not RUN_KRONA_EXPORT:
        print("\n🌐 KRONA EXPORT skipped (RUN_KRONA_EXPORT=False).")
        return

    print("\n🌐 KRONA EXPORT MODULE")
    print("=" * 70)

    krona_dir = output_dir / "Krona"
    krona_dir.mkdir(parents=True, exist_ok=True)

    krona_table = build_krona_total_table(final_df, site_cols, krona_dir)

    # Reuse the same grouped table to render HTML: no biological recalculation.
    html_path = build_krona_html(krona_table, krona_dir)

    per_sample_files = 0
    if RUN_KRONA_PER_SAMPLE_EXPORT:
        per_sample_files = build_krona_per_sample_tables(final_df, site_cols, krona_dir)

    build_krona_info_file(
        krona_dir, krona_table=krona_table, per_sample_files=per_sample_files, html_path=html_path
    )

    report_lines.append("\n  Krona Export:")
    report_lines.append(f"    Krona folder: {krona_dir.resolve()}")
    report_lines.append(f"    Collapse strategy suffix: {str(COLLAPSE_STRATEGY).strip().lower()}")
    report_lines.append(f"    KRONA_MODE: {KRONA_MODE}")
    report_lines.append(f"    Terminal rank: {krona_table.attrs.get('used_ranks', get_krona_ranks())[-1]}")
    report_lines.append(f"    Grouped rows before row-limit: {krona_table.attrs.get('grouped_rows_before_limit', 'NA')}")
    report_lines.append(f"    Total Krona rows exported: {len(krona_table):,}")
    report_lines.append(f"    Krona HTML: {html_path.resolve()}")
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
    RAM-safe For_R export.

    Reads FINAL_DB in row chunks and writes abundance/taxonomy/FASTA outputs
    incrementally. It never loads the complete wide final database into memory.
    """
    final_csv_path = Path(final_csv_path)
    print("\n" + "=" * 70)
    print(f"📦 FOR_R EXPORT: {suffix} (RAM-SAFE)")
    print("=" * 70)

    if not final_csv_path.exists():
        print(f"   ⚠️ Final database not found: {final_csv_path}")
        return

    header = pd.read_csv(final_csv_path, nrows=0)
    columns = list(header.columns)
    id_col = "SPPN" if "SPPN" in columns else ("OTU_XX" if "OTU_XX" in columns else None)
    if id_col is None:
        raise ValueError("For_R export requires either SPPN or OTU_XX column.")

    valid_site_cols = [c for c in site_cols if c in columns]
    sample_names = [str(c) for c in valid_site_cols]

    outdir = build_for_r_output_dir(suffix)
    outdir.mkdir(parents=True, exist_ok=True)

    abundance_path = outdir / "abundance_table.csv"
    taxonomy_table_path = outdir / "taxonomy_table.csv"
    metadata_path = outdir / "sample_metadata.csv"
    fasta_path = outdir / "sequences.fasta"

    tax_cols = [c for c in TAX_RANKS if c in columns]
    ecological_cols = [
        c for c in [
            "primary_lifestyle", "Secondary_lifestyle",
            "general_lifestyle", "secondary_lifestyle", "functional_group",
            "functional_confidence", "functional_source", "functional_match_rank",
        ] if c in columns
    ]
    taxonomy_cols = [id_col] + tax_cols + ecological_cols

    # Abundance export
    first = True
    feature_count = 0
    for chunk in pd.read_csv(
        final_csv_path,
        usecols=[id_col] + valid_site_cols,
        chunksize=RAM_SAFE_CHUNK_ROWS,
        low_memory=False
    ):
        if valid_site_cols:
            numeric = (
                chunk.loc[:, valid_site_cols]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(0)
                .astype(RAM_SAFE_ABUNDANCE_DTYPE, copy=False)
            )
            for c in valid_site_cols:
                chunk[c] = numeric[c]
            del numeric

        chunk.to_csv(
            abundance_path,
            mode="w" if first else "a",
            header=first,
            index=False,
            encoding="utf-8"
        )
        feature_count += len(chunk)
        first = False
        del chunk
        gc.collect()

    # Taxonomy export
    first = True
    for chunk in pd.read_csv(
        final_csv_path,
        usecols=taxonomy_cols,
        chunksize=RAM_SAFE_CHUNK_ROWS,
        low_memory=False
    ):
        chunk.to_csv(
            taxonomy_table_path,
            mode="w" if first else "a",
            header=first,
            index=False,
            encoding="utf-8"
        )
        first = False
        del chunk
        gc.collect()

    # Sample metadata template.
    metadata = pd.DataFrame({
        "sample_id": sample_names,
        "Group": "",
        "Site": "",
        "Treatment": "",
        "Latitude": "",
        "Longitude": "",
    })
    metadata.to_csv(metadata_path, index=False, encoding="utf-8")

    # FASTA export.
    nseq = 0
    with open(fasta_path, "w", encoding="utf-8") as handle:
        if "sequence" in columns:
            for chunk in pd.read_csv(
                final_csv_path,
                usecols=[id_col, "sequence"],
                chunksize=max(RAM_SAFE_CHUNK_ROWS, 10000),
                dtype=str
            ):
                for feature_id, seq in zip(
                    chunk[id_col].fillna("").astype(str),
                    chunk["sequence"].fillna("").astype(str)
                ):
                    feature_id = feature_id.strip()
                    seq = re.sub(r"\s+", "", seq.strip())
                    if not feature_id or feature_id.lower() in {"nan", "none"}:
                        continue
                    if not seq or seq.lower() in {"nan", "none"}:
                        continue
                    handle.write(f">{feature_id}\n")
                    for i in range(0, len(seq), 80):
                        handle.write(seq[i:i+80] + "\n")
                    nseq += 1
                del chunk
                gc.collect()

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
Features exported: {feature_count:,}
Samples exported: {len(valid_site_cols):,}

Files generated:
- abundance_table.csv
- taxonomy_table.csv
- sample_metadata.csv
- sequences.fasta

Important note about sample names:
For_R keeps the original sample names exactly as they appear in FINAL_DB.
When reading the tables in R, use check.names = FALSE.

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

    print(f"   ✅ abundance_table.csv: {feature_count:,} features x {len(valid_site_cols):,} samples")
    print(f"   ✅ taxonomy_table.csv: {feature_count:,} features")
    print(f"   ✅ sample_metadata.csv: {len(metadata):,} samples")
    print(f"   ✅ sequences.fasta: {nseq:,} sequences")
    print(f"   📁 For_R folder: {outdir.resolve()}")

def subdir_hint() -> str:
    """Return a human-readable description of the For_R folder scheme."""
    return f"For_R/<subset>/{build_for_r_tag(SUBSET_NAME)}/"


# ============================================================================
# REINPUT EXPORT PACKAGE
# ============================================================================

def _format_sintax_taxon(value) -> str:
    """Clean a taxon value for SINTAX-like export."""
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if value.lower() in {"", "nan", "none", "na", "null", "unclassified"}:
        return ""
    return value


def _format_sintax_pvalue(value, default: float = 1.0) -> float:
    """Return a safe probability value for SINTAX-like export."""
    try:
        if pd.isna(value):
            return default
        value = float(value)
        if value < 0:
            return default
        if value > 1:
            return 1.0
        return value
    except Exception:
        return default


def build_reinput_sintax_string(row: pd.Series) -> str:
    """
    Rebuild a SINTAX-like taxonomy string from FINAL_DB rank and p-value columns.

    Output example:
    d:Bacteria(1.0000),p:Proteobacteria(0.9800),...
    """
    rank_to_prefix = {
        "domain": "d",
        "phylum": "p",
        "class": "c",
        "order": "o",
        "family": "f",
        "genus": "g",
        "species": "s",
    }

    parts = []
    for rank in TAX_RANKS:
        taxon = _format_sintax_taxon(row.get(rank, ""))
        if not taxon:
            continue
        pval = _format_sintax_pvalue(row.get(f"{rank}_pvalue", 1.0), default=1.0)
        parts.append(f"{rank_to_prefix[rank]}:{taxon}({pval:.4f})")

    return ",".join(parts)


def _write_fasta_wrapped(records: pd.DataFrame, id_col: str, sequence_col: str, fasta_path: Path, wrap: int = 80) -> int:
    """Write sequences to FASTA and return the number of exported records."""
    nseq = 0
    with open(fasta_path, "w", encoding="utf-8") as f:
        for _, row in records.iterrows():
            seq_id = str(row[id_col]).strip()
            seq = row.get(sequence_col, "")
            if pd.isna(seq):
                seq = ""
            seq = re.sub(r"\s+", "", str(seq).strip())
            if not seq_id or not seq:
                continue
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i + wrap] + "\n")
            nseq += 1
    return nseq



def export_reinput_package(final_csv_path: Path, suffix: str, site_cols: list, report_lines: list | None = None):
    """
    RAM-safe ReInput export.

    Streams FINAL_DB in row chunks and writes abundance, SINTAX taxonomy and
    FASTA files without loading the complete wide database into memory.
    """
    if not RUN_REINPUT_EXPORT:
        return
    if not final_csv_path.exists():
        print(f"   ⚠️ ReInput skipped; FINAL_DB file not found: {final_csv_path}")
        return

    print(f"\n🔁 Exporting ReInput package for: {suffix} (RAM-SAFE)")

    header = pd.read_csv(final_csv_path, nrows=0)
    columns = list(header.columns)
    id_col = "SPPN" if "SPPN" in columns else columns[0]
    valid_site_cols = [c for c in site_cols if c in columns]
    if not valid_site_cols:
        print("   ⚠️ ReInput skipped; no sample abundance columns found.")
        return

    outdir = REINPUT_DIR / suffix / build_for_r_tag(suffix)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = f"MetaDivMaster_{MODE}_{build_for_r_tag(suffix)}"

    abundance_path = outdir / f"{prefix}_abundance.csv"
    taxonomy_path = outdir / f"{prefix}_taxonomy.sintax"
    sequences_path = outdir / f"{prefix}_sequences.fasta"
    info_path = outdir / "REINPUT_INFO.txt"

    first_abundance = True
    nrows = 0
    nseq = 0

    with open(taxonomy_path, "w", encoding="utf-8") as tax_handle, \
         open(sequences_path, "w", encoding="utf-8") as fasta_handle:

        for chunk in pd.read_csv(
            final_csv_path,
            chunksize=RAM_SAFE_CHUNK_ROWS,
            low_memory=False
        ):
            # Abundance
            abundance = chunk.loc[:, [id_col] + valid_site_cols].copy()
            abundance.rename(columns={id_col: "OTU"}, inplace=True)
            numeric = (
                abundance.loc[:, valid_site_cols]
                .apply(pd.to_numeric, errors="coerce")
                .fillna(0)
                .astype(RAM_SAFE_ABUNDANCE_DTYPE, copy=False)
            )
            for c in valid_site_cols:
                abundance[c] = numeric[c]
            abundance.to_csv(
                abundance_path,
                mode="w" if first_abundance else "a",
                header=first_abundance,
                index=False,
                encoding="utf-8"
            )
            first_abundance = False

            # Taxonomy and FASTA
            for _, row in chunk.iterrows():
                feature_id = str(row.get(id_col, "") or "").strip()
                if not feature_id or feature_id.lower() in {"nan", "none"}:
                    continue
                sintax = build_reinput_sintax_string(row)
                tax_handle.write(f"{feature_id}\t{sintax}\n")

                if "sequence" in chunk.columns:
                    seq = re.sub(r"\s+", "", str(row.get("sequence", "") or "").strip())
                    if seq and seq.lower() not in {"nan", "none"}:
                        fasta_handle.write(f">{feature_id}\n")
                        for i in range(0, len(seq), 80):
                            fasta_handle.write(seq[i:i+80] + "\n")
                        nseq += 1

            nrows += len(chunk)
            del chunk, abundance, numeric
            gc.collect()

    info = f"""# MetaDiv Builder - ReInput Package

Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
MODE: {MODE}
SUBSET_MODE: {SUBSET_MODE}
COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}
P_VALUE_THRESHOLD: {P_VALUE_THRESHOLD}
SPPN_P_THRESHOLD: {SPPN_P_THRESHOLD}

Source final database:
{final_csv_path}

Files generated:
- {abundance_path.name}
- {taxonomy_path.name}
- {sequences_path.name}

How to reuse:
Copy the three files above into input/{MODE}/ together with any new dataset files.
"""
    with open(info_path, "w", encoding="utf-8") as f:
        f.write(info)

    print(f"   ✅ ReInput abundance: {nrows:,} features x {len(valid_site_cols):,} samples")
    print(f"   ✅ ReInput taxonomy: {nrows:,} records")
    print(f"   ✅ ReInput sequences: {nseq:,} sequences")
    print(f"   📁 ReInput folder: {outdir.resolve()}")

    if report_lines is not None:
        report_lines.append(f"ReInput package exported for {suffix}: {outdir.resolve()}")

def _write_masked_dataframe_in_row_chunks(
    df: pd.DataFrame,
    mask,
    path: Path,
    columns: list,
    chunk_rows: int = None
):
    """Write a large row subset without materializing the entire subset."""
    if chunk_rows is None:
        chunk_rows = RAM_SAFE_CHUNK_ROWS

    mask_array = np.asarray(mask, dtype=bool)
    first = True
    written = 0

    for start in range(0, len(df), chunk_rows):
        stop = min(start + chunk_rows, len(df))
        local_mask = mask_array[start:stop]
        if not local_mask.any():
            continue
        chunk = df.iloc[start:stop]
        selected = chunk.loc[local_mask]
        selected.to_csv(
            path,
            mode="w" if first else "a",
            header=first,
            index=False,
            encoding="utf-8",
            columns=columns
        )
        written += len(selected)
        first = False
        del chunk, selected
        gc.collect()

    if first:
        pd.DataFrame(columns=columns).to_csv(path, index=False, encoding="utf-8")
    return written

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
    if MODE == "16S":
        # Keep only the two compact, user-relevant FAPROTAX annotations in
        # FINAL_DB. Diagnostic matching fields remain internal and are not
        # exported.
        lifestyle_cols = [
            col for col in ("faprotax_functions", "faprotax_function_count")
            if col in final_df.columns
        ]
    first_cols = (
        ["SPPN", "OTU_XX", "Original_ID"] + rank_cols + pval_cols +
        lifestyle_cols + ["sequence", "sequence_lenght", "Total_Abundance"]
    )
    existing = [c for c in first_cols if c in final_df.columns]
    # Defensive guard: never export the legacy exact-case ``Sequence`` field as
    # an abundance/sample column. The canonical lowercase ``sequence`` column
    # above is preserved unchanged.
    final_order = existing + [
        c for c in site_cols
        if c in final_df.columns and c != "Sequence"
    ]

    abundance_order = _abundance_order_positions(final_df)
    _write_dataframe_in_order_chunks(
        final_df,
        abundance_order,
        out_csv_main,
        final_order
    )
    print(f"💾 Final database saved (RAM-safe abundance order): {out_csv_main}")
    produced[suffix_main] = out_csv_main

    # For 16S with all_prokaryotes mode, also export only_bacteria subset
    if RUN_SECONDARY_SUBSET_EXPORTS and MODE == "16S" and SUBSET_MODE == "all_prokaryotes":
        bacteria_out = FINAL_DB_DIR / build_final_db_filename("only_bacteria")
        domain_series = final_df["domain"].astype(str).str.strip().str.lower()
        bacteria_mask = domain_series.eq("bacteria").to_numpy()
        bacteria_positions = abundance_order[bacteria_mask[abundance_order]]
        bacteria_n = _write_dataframe_in_order_chunks(
            final_df, bacteria_positions, bacteria_out, final_order
        )
        del bacteria_positions
        print(f"🦠 Bacteria-only subset saved: {bacteria_out} | Rows: {bacteria_n}")
        produced["only_bacteria"] = bacteria_out
        report_lines.append(f"\n  Bacteria-only subset: {bacteria_n:,} features")

    # For ITS with all_eukaryotes mode, also export only_fungi subset
    if RUN_SECONDARY_SUBSET_EXPORTS and MODE == "ITS" and SUBSET_MODE == "all_eukaryotes":
        fungi_out = FINAL_DB_DIR / build_final_db_filename("only_fungi")
        domain_series = final_df["domain"].astype(str).str.strip().str.lower()
        fungi_mask = domain_series.eq("fungi").to_numpy()
        fungi_positions = abundance_order[fungi_mask[abundance_order]]
        fungi_n = _write_dataframe_in_order_chunks(
            final_df, fungi_positions, fungi_out, final_order
        )
        del fungi_positions
        print(f"🍄 Fungi-only subset saved: {fungi_out} | Rows: {fungi_n}")
        produced["only_fungi"] = fungi_out
        report_lines.append(f"\n  Fungi-only subset: {fungi_n:,} features")


    # For CO1 with all_eukaryotes mode, also export only_metazoa subset
    if RUN_SECONDARY_SUBSET_EXPORTS and MODE == "CO1" and SUBSET_MODE == "all_eukaryotes":
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
        metazoa_mask_array = metazoa_mask.to_numpy()
        metazoa_positions = abundance_order[metazoa_mask_array[abundance_order]]
        metazoa_n = _write_dataframe_in_order_chunks(
            final_df, metazoa_positions, metazoa_out, final_order
        )
        del metazoa_positions, metazoa_mask_array
        print(f"🧬 CO1 Metazoa-only subset saved: {metazoa_out} | Rows: {metazoa_n}")
        produced["only_metazoa"] = metazoa_out
        report_lines.append(f"\n  CO1 Metazoa-only subset: {metazoa_n:,} features")

    del abundance_order
    gc.collect()

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

def _report_stage_time(report_lines, stage_name, started_at):
    """Record and print elapsed wall time for one pipeline stage."""
    elapsed = time.perf_counter() - started_at
    if RUN_STAGE_TIMING:
        message = f"STAGE TIME — {stage_name}: {elapsed:.2f} seconds ({elapsed/3600:.3f} h)"
        report_lines.append(message)
        print(f"⏱️ {message}")
    return elapsed


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    pipeline_start = time.perf_counter()
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
    report_lines.append(f"RAM_SAFE_MODE: {RAM_SAFE_MODE}")
    report_lines.append(f"RAM_SAFE_CHUNK_ROWS: {RAM_SAFE_CHUNK_ROWS}")
    report_lines.append(f"RAM_SAFE_SITE_CHUNK: {RAM_SAFE_SITE_CHUNK}")
    report_lines.append(f"RAM_SAFE_GROUP_ROWS: {RAM_SAFE_GROUP_ROWS}")
    report_lines.append(f"RUN_STAGE_TIMING: {RUN_STAGE_TIMING}")
    report_lines.append(f"DEV_MODE: {DEV_MODE} (temporary folder will be kept only if True)")
    report_lines.append(f"RUN_FUNGALTRAITS_ITS: {RUN_FUNGALTRAITS_ITS}")
    report_lines.append(f"RUN_FAPROTAX_16S: {RUN_FAPROTAX_16S}")
    report_lines.append("RUN_CO1_FUNCTIONAL_PLACEHOLDER: True if MODE == 'CO1' else False")
    report_lines.append(f"RUN_FOR_R_EXPORT: {RUN_FOR_R_EXPORT}")
    report_lines.append(f"RUN_SECONDARY_SUBSET_EXPORTS: {RUN_SECONDARY_SUBSET_EXPORTS}")
    report_lines.append(f"RUN_REINPUT_EXPORT: {RUN_REINPUT_EXPORT}")
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

    # Remove stale active-target intermediates from previous runs. This is
    # restricted to TARGET_DIR/PATTERN and prevents Part 2 from mixing old
    # datasets with the current run.
    stale_files = list(TARGET_DIR.glob(PATTERN))
    if stale_files:
        print(f"\n🧹 Removing {len(stale_files)} stale Part 2 input file(s)...")
        for stale_file in stale_files:
            stale_file.unlink()

    all_stats = []

    for prefix, files in datasets.items():
        try:
            stats = process_single_dataset(prefix, files, report_lines)
            all_stats.append(stats)
        except Exception as e:
            print(f"❌ Error processing dataset {prefix}: {e}")
            report_lines.append(f"\n  ERROR processing {prefix}: {e}")
            # Fail fast: never allow Part 2 to reuse stale intermediate files
            # after a failed Part 1 dataset.
            raise RuntimeError(
                f"Part 1 failed for dataset '{prefix}'. Pipeline stopped to "
                "prevent reuse of stale intermediate files."
            ) from e

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
    stage_started = time.perf_counter()
    combined, site_cols, part2_input_dir, part2_pattern, abundance_matrix = load_datasets_with_stats_part2(report_lines)
    _report_stage_time(report_lines, "Part 2 load", stage_started)

    if combined.empty:
        print("❌ No data to collapse. Check that Part 1 produced files.")
        return

    stage_started = time.perf_counter()
    final_df, logs = collapse_otus_simple(
        combined,
        site_cols,
        report_lines,
        abundance_matrix=abundance_matrix
    )
    _report_stage_time(report_lines, "Taxonomic collapse", stage_started)
    # The combined pre-collapse matrix is no longer needed. Releasing it here
    # is essential before sorting, taxonomy sanitation and functional ecology.
    del combined, abundance_matrix
    gc.collect()

    stage_started = time.perf_counter()
    final_df = compute_total_abundance_and_sort(final_df, site_cols)
    final_df = sanitize_taxonomy_columns_by_threshold(final_df, threshold=SPPN_P_THRESHOLD, cascade_blank=True)

    # Integrate fungal lifestyles if MODE is ITS.
    # The compact FINAL_DB receives only primary_lifestyle and Secondary_lifestyle.
    # A complete FungalTraits export is created later in Functional_Ecology.
    traits_df = pd.DataFrame()
    if MODE == "ITS" and FUNGAL_TRAITS_DB is not None:
        traits_df = load_fungal_traits_db(FUNGAL_TRAITS_DB)
        if not traits_df.empty and "genus" in final_df.columns:
            # RAM-safe mapping avoids merging/copying the complete wide abundance matrix.
            genus_clean = final_df["genus"].astype(str).str.strip().str.lower()
            traits_unique = (
                traits_df[["genus_clean", "primary_lifestyle", "Secondary_lifestyle"]]
                .drop_duplicates(subset=["genus_clean"], keep="first")
                .set_index("genus_clean")
            )
            final_df["primary_lifestyle"] = genus_clean.map(
                traits_unique["primary_lifestyle"]
            ).fillna("")
            final_df["Secondary_lifestyle"] = genus_clean.map(
                traits_unique["Secondary_lifestyle"]
            ).fillna("")
            del genus_clean, traits_unique
            gc.collect()
            matched = final_df["primary_lifestyle"].ne("").sum()
            print(f"   • Lifestyle info added for {matched} rows ({matched/len(final_df)*100:.1f}%)")
        else:
            # Add empty columns if no traits available or no genus column
            final_df["primary_lifestyle"] = ""
            final_df["Secondary_lifestyle"] = ""
            print("   • No lifestyle data added (missing traits or genus column)")
    else:
        # 16S functional ecology is added after zero pruning and SPPN assignment
        # so module-specific outputs include final stable identifiers.
        pass

    final_df = prune_zero_abundance(final_df, site_cols)
    final_df = assign_sppn_unique(final_df)
    _report_stage_time(report_lines, "Post-collapse totals/taxonomy/SPPN", stage_started)

    stage_started = time.perf_counter()
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
    _report_stage_time(report_lines, "Functional ecology", stage_started)

    stage_started = time.perf_counter()
    produced_csvs = save_outputs(final_df, logs, site_cols, part2_input_dir, part2_pattern, report_lines)
    _report_stage_time(report_lines, "FINAL_DB exports", stage_started)

    # KRONA: Export Krona-compatible visualization tables
    stage_started = time.perf_counter()
    export_krona_package(final_df, site_cols, OUTPUT_DIR, report_lines)
    _report_stage_time(report_lines, "Krona export", stage_started)

    # All remaining exports read the saved FINAL_DB in RAM-safe row chunks.
    # Release the large in-memory database before those steps.
    final_feature_count = len(final_df)
    del final_df
    gc.collect()
    _cleanup_ram_safe_temp_files()

    # PART 3: Export For_R packages
    print("\n🧬 PART 3 — EXPORTING For_R PACKAGES")
    stage_started = time.perf_counter()

    if RUN_FOR_R_EXPORT:
        if SUBSET_NAME in produced_csvs:
            export_for_r_package(produced_csvs[SUBSET_NAME], suffix=SUBSET_NAME, site_cols=site_cols)

        if MODE == "16S" and "only_bacteria" in produced_csvs:
            export_for_r_package(produced_csvs["only_bacteria"], suffix="only_bacteria", site_cols=site_cols)

        if MODE == "ITS" and "only_fungi" in produced_csvs:
            export_for_r_package(produced_csvs["only_fungi"], suffix="only_fungi", site_cols=site_cols)

        if MODE == "CO1" and "only_metazoa" in produced_csvs:
            export_for_r_package(produced_csvs["only_metazoa"], suffix="only_metazoa", site_cols=site_cols)
    else:
        print("   ⏭️ For_R export disabled by USER CONFIGURATION.")

    _report_stage_time(report_lines, "For_R exports", stage_started)

    # PART 4: Export ReInput packages
    print("\n🔁 PART 4 — EXPORTING ReInput PACKAGES")
    stage_started = time.perf_counter()
    for suffix, csv_path in produced_csvs.items():
        export_reinput_package(csv_path, suffix=suffix, site_cols=site_cols, report_lines=report_lines)
    _report_stage_time(report_lines, "ReInput exports", stage_started)

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
    report_lines.append(f"Final features after collapse and pruning: {final_feature_count:,}")
    report_lines.append(f"FINAL_DB folder: {FINAL_DB_DIR.resolve()}")
    report_lines.append(f"For_R folder: {FOR_R_DIR.resolve()}")
    if RUN_REINPUT_EXPORT:
        report_lines.append(f"ReInput folder: {REINPUT_DIR.resolve()}")
    if (MODE == "ITS" and RUN_FUNGALTRAITS_ITS) or (MODE == "16S" and RUN_FAPROTAX_16S) or MODE == "CO1":
        report_lines.append(f"Functional_Ecology folder: {(OUTPUT_DIR / 'Functional_Ecology').resolve()}")
    report_lines.append(f"Total execution time: {time.perf_counter() - pipeline_start:.2f} seconds")
    report_lines.append("=" * 70)
    report_lines.append("\n✅ MetaDiv Builder pipeline completed successfully")

    report_path = OUTPUT_DIR / build_report_filename(SUBSET_NAME)
    write_report(report_lines, report_path)

    print("\n🎉 ALL DONE")
    print(f"   MODE: {MODE}")
    print(f"   SUBSET_MODE: {SUBSET_MODE}")
    print(f"   COLLAPSE_STRATEGY: {COLLAPSE_STRATEGY}")
    print(f"   Final features: {final_feature_count:,}")
    print(f"   FINAL_DB folder: {FINAL_DB_DIR.resolve()}")
    print(f"   For_R folder: {FOR_R_DIR.resolve()}")
    print(f"   Report saved: {report_path.resolve()}")

if __name__ == "__main__":
    main()

