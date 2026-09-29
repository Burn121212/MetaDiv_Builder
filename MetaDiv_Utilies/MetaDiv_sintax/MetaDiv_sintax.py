#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MetaDiv_sintax
VSEARCH/SINTAX taxonomic classification utility for MetaDiv Builder.

This CLI version is an alternative to the manual user-settings block with command-line
arguments. VSEARCH can be run through Docker (default) or through a Conda
environment such as MetaDiv_sintax.
"""

from pathlib import Path
from datetime import datetime
import argparse
import json
import os
import shlex
import shutil
import subprocess

# ============================================================
# INTERNAL HELPERS
# ============================================================

FASTA_EXTENSIONS = {
    ".fasta",
    ".fa",
    ".fas",
    ".fna",
    ".ffn",
}


def normalize_mode(mode):
    """
    Normalize the marker name used by MetaDiv Builder.
    """
    value = str(mode).strip().upper()

    aliases = {
        "ITS": "ITS",
        "16S": "16S",
        "COI": "COI",
        "CO1": "COI",
    }

    if value not in aliases:
        raise ValueError(
            "MODE must be one of: ITS, 16S, COI."
        )

    return aliases[value]


def find_project_root(configured_root=None):
    """
    Detect the MetaDiv Builder project root automatically, or use a user path.

    The search starts from the current working directory and moves upward until
    both of the following are found:

        input/
        databases/taxonomic_reference_db/
    """
    if configured_root is not None:
        project_root = Path(configured_root).expanduser().resolve()
        if not project_root.is_dir():
            raise FileNotFoundError(
                f"\nProject root does not exist or is not a directory:\n{project_root}"
            )
        return project_root

    current = Path.cwd().resolve()

    for candidate in [current] + list(current.parents):

        input_dir = candidate / "input"

        reference_dir = (
            candidate
            / "databases"
            / "taxonomic_reference_db"
        )

        if input_dir.is_dir() and reference_dir.is_dir():
            return candidate

    raise FileNotFoundError(
        "\nMetaDiv Builder project root could not be detected.\n\n"
        "Run from inside a project containing:\n\n"
        "  input/\n"
        "  databases/taxonomic_reference_db/\n\n"
        "or provide --project-root explicitly.\n"
    )


def locate_marker_input(project_root, mode):
    """
    Locate the normal MetaDiv input folder for the selected marker.
    """
    candidates = [
        project_root / "input" / mode,
    ]

    # Backward compatibility with repositories using CO1.
    if mode == "COI":
        candidates.append(
            project_root / "input" / "CO1"
        )

    for directory in candidates:
        if directory.is_dir():
            return directory

    expected = project_root / "input" / mode

    raise FileNotFoundError(
        f"\nMarker input directory not found:\n{expected}\n"
    )


def resolve_input_fastas(input_dir, file_names):
    """
    Resolve and validate the FASTA filenames selected by the user.
    """
    if not isinstance(file_names, (list, tuple)) or len(file_names) == 0:
        raise ValueError(
            "At least one FASTA file must be provided with --fasta-files, or use --all-fastas."
        )

    resolved = []

    for name in file_names:

        fasta = (input_dir / name).resolve()

        try:
            fasta.relative_to(input_dir.resolve())
        except ValueError:
            raise ValueError(
                f"FASTA file must be located inside:\n{input_dir}\n\n"
                f"Invalid path: {name}"
            )

        if not fasta.is_file():
            raise FileNotFoundError(
                f"\nFASTA file not found:\n{fasta}"
            )

        if fasta.suffix.lower() not in FASTA_EXTENSIONS:
            raise ValueError(
                f"\nUnsupported FASTA extension:\n{fasta.name}\n\n"
                f"Accepted extensions: {sorted(FASTA_EXTENSIONS)}"
            )

        resolved.append(fasta)

    return resolved


def output_path_for_fasta(fasta_file):
    """
    Create the companion SINTAX output filename.

    Example:
        ATLASMXB.fasta
        -> ATLASMXB_taxonomy.sintax
    """
    if fasta_file.stem.endswith('_sequences'):
        newname = fasta_file.with_name(
        f"{fasta_file.stem}".removesuffix('_sequences')
    )
    else:
        newname = fasta_file

        
    return newname.with_name(
        f"{newname.stem}_taxonomy.sintax"
    )


def pretty_command(command):
    """
    Format a command for the processing report.
    """
    try:
        return shlex.join([str(x) for x in command])
    except Exception:
        return " ".join(str(x) for x in command)


def run_command(
    command,
    check=True,
    capture_output=True,
):
    """
    Execute a system command without using a shell.
    """
    result = subprocess.run(
        [str(x) for x in command],
        check=False,
        text=True,
        capture_output=capture_output,
    )

    if check and result.returncode != 0:

        message = [
            f"Command failed with exit code {result.returncode}.",
            "",
            "Command:",
            pretty_command(command),
        ]

        if result.stdout:
            message.extend(
                [
                    "",
                    "STDOUT:",
                    result.stdout,
                ]
            )

        if result.stderr:
            message.extend(
                [
                    "",
                    "STDERR:",
                    result.stderr,
                ]
            )

        raise RuntimeError(
            "\n".join(message)
        )

    return result


def count_fasta_sequences(fasta_file):
    """
    Count FASTA records using header lines beginning with '>'.
    """
    count = 0

    with fasta_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith(">"):
                count += 1

    return count


def validate_sintax_output(output_file):
    """
    Summarize the raw VSEARCH --tabbedout output.

    Column 1: query identifier
    Column 2: taxonomic prediction with bootstrap support
    Column 3: strand
    Column 4: cutoff-filtered taxonomy when available
    """
    total_records = 0
    taxonomy_predictions = 0
    cutoff_predictions = 0

    with output_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            line = line.rstrip("\r\n")

            if not line:
                continue

            fields = line.split("\t")

            total_records += 1

            if (
                len(fields) >= 2
                and fields[1].strip() not in {"", "*"}
            ):
                taxonomy_predictions += 1

            if (
                len(fields) >= 4
                and fields[3].strip() not in {"", "*"}
            ):
                cutoff_predictions += 1

    return {
        "records": total_records,
        "predictions": taxonomy_predictions,
        "cutoff_predictions": cutoff_predictions,
    }



# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================


def probability(value):
    """Parse a float in the closed interval [0, 1]."""
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Expected a number between 0 and 1, got: {value}")
    if not 0.0 <= value <= 1.0:
        raise argparse.ArgumentTypeError(f"Expected a number between 0 and 1, got: {value}")
    return value


def positive_int(value):
    """Parse a positive integer."""
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Expected a positive integer, got: {value}")
    if value < 1:
        raise argparse.ArgumentTypeError(f"Expected a positive integer, got: {value}")
    return value


def path_argument(value):
    """Expand environment variables and '~' in path-like arguments."""
    return Path(os.path.expandvars(str(value))).expanduser()


def discover_all_fastas(input_dir):
    """Return all supported FASTA files under input_dir as paths relative to input_dir."""
    fastas = []
    for fasta in sorted(input_dir.rglob("*")):
        if fasta.is_file() and fasta.suffix.lower() in FASTA_EXTENSIONS:
            fastas.append(str(fasta.relative_to(input_dir)))
    if not fastas:
        raise FileNotFoundError(
            f"\nNo FASTA files with extensions {sorted(FASTA_EXTENSIONS)} found in:\n{input_dir}"
        )
    return fastas


class ArgumentFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter,
):
    """Keep example formatting and display default values in --help."""
    pass


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="MetaDiv_sintax.py",
        formatter_class=ArgumentFormatter,
        description=(
            "Run VSEARCH --sintax on MetaDiv Builder FASTA files and write "
            "<sample>_sintax.txt outputs beside each input FASTA."
        ),
        epilog="""
Examples:

  Docker runtime, one ITS FASTA file:
    python MetaDiv_sintax.py \
      --mode ITS \
      --fasta-files ATLASMXBC_red_sequences.fasta \
      --reference-db SINTAX_EUKARYOME_ITS_v2.0.fasta \
      --runtime docker \
      --threads 16 \
      --overwrite-existing

  Conda runtime, automatically create/use environment MetaDiv_sintax:
    python MetaDiv_sintax.py \
      --mode ITS \
      --fasta-files ATLASMXBC_red_sequences.fasta \
      --reference-db SINTAX_EUKARYOME_ITS_v2.0.fasta \
      --runtime conda \
      --conda-env-name MetaDiv_sintax \
      --auto-create-conda-env

  Classify every FASTA under input/COI using Docker:
    python MetaDiv_sintax.py \
      --mode COI \
      --all-fastas \
      --reference-db SINTAX_MIDORI2_LONGEST_NUC_GB271_CO1.udb \
      --runtime docker

Notes:
  - FASTA paths supplied with --fasta-files are interpreted relative to input/<MODE>/
    unless --input-dir is provided.
  - The reference database is interpreted relative to
    databases/taxonomic_reference_db/ unless --reference-db-dir is provided.
  - Use --runtime conda when Docker is unavailable or when you prefer VSEARCH from
    a Conda environment.
""",
    )

    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        type=normalize_mode,
        help="Genetic marker. Options: ITS, 16S, COI. CO1 is accepted as an alias for COI.",
    )

    fasta_group = parser.add_mutually_exclusive_group(required=True)
    fasta_group.add_argument(
        "-f",
        "--fasta-files",
        nargs="+",
        help=(
            "One or more FASTA filenames to classify. Paths are relative to the "
            "selected input directory. Subfolders are accepted."
        ),
    )
    fasta_group.add_argument(
        "--all-fastas",
        action="store_true",
        help="Classify every supported FASTA file found recursively in the selected input directory.",
    )

    parser.add_argument(
        "-d",
        "--reference-db",
        required=True,
        help=(
            "Reference SINTAX database filename, usually stored inside "
            "databases/taxonomic_reference_db/."
        ),
    )

    parser.add_argument(
        "--project-root",
        type=path_argument,
        default=None,
        help=(
            "MetaDiv Builder project root. If omitted, the script searches upward "
            "from the current directory."
        ),
    )

    parser.add_argument(
        "--input-dir",
        type=path_argument,
        default=None,
        help="Custom input directory. If omitted, input/<MODE>/ is used.",
    )

    parser.add_argument(
        "--reference-db-dir",
        type=path_argument,
        default=None,
        help="Custom taxonomic reference database directory.",
    )

    parser.add_argument(
        "--sintax-cutoff",
        type=probability,
        default=0.8,
        help="VSEARCH --sintax_cutoff value.",
    )

    parser.add_argument(
        "--strand",
        choices=["plus", "both"],
        default="both",
        help="VSEARCH --strand option.",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        default=16,
        help="Number of VSEARCH threads.",
    )

    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Overwrite existing <sample>_sintax.txt files.",
    )

    parser.add_argument(
        "--skip-existing",
        dest="overwrite_existing",
        action="store_false",
        help="Skip FASTA files that already have a <sample>_sintax.txt output.",
    )
    parser.set_defaults(overwrite_existing=True)

    parser.add_argument(
        "--runtime",
        choices=["docker", "conda"],
        default="docker",
        help="Runtime used to execute VSEARCH.",
    )

    parser.add_argument(
        "--docker-image",
        default="pipecraft/vsearch:2.30.4-pc1.2.0",
        help="Docker image used when --runtime docker.",
    )

    parser.add_argument(
        "--auto-pull-image",
        dest="auto_pull_image",
        action="store_true",
        help="Automatically pull the Docker image if it is not available locally.",
    )

    parser.add_argument(
        "--no-auto-pull-image",
        dest="auto_pull_image",
        action="store_false",
        help="Do not pull a missing Docker image automatically.",
    )
    parser.set_defaults(auto_pull_image=True)

    parser.add_argument(
        "--conda-env-name",
        default="MetaDiv_sintax",
        help="Conda environment used when --runtime conda.",
    )

    parser.add_argument(
        "--auto-create-conda-env",
        action="store_true",
        help=(
            "When --runtime conda, create or repair the Conda environment with "
            "VSEARCH if it is missing or VSEARCH cannot be executed."
        ),
    )

    parser.add_argument(
        "--vsearch-package",
        default="vsearch",
        help=(
            "Conda package spec used with --auto-create-conda-env. Examples: "
            "vsearch, vsearch=2.30.0."
        ),
    )

    parser.add_argument(
        "--report-file",
        type=path_argument,
        default=None,
        help="Processing report path. Default: ./MetaDiv_sintax_report.txt.",
    )

    return parser.parse_args()


def conda_env_exists(env_name):
    """Return True if a named Conda environment exists."""
    result = run_command(
        ["conda", "env", "list", "--json"],
        check=True,
    )
    data = json.loads(result.stdout)
    for env_path in data.get("envs", []):
        if Path(env_path).name == env_name:
            return True
    return False


def install_vsearch_conda_env(env_name, package_spec):
    """Create or update the Conda environment used for VSEARCH."""
    if conda_env_exists(env_name):
        print(f"\nConda environment found: {env_name}")
        print(f"Installing/updating {package_spec} in {env_name} ...")
        command = [
            "conda", "install", "-y",
            "-n", env_name,
            "-c", "conda-forge",
            "-c", "bioconda",
            package_spec,
        ]
    else:
        print(f"\nConda environment not found: {env_name}")
        print(f"Creating {env_name} with {package_spec} ...")
        command = [
            "conda", "create", "-y",
            "-n", env_name,
            "-c", "conda-forge",
            "-c", "bioconda",
            package_spec,
        ]

    run_command(command, capture_output=False)


def verify_or_prepare_conda_vsearch(env_name, auto_create, package_spec):
    """Verify VSEARCH through Conda, optionally creating/updating the env."""
    if shutil.which("conda") is None:
        raise RuntimeError(
            "\nConda was not detected in PATH. Install/configure Conda or use --runtime docker."
        )

    test = run_command(
        ["conda", "run", "-n", env_name, "vsearch", "--version"],
        check=False,
    )

    if test.returncode == 0:
        return test

    if not auto_create:
        raise RuntimeError(
            f"\nVSEARCH could not be executed from Conda environment '{env_name}'.\n\n"
            "Either create the environment manually, for example:\n\n"
            f"  conda create -n {env_name} -c conda-forge -c bioconda {package_spec}\n\n"
            "or rerun with:\n\n"
            f"  --runtime conda --conda-env-name {env_name} --auto-create-conda-env\n"
        )

    install_vsearch_conda_env(env_name, package_spec)

    return run_command(
        ["conda", "run", "-n", env_name, "vsearch", "--version"],
        check=True,
    )


def check_docker_runtime(docker_image, auto_pull_image):
    """Validate Docker and ensure the requested image is available."""
    if shutil.which("docker") is None:
        raise RuntimeError(
            "\nDocker CLI was not detected.\n\n"
            "Install/start Docker Desktop or rerun with --runtime conda."
        )

    docker_info = run_command(
        ["docker", "info"],
        check=False,
    )

    if docker_info.returncode != 0:
        raise RuntimeError(
            "\nDocker was detected, but Docker Desktop/Engine does not appear to be running.\n"
            "Start Docker or rerun with --runtime conda."
        )

    image_check = run_command(
        ["docker", "image", "inspect", docker_image],
        check=False,
    )

    if image_check.returncode != 0:
        if not auto_pull_image:
            raise RuntimeError(
                f"\nDocker image not found locally:\n{docker_image}\n\n"
                "Rerun with --auto-pull-image or pull it manually."
            )

        print(f"\nDocker image not found locally. Pulling {docker_image} ...")
        run_command(["docker", "pull", docker_image], capture_output=False)

    return run_command(
        [
            "docker", "run", "--rm",
            "--entrypoint", "vsearch",
            docker_image,
            "--version",
        ],
        check=True,
    )


def build_vsearch_command(fasta_file, output_file):
    """Build the VSEARCH --sintax command for the selected runtime."""
    if RUNTIME == "docker":
        working_directory = fasta_file.parent
        container_fasta = f"/work/{fasta_file.name}"
        container_output = f"/work/{output_file.name}"
        container_database = f"/db/{REFERENCE_DB.name}"

        return [
            "docker", "run", "--rm",
            "--mount", f"type=bind,source={working_directory},target=/work",
            "--mount", f"type=bind,source={REFERENCE_DB_DIR},target=/db,readonly",
            "--entrypoint", "vsearch",
            DOCKER_IMAGE,
            "--sintax", container_fasta,
            "--db", container_database,
            "--tabbedout", container_output,
            "--sintax_cutoff", str(SINTAX_CUTOFF),
            "--strand", STRAND,
            "--threads", str(THREADS),
        ]

    if RUNTIME == "conda":
        return [
            "conda", "run", "-n", CONDA_ENV_NAME,
            "vsearch",
            "--sintax", str(fasta_file),
            "--db", str(REFERENCE_DB),
            "--tabbedout", str(output_file),
            "--sintax_cutoff", str(SINTAX_CUTOFF),
            "--strand", STRAND,
            "--threads", str(THREADS),
        ]

    raise ValueError(f"Unsupported runtime: {RUNTIME}")


# ============================================================
# PARSE / VALIDATE COMMAND-LINE SETTINGS
# ============================================================

ARGS = parse_arguments()

MODE = ARGS.mode
FASTA_FILES = ARGS.fasta_files
REFERENCE_DB_NAME = ARGS.reference_db
SINTAX_CUTOFF = ARGS.sintax_cutoff
STRAND = ARGS.strand
THREADS = ARGS.threads
OVERWRITE_EXISTING = ARGS.overwrite_existing
RUNTIME = ARGS.runtime
DOCKER_IMAGE = ARGS.docker_image
AUTO_PULL_IMAGE = ARGS.auto_pull_image
CONDA_ENV_NAME = ARGS.conda_env_name
AUTO_CREATE_CONDA_ENV = ARGS.auto_create_conda_env
VSEARCH_PACKAGE = ARGS.vsearch_package

# The type/choices in argparse already validate these, but keep explicit
# safeguards close to the legacy validation messages.
MODE = normalize_mode(MODE)

if not 0.0 <= float(SINTAX_CUTOFF) <= 1.0:
    raise ValueError("SINTAX_CUTOFF must be between 0.0 and 1.0.")

if STRAND not in {"plus", "both"}:
    raise ValueError("STRAND must be 'plus' or 'both'.")

if int(THREADS) < 1:
    raise ValueError("THREADS must be >= 1.")


# ============================================================
# DETECT METADIV PATHS
# ============================================================

PROJECT_ROOT = find_project_root(ARGS.project_root)

if ARGS.input_dir is not None:
    INPUT_DIR = ARGS.input_dir.expanduser().resolve()
    if not INPUT_DIR.is_dir():
        raise FileNotFoundError(f"\nInput directory not found:\n{INPUT_DIR}")
else:
    INPUT_DIR = locate_marker_input(PROJECT_ROOT, MODE)

if ARGS.reference_db_dir is not None:
    REFERENCE_DB_DIR = ARGS.reference_db_dir.expanduser().resolve()
    if not REFERENCE_DB_DIR.is_dir():
        raise FileNotFoundError(f"\nReference database directory not found:\n{REFERENCE_DB_DIR}")
else:
    REFERENCE_DB_DIR = (
        PROJECT_ROOT
        / "databases"
        / "taxonomic_reference_db"
    )

REFERENCE_DB = (
    REFERENCE_DB_DIR
    / REFERENCE_DB_NAME
).resolve()

if not REFERENCE_DB.is_file():
    raise FileNotFoundError(
        "\nReference database not found:\n"
        f"{REFERENCE_DB}\n\n"
        "Check --reference-db and --reference-db-dir."
    )

if ARGS.all_fastas:
    FASTA_FILES = discover_all_fastas(INPUT_DIR)

INPUT_FASTAS = resolve_input_fastas(
    INPUT_DIR,
    FASTA_FILES,
)

REPORT_FILE = (
    ARGS.report_file.expanduser().resolve()
    if ARGS.report_file is not None
    else Path.cwd() / "MetaDiv_sintax_report.txt"
)


# ============================================================
# PREPARE / VERIFY VSEARCH RUNTIME
# ============================================================

if RUNTIME == "docker":
    version_result = check_docker_runtime(
        docker_image=DOCKER_IMAGE,
        auto_pull_image=AUTO_PULL_IMAGE,
    )
elif RUNTIME == "conda":
    version_result = verify_or_prepare_conda_vsearch(
        env_name=CONDA_ENV_NAME,
        auto_create=AUTO_CREATE_CONDA_ENV,
        package_spec=VSEARCH_PACKAGE,
    )
else:
    raise ValueError(f"Unsupported runtime: {RUNTIME}")

VSEARCH_VERSION = (
    version_result.stdout
    or version_result.stderr
    or "unknown"
).strip()


# ============================================================
# RUN SINTAX
# ============================================================

print("=" * 72)
print("MetaDiv_sintax")
print("=" * 72)
print(f"MODE: {MODE}")
print(f"Runtime: {RUNTIME}")
if RUNTIME == "docker":
    print(f"Docker image: {DOCKER_IMAGE}")
if RUNTIME == "conda":
    print(f"Conda environment: {CONDA_ENV_NAME}")
print(f"Reference database: {REFERENCE_DB.name}")
print(f"SINTAX cutoff: {SINTAX_CUTOFF}")
print(f"Strand: {STRAND}")
print(f"Threads: {THREADS}")
print(f"Input directory: {INPUT_DIR}")
print(f"FASTA files selected: {len(INPUT_FASTAS)}")
print("=" * 72)


run_started = datetime.now()

RUN_SUMMARY = []

for fasta_file in INPUT_FASTAS:

    output_file = output_path_for_fasta(
        fasta_file
    )

    # --------------------------------------------------------
    # Skip existing outputs unless overwrite is enabled
    # --------------------------------------------------------

    if output_file.exists() and not OVERWRITE_EXISTING:

        print("\n" + "-" * 72)
        print(f"SKIPPED: {fasta_file.name}")
        print(f"Existing output: {output_file.name}")
        print(
            "Rerun with --overwrite-existing to classify this file again."
        )
        print("-" * 72)

        RUN_SUMMARY.append(
            {
                "input": fasta_file,
                "output": output_file,
                "status": "SKIPPED_EXISTING",
                "sequences": count_fasta_sequences(fasta_file),
                "seconds": 0.0,
                "validation": validate_sintax_output(output_file),
                "command": None,
            }
        )

        continue


    # --------------------------------------------------------
    # Build VSEARCH command for the selected runtime
    # --------------------------------------------------------

    command = build_vsearch_command(
        fasta_file=fasta_file,
        output_file=output_file,
    )


    sequence_count = count_fasta_sequences(
        fasta_file
    )


    print("\n" + "-" * 72)
    print(f"Classifying: {fasta_file.name}")
    print(f"Representative sequences: {sequence_count:,}")
    print(f"Output: {output_file.name}")
    print("-" * 72)


    start_time = datetime.now()

    run_command(
        command,
        check=True,
    )

    end_time = datetime.now()

    elapsed = (
        end_time - start_time
    ).total_seconds()


    if not output_file.is_file():
        raise RuntimeError(
            "\nVSEARCH completed but the expected SINTAX output "
            "was not created:\n"
            f"{output_file}"
        )


    validation = validate_sintax_output(
        output_file
    )


    print(f"SINTAX records: {validation['records']:,}")
    print(
        f"Taxonomy predictions: "
        f"{validation['predictions']:,}"
    )
    print(
        f"Cutoff-filtered predictions: "
        f"{validation['cutoff_predictions']:,}"
    )
    print(f"Execution time: {elapsed:.2f} s")


    RUN_SUMMARY.append(
        {
            "input": fasta_file,
            "output": output_file,
            "status": "CLASSIFIED",
            "sequences": sequence_count,
            "seconds": elapsed,
            "validation": validation,
            "command": command,
        }
    )


run_finished = datetime.now()

TOTAL_SECONDS = (
    run_finished - run_started
).total_seconds()


# ============================================================
# WRITE PROCESSING REPORT
# ============================================================

report = [
    "=" * 72,
    "MetaDiv_sintax - Processing Report",
    "Author: Bernardo Águila, UNAM",
    f"Date: {run_finished.isoformat(timespec='seconds')}",
    "=" * 72,
    f"MODE: {MODE}",
    f"Runtime: {RUNTIME}",
    f"Docker image: {DOCKER_IMAGE}" if RUNTIME == "docker" else f"Conda environment: {CONDA_ENV_NAME}",
    f"VSEARCH version: {VSEARCH_VERSION}",
    f"Reference database: {REFERENCE_DB}",
    f"SINTAX_CUTOFF: {SINTAX_CUTOFF}",
    f"STRAND: {STRAND}",
    f"THREADS: {THREADS}",
    f"Input directory: {INPUT_DIR}",
    f"FASTA files selected: {len(INPUT_FASTAS)}",
    f"Total execution time (s): {TOTAL_SECONDS:.2f}",
    "=" * 72,
]


for number, item in enumerate(
    RUN_SUMMARY,
    start=1,
):

    report.extend(
        [
            "",
            f"FILE {number}",
            f"Input FASTA: {item['input']}",
            f"Output taxonomy: {item['output']}",
            f"Status: {item['status']}",
            f"Representative sequences: {item['sequences']}",
            (
                "SINTAX records: "
                f"{item['validation']['records']}"
            ),
            (
                "Taxonomy predictions: "
                f"{item['validation']['predictions']}"
            ),
            (
                "Cutoff-filtered predictions: "
                f"{item['validation']['cutoff_predictions']}"
            ),
            f"Execution time (s): {item['seconds']:.2f}",
        ]
    )

    if item["command"] is not None:

        report.extend(
            [
                "Command:",
                pretty_command(
                    item["command"]
                ),
            ]
        )


report.extend(
    [
        "",
        "=" * 72,
        "MetaDiv_sintax completed successfully",
        "=" * 72,
    ]
)


REPORT_FILE.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8",
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 72)
print("MetaDiv_sintax completed successfully")
print("=" * 72)

for item in RUN_SUMMARY:

    print(
        f"{item['input'].name}"
        f"  ->  "
        f"{item['output'].name}"
        f"  [{item['status']}]"
    )

print(
    f"\nTotal execution time: "
    f"{TOTAL_SECONDS:.2f} s"
)

print(
    f"\nProcessing report:\n"
    f"{REPORT_FILE}"
)

print("=" * 72)
