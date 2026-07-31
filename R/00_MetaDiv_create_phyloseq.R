# =============================================================================
# Portable startup shared by MetaDiv R modules
# =============================================================================

SCRIPT_DIR <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = TRUE
    ))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    dirname(normalizePath(
      sys.frames()[[1]]$ofile,
      winslash = "/",
      mustWork = TRUE
    ))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
})

CHECK_ENV_FILE <- file.path(SCRIPT_DIR, "check_environment.R")

if (!file.exists(CHECK_ENV_FILE)) {
  stop(
    "check_environment.R was not found in the same directory as this script:\n",
    CHECK_ENV_FILE,
    call. = FALSE
  )
}

source(CHECK_ENV_FILE)
check_metadiv_environment()

windows_path_to_wsl <- function(path) {
  path <- trimws(path)

  if (.Platform$OS.type == "windows") {
    return(path)
  }

  path <- gsub("\\\\", "/", path)

  if (grepl("^[A-Za-z]:/", path)) {
    drive_letter <- tolower(substr(path, 1, 1))
    path <- paste0("/mnt/", drive_letter, substring(path, 3))
    message("[MetaDiv] Windows path converted for WSL: ", path)
  }

  path
}

resolve_for_r_dir <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) >= 1 && nzchar(args[1])) {
    selected_dir <- args[1]
  } else if (
    exists("FOR_R_DIR", inherits = TRUE) &&
    is.character(FOR_R_DIR) &&
    length(FOR_R_DIR) == 1 &&
    nzchar(FOR_R_DIR)
  ) {
    selected_dir <- FOR_R_DIR
  } else if (interactive() && .Platform$OS.type == "windows") {
    selected_dir <- utils::choose.dir(
      caption = "Select the MetaDiv For_R analysis directory"
    )

    if (is.na(selected_dir) || !nzchar(selected_dir)) {
      stop("No input directory was selected.", call. = FALSE)
    }
  } else if (interactive()) {
    selected_dir <- readline(
      prompt = "Enter the MetaDiv For_R directory: "
    )
  } else {
    stop(
      paste(
        "No input directory was provided.",
        "Run the script as:",
        "Rscript script.R /path/to/For_R_directory"
      ),
      call. = FALSE
    )
  }

  selected_dir <- windows_path_to_wsl(selected_dir)

  normalizePath(
    selected_dir,
    winslash = "/",
    mustWork = TRUE
  )
}

progress_message <- function(step, text) {
  cat(
    sprintf(
      "\n[%s] %s | %s\n",
      step,
      text,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  )
  flush.console()
}

# =============================================================================
# 00_MetaDiv_create_phyloseq.R
# =============================================================================
# Creates the MetaDiv phyloseq object.
#
# Sample-selection rule:
#   1. If sample_metadata.csv exists and contains valid SampleID values,
#      metadata is the master sample list.
#   2. The metadata may contain only one column: SampleID.
#   3. Only SampleID values shared by metadata and abundance_table.csv are kept.
#   4. If metadata is absent, empty, or has no valid SampleID values, all
#      abundance-table samples are retained and minimal metadata is generated.
#
# Author: Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================

FOR_R_DIR <- resolve_for_r_dir()

OUTPUT_DIR <- file.path(FOR_R_DIR, "R_analysis_output")
PS_FILE <- file.path(OUTPUT_DIR, "ps_metadiv_builder.rds")

# Set TRUE to overwrite an existing phyloseq object.
FORCE_REBUILD_PS <- FALSE

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(data.table)
})

abundance_path <- file.path(FOR_R_DIR, "abundance_table.csv")
taxonomy_path  <- file.path(FOR_R_DIR, "taxonomy_table.csv")
metadata_path  <- file.path(FOR_R_DIR, "sample_metadata.csv")
fasta_path     <- file.path(FOR_R_DIR, "sequences.fasta")

required_files <- c(abundance_path, taxonomy_path)
missing_required <- required_files[!file.exists(required_files)]

if (length(missing_required) > 0) {
  stop(
    "Required input files were not found:\n",
    paste0(" - ", missing_required, collapse = "\n"),
    call. = FALSE
  )
}

make_minimal_metadata <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))

  metadata <- data.frame(
    SampleID = sample_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  rownames(metadata) <- metadata$SampleID
  metadata
}

load_metadata_as_master <- function(metadata_path, abundance_sample_ids) {
  abundance_sample_ids <- trimws(as.character(abundance_sample_ids))

  if (!file.exists(metadata_path)) {
    message(
      "[MetaDiv] sample_metadata.csv was not found. ",
      "All abundance-table samples will be retained."
    )

    return(list(
      metadata = make_minimal_metadata(abundance_sample_ids),
      metadata_controls_samples = FALSE,
      requested_sample_ids = abundance_sample_ids
    ))
  }

  metadata <- fread(
    metadata_path,
    data.table = FALSE,
    check.names = FALSE
  )

  if (nrow(metadata) == 0 || ncol(metadata) == 0) {
    message(
      "[MetaDiv] sample_metadata.csv is empty. ",
      "All abundance-table samples will be retained."
    )

    return(list(
      metadata = make_minimal_metadata(abundance_sample_ids),
      metadata_controls_samples = FALSE,
      requested_sample_ids = abundance_sample_ids
    ))
  }

  if (!"SampleID" %in% colnames(metadata)) {
    message(
      "[MetaDiv] SampleID was not found. ",
      "The first metadata column will be used as SampleID."
    )
    colnames(metadata)[1] <- "SampleID"
  }

  metadata$SampleID <- trimws(as.character(metadata$SampleID))

  metadata <- metadata[
    !is.na(metadata$SampleID) & metadata$SampleID != "",
    ,
    drop = FALSE
  ]

  duplicated_ids <- unique(
    metadata$SampleID[duplicated(metadata$SampleID)]
  )

  if (length(duplicated_ids) > 0) {
    write.csv(
      data.frame(SampleID = duplicated_ids),
      file.path(OUTPUT_DIR, "metadata_duplicated_SampleID_removed.csv"),
      row.names = FALSE
    )

    message(
      "[MetaDiv] Duplicate SampleID values were found. ",
      "Only the first occurrence of each SampleID will be used."
    )

    metadata <- metadata[
      !duplicated(metadata$SampleID),
      ,
      drop = FALSE
    ]
  }

  if (nrow(metadata) == 0) {
    message(
      "[MetaDiv] No valid SampleID values remained after cleaning. ",
      "All abundance-table samples will be retained."
    )

    return(list(
      metadata = make_minimal_metadata(abundance_sample_ids),
      metadata_controls_samples = FALSE,
      requested_sample_ids = abundance_sample_ids
    ))
  }

  rownames(metadata) <- metadata$SampleID

  list(
    metadata = metadata,
    metadata_controls_samples = TRUE,
    requested_sample_ids = metadata$SampleID
  )
}

clean_taxonomy_table <- function(tax_df) {
  tax_df[] <- lapply(tax_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })

  tax_df
}

load_sequences_safely <- function(fasta_path, taxa_ids) {
  if (!file.exists(fasta_path)) {
    message(
      "[MetaDiv] sequences.fasta was not found. ",
      "The phyloseq object will be created without refseq."
    )
    return(NULL)
  }

  seqs <- readDNAStringSet(fasta_path)
  seqs <- seqs[names(seqs) %in% taxa_ids]

  if (length(seqs) == 0) {
    message(
      "[MetaDiv] No compatible FASTA IDs were found. ",
      "The phyloseq object will be created without refseq."
    )
    return(NULL)
  }

  seqs
}

if (file.exists(PS_FILE) && !FORCE_REBUILD_PS) {
  message(
    "[MetaDiv] Existing phyloseq object found. ",
    "Set FORCE_REBUILD_PS <- TRUE to rebuild it."
  )

  ps <- readRDS(PS_FILE)

  cat("\nExisting phyloseq object:\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")

} else {
  progress_message("00.1", "Loading abundance_table.csv")

  otu_dt <- fread(
    abundance_path,
    data.table = FALSE,
    check.names = FALSE
  )

  if (ncol(otu_dt) < 2) {
    stop(
      "abundance_table.csv must contain an ID column and at least one sample column.",
      call. = FALSE
    )
  }

  otu_ids <- trimws(as.character(otu_dt[[1]]))

  if (anyNA(otu_ids) || any(otu_ids == "") || anyDuplicated(otu_ids)) {
    stop(
      "The first abundance-table column must contain unique, non-empty taxon IDs.",
      call. = FALSE
    )
  }

  otu_df <- otu_dt[, -1, drop = FALSE]
  rownames(otu_df) <- otu_ids

  rm(otu_dt)
  gc()

  otu_mat <- as.matrix(otu_df)
  rm(otu_df)
  gc()

  storage.mode(otu_mat) <- "numeric"
  otu_mat[is.na(otu_mat)] <- 0

  abundance_sample_ids <- trimws(as.character(colnames(otu_mat)))
  colnames(otu_mat) <- abundance_sample_ids

  if (any(abundance_sample_ids == "") || anyDuplicated(abundance_sample_ids)) {
    stop(
      "Abundance-table sample names must be unique and non-empty.",
      call. = FALSE
    )
  }

  progress_message("00.2", "Loading sample_metadata.csv")

  metadata_result <- load_metadata_as_master(
    metadata_path,
    abundance_sample_ids
  )

  metadata_current <- metadata_result$metadata
  metadata_controls_samples <- metadata_result$metadata_controls_samples
  requested_sample_ids <- metadata_result$requested_sample_ids

  included_samples <- requested_sample_ids[
    requested_sample_ids %in% abundance_sample_ids
  ]

  excluded_abundance_samples <- setdiff(
    abundance_sample_ids,
    included_samples
  )

  metadata_samples_not_in_abundance <- setdiff(
    requested_sample_ids,
    abundance_sample_ids
  )

  if (length(included_samples) == 0) {
    stop(
      "No SampleID values are shared between sample_metadata.csv and abundance_table.csv.",
      call. = FALSE
    )
  }

  otu_mat <- otu_mat[, included_samples, drop = FALSE]
  metadata_current <- metadata_current[included_samples, , drop = FALSE]

  write.csv(
    data.frame(SampleID = abundance_sample_ids),
    file.path(OUTPUT_DIR, "samples_in_abundance.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(SampleID = requested_sample_ids),
    file.path(OUTPUT_DIR, "samples_requested_by_metadata.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(SampleID = included_samples),
    file.path(OUTPUT_DIR, "samples_included_phyloseq.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(SampleID = excluded_abundance_samples),
    file.path(OUTPUT_DIR, "samples_excluded_not_in_metadata.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(SampleID = metadata_samples_not_in_abundance),
    file.path(OUTPUT_DIR, "metadata_samples_not_in_abundance.csv"),
    row.names = FALSE
  )

  cat("\nSample-selection diagnostics:\n")
  cat("Samples in abundance table:", length(abundance_sample_ids), "\n")
  cat("Valid SampleID values in metadata:", length(requested_sample_ids), "\n")
  cat("Metadata controls sample inclusion:", metadata_controls_samples, "\n")
  cat("Shared samples included:", length(included_samples), "\n")
  cat(
    "Abundance samples excluded because they are absent from metadata:",
    length(excluded_abundance_samples),
    "\n"
  )
  cat(
    "Metadata samples absent from abundance table:",
    length(metadata_samples_not_in_abundance),
    "\n"
  )

  progress_message("00.3", "Loading taxonomy_table.csv")

  tax_dt <- fread(
    taxonomy_path,
    data.table = FALSE,
    check.names = FALSE
  )

  if (ncol(tax_dt) < 2) {
    stop(
      "taxonomy_table.csv must contain an ID column and at least one taxonomy column.",
      call. = FALSE
    )
  }

  tax_ids <- trimws(as.character(tax_dt[[1]]))
  tax_df <- tax_dt[, -1, drop = FALSE]
  rownames(tax_df) <- tax_ids

  rm(tax_dt)
  gc()

  tax_df <- clean_taxonomy_table(tax_df)

  common_taxa <- intersect(
    rownames(otu_mat),
    rownames(tax_df)
  )

  if (length(common_taxa) == 0) {
    stop(
      "No shared taxa IDs were found between abundance_table.csv and taxonomy_table.csv.",
      call. = FALSE
    )
  }

  taxa_missing_taxonomy <- setdiff(
    rownames(otu_mat),
    rownames(tax_df)
  )

  taxonomy_taxa_not_in_abundance <- setdiff(
    rownames(tax_df),
    rownames(otu_mat)
  )

  write.csv(
    data.frame(TaxonID = common_taxa),
    file.path(OUTPUT_DIR, "taxa_included_phyloseq.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(TaxonID = taxa_missing_taxonomy),
    file.path(OUTPUT_DIR, "abundance_taxa_missing_taxonomy.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(TaxonID = taxonomy_taxa_not_in_abundance),
    file.path(OUTPUT_DIR, "taxonomy_taxa_not_in_abundance.csv"),
    row.names = FALSE
  )

  otu_mat <- otu_mat[common_taxa, , drop = FALSE]
  tax_mat <- as.matrix(tax_df[common_taxa, , drop = FALSE])

  rm(tax_df)
  gc()

  progress_message("00.4", "Loading representative sequences")

  seqs <- load_sequences_safely(
    fasta_path,
    common_taxa
  )

  progress_message("00.5", "Creating phyloseq object")

  ps_components <- list(
    otu_table(otu_mat, taxa_are_rows = TRUE),
    tax_table(tax_mat),
    sample_data(metadata_current)
  )

  if (!is.null(seqs)) {
    ps_components <- c(ps_components, list(refseq(seqs)))
  }

  ps <- do.call(phyloseq, ps_components)

  taxa_before_pruning <- ntaxa(ps)
  samples_before_pruning <- nsamples(ps)

  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)

  saveRDS(ps, PS_FILE)

  metadata_used <- data.frame(
    sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"SampleID" %in% colnames(metadata_used)) {
    metadata_used$SampleID <- rownames(metadata_used)
  }

  sample_depth <- data.frame(
    SampleID = sample_names(ps),
    Reads = sample_sums(ps),
    stringsAsFactors = FALSE
  )

  write.csv(
    sample_depth,
    file.path(OUTPUT_DIR, "sample_read_depth.csv"),
    row.names = FALSE
  )

  write.csv(
    metadata_used,
    file.path(OUTPUT_DIR, "clean_metadata_used.csv"),
    row.names = FALSE
  )

  summary_lines <- c(
    "MetaDiv phyloseq summary",
    "========================",
    "",
    paste("Input directory:", FOR_R_DIR),
    paste("Output directory:", OUTPUT_DIR),
    paste("Saved phyloseq object:", PS_FILE),
    "",
    paste("Metadata controls sample inclusion:", metadata_controls_samples),
    paste("Samples in original abundance table:", length(abundance_sample_ids)),
    paste("Samples requested by metadata:", length(requested_sample_ids)),
    paste("Shared samples before zero-read pruning:", length(included_samples)),
    paste("Samples before zero-read pruning:", samples_before_pruning),
    paste("Samples after zero-read pruning:", nsamples(ps)),
    paste("Taxa before zero-abundance pruning:", taxa_before_pruning),
    paste("Taxa after zero-abundance pruning:", ntaxa(ps)),
    paste("Total reads:", sum(sample_sums(ps))),
    "",
    "Taxonomic ranks:",
    paste(rank_names(ps), collapse = ", "),
    "",
    "Sample variables:",
    paste(sample_variables(ps), collapse = ", ")
  )

  writeLines(
    summary_lines,
    file.path(OUTPUT_DIR, "ps_summary.txt")
  )

  capture.output(
    sessionInfo(),
    file = file.path(OUTPUT_DIR, "R_sessionInfo_module00_phyloseq.txt")
  )

  progress_message("00.6", "Phyloseq object created successfully")

  cat("\nFinal phyloseq object:\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")
  cat("Saved in:", PS_FILE, "\n")
}
