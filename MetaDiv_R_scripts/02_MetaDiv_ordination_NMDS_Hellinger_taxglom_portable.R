# =============================================================================
# 02_MetaDiv_ordination_NMDS_Hellinger_taxglom_portable.R
# =============================================================================
# MetaDiv Explorer
#
# Portable ordination module using tax_glom, Hellinger transformation and NMDS.
#
# This module assumes that the phyloseq object was already created by the
# main MetaDiv R workflow and saved as:
#   ps_metadiv_builder.rds
#
# Workflow:
# 1. Resolve the MetaDiv For_R run directory
# 2. Load ps_metadiv_builder.rds
# 3. Collapse taxonomy with tax_glom
# 4. Apply Hellinger transformation
# 5. Run NMDS with Bray-Curtis distance
# 6. Join NMDS scores with metadata without duplicating SampleID
# 7. Export reusable objects, tables, summaries and plots
#
# Author: Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================


# =============================================================================
# 1. PORTABLE STARTUP
# =============================================================================

get_script_dir <- function() {

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = TRUE
    )))
  }

  source_file <- tryCatch(
    sys.frames()[[1]]$ofile,
    error = function(e) NULL
  )

  if (!is.null(source_file)) {
    return(dirname(normalizePath(
      source_file,
      winslash = "/",
      mustWork = TRUE
    )))
  }

  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}

SCRIPT_DIR <- get_script_dir()

CHECK_ENV_FILE <- file.path(SCRIPT_DIR, "check_environment.R")

if (!file.exists(CHECK_ENV_FILE)) {
  stop(
    "The environment checker was not found:\n",
    CHECK_ENV_FILE,
    "\nPlace check_environment.R in the same directory as this script."
  )
}

source(CHECK_ENV_FILE)
check_metadiv_environment()


is_wsl <- function() {

  if (.Platform$OS.type == "windows") {
    return(FALSE)
  }

  version_text <- tryCatch(
    paste(readLines("/proc/version", warn = FALSE), collapse = " "),
    error = function(e) ""
  )

  grepl("microsoft|wsl", version_text, ignore.case = TRUE)
}


windows_path_to_wsl <- function(path) {

  is_windows_path <- grepl("^[A-Za-z]:[/\\\\]", path)

  if (!is_windows_path) {
    return(path)
  }

  if (.Platform$OS.type == "windows") {
    return(path)
  }

  if (!is_wsl()) {
    stop(
      "A Windows path was provided, but this is not Windows or WSL:\n",
      path,
      "\nUse a native Linux/macOS path on this system."
    )
  }

  drive <- tolower(substr(path, 1, 1))
  remainder <- substring(path, 3)
  remainder <- gsub("\\\\", "/", remainder)

  paste0("/mnt/", drive, remainder)
}


resolve_for_r_dir <- function() {

  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) >= 1 && nzchar(args[1])) {

    selected_dir <- args[1]

  } else if (
    exists("FOR_R_DIR", envir = .GlobalEnv) &&
    nzchar(get("FOR_R_DIR", envir = .GlobalEnv))
  ) {

    selected_dir <- get("FOR_R_DIR", envir = .GlobalEnv)

  } else if (.Platform$OS.type == "windows" && interactive()) {

    selected_dir <- choose.dir(
      caption = "Select the MetaDiv For_R dataset directory"
    )

    if (is.na(selected_dir) || !nzchar(selected_dir)) {
      stop("No dataset directory was selected.")
    }

  } else if (interactive()) {

    selected_dir <- readline(
      "Enter the MetaDiv For_R dataset directory: "
    )

  } else {

    stop(
      "Dataset directory not provided.\n",
      "Usage:\n",
      "Rscript 02_MetaDiv_ordination_NMDS_Hellinger_taxglom_portable.R ",
      "/path/to/For_R/run"
    )
  }

  selected_dir <- windows_path_to_wsl(selected_dir)

  normalizePath(
    selected_dir,
    winslash = "/",
    mustWork = TRUE
  )
}


# =============================================================================
# 2. CONFIGURATION
# =============================================================================

FOR_R_DIR <- resolve_for_r_dir()

# All R modules use the same output directory created by Module 01.
BASE_OUTPUT_DIR <- file.path(FOR_R_DIR, "R_analysis_output")
EXPLORER_DIR <- BASE_OUTPUT_DIR

PS_FILE <- file.path(EXPLORER_DIR, "ps_metadiv_builder.rds")
NMDS_DIR <- file.path(EXPLORER_DIR, "02_Ordination_NMDS_Hellinger")

dir.create(
  NMDS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Analysis parameters
tax_rank <- "genus"
distance_method <- "bray"
transformation_method <- "hellinger"
trymax_value <- 300
seed_value <- 123
FORCE_RERUN_NMDS <- FALSE

# The first existing and informative variable will be used for coloring.
# Add project-specific metadata fields here when needed.
candidate_color_variables <- c(
  "vegetacion_CONABIO",
  "ecoregion_WWF",
  "wwf_ecoregions_REALM",
  "Ecosystem",
  "sequencing_platform",
  "sequencing_technology",
  "year_atlas",
  "Group",
  "Site",
  "biome",
  "country"
)

cat("\nMetaDiv NMDS configuration loaded.\n")
cat("Input For_R directory:\n", FOR_R_DIR, "\n")
cat("Explorer directory:\n", EXPLORER_DIR, "\n")
cat("Phyloseq object:\n", PS_FILE, "\n")
cat("NMDS output directory:\n", NMDS_DIR, "\n")


# =============================================================================
# 3. LOAD PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tibble)
})


# =============================================================================
# 4. LOAD AND VALIDATE PHYLOSEQ OBJECT
# =============================================================================

if (!file.exists(PS_FILE)) {
  stop(
    "Phyloseq object not found:\n",
    PS_FILE,
    "\nRun the main MetaDiv R script first for this dataset."
  )
}

ps <- readRDS(PS_FILE)

if (!inherits(ps, "phyloseq")) {
  stop("The loaded RDS object is not a phyloseq object: ", PS_FILE)
}

ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

if (ntaxa(ps) == 0) {
  stop("No taxa with positive abundance remain after pruning.")
}

if (nsamples(ps) < 2) {
  stop("NMDS requires at least two non-empty samples.")
}

cat("\nPhyloseq object loaded.\n")
cat("Taxa:", ntaxa(ps), "\n")
cat("Samples:", nsamples(ps), "\n")
cat("Total reads:", sum(sample_sums(ps)), "\n")


# =============================================================================
# 5. VALIDATE TAXONOMIC RANK
# =============================================================================

available_ranks <- rank_names(ps)

cat("\nAvailable taxonomic ranks:\n")
print(available_ranks)

if (!tax_rank %in% available_ranks) {
  stop(
    "The taxonomic rank '", tax_rank,
    "' does not exist in tax_table(ps).\nAvailable ranks: ",
    paste(available_ranks, collapse = ", ")
  )
}


# =============================================================================
# 6. DEFINE OUTPUT FILES
# =============================================================================

base_suffix <- paste0(
  "taxglom_",
  tax_rank,
  "_",
  transformation_method,
  "_",
  distance_method
)

PS_HELLINGER_FILE <- file.path(
  NMDS_DIR,
  paste0("ps_", base_suffix, ".rds")
)

NMDS_OBJECT_FILE <- file.path(
  NMDS_DIR,
  paste0("nmds_object_", base_suffix, ".rds")
)

NMDS_TABLE_FILE <- file.path(
  NMDS_DIR,
  paste0("nmds_scores_metadata_", base_suffix, ".rds")
)

NMDS_CSV_FILE <- file.path(
  NMDS_DIR,
  paste0("NMDS_scores_metadata_", base_suffix, ".csv")
)


# =============================================================================
# 7. HELPER FUNCTIONS
# =============================================================================

add_rownames_safe <- function(df, colname = "SampleID") {

  df <- as.data.frame(df, check.names = FALSE)

  if (colname %in% names(df)) {

    replacement <- paste0(colname, "_original")

    while (replacement %in% names(df)) {
      replacement <- paste0(replacement, "_copy")
    }

    names(df)[names(df) == colname] <- replacement

    message(
      "Existing metadata column '", colname,
      "' was renamed to '", replacement,
      "' to avoid duplication."
    )
  }

  tibble::rownames_to_column(df, var = colname)
}


select_color_variable <- function(df, candidates) {

  available <- candidates[candidates %in% names(df)]

  if (length(available) == 0) {
    return(NULL)
  }

  informative <- available[vapply(
    df[available],
    function(x) {
      values <- unique(x[!is.na(x) & trimws(as.character(x)) != ""])
      length(values) > 1
    },
    logical(1)
  )]

  if (length(informative) == 0) {
    return(NULL)
  }

  informative[1]
}


# =============================================================================
# 8. HEAVY STEP: TAX_GLOM + HELLINGER + NMDS
# =============================================================================

cached_files_exist <- all(file.exists(c(
  PS_HELLINGER_FILE,
  NMDS_OBJECT_FILE,
  NMDS_TABLE_FILE
)))

if (cached_files_exist && !FORCE_RERUN_NMDS) {

  message("Existing NMDS objects found. Loading cached results.")

  ps_hellinger <- readRDS(PS_HELLINGER_FILE)
  nmds <- readRDS(NMDS_OBJECT_FILE)
  nmds_df <- readRDS(NMDS_TABLE_FILE)

} else {

  message("Running NMDS heavy step.")

  ps_glom <- tax_glom(
    ps,
    taxrank = tax_rank,
    NArm = FALSE
  )

  ps_glom <- prune_taxa(taxa_sums(ps_glom) > 0, ps_glom)
  ps_glom <- prune_samples(sample_sums(ps_glom) > 0, ps_glom)

  if (ntaxa(ps_glom) < 2) {
    stop(
      "Fewer than two taxa remain after tax_glom at rank '",
      tax_rank,
      "'. NMDS cannot be calculated reliably."
    )
  }

  cat("\nAfter tax_glom at", tax_rank, "level:\n")
  cat("Taxa:", ntaxa(ps_glom), "\n")
  cat("Samples:", nsamples(ps_glom), "\n")
  cat("Total reads:", sum(sample_sums(ps_glom)), "\n")

  ps_hellinger <- transform_sample_counts(
    ps_glom,
    function(x) {
      total <- sum(x)
      if (total == 0) return(x)
      sqrt(x / total)
    }
  )

  set.seed(seed_value)

  nmds <- ordinate(
    ps_hellinger,
    method = "NMDS",
    distance = distance_method,
    trymax = trymax_value,
    trace = FALSE
  )

  stress_value <- round(nmds$stress, 4)
  cat("\nNMDS stress value:", stress_value, "\n")

  nmds_scores <- scores(nmds, display = "sites")
  nmds_scores <- add_rownames_safe(nmds_scores, "SampleID")

  if (!is.null(sample_data(ps_hellinger, errorIfNULL = FALSE))) {

    metadata_df <- as(sample_data(ps_hellinger), "data.frame")
    metadata_df <- add_rownames_safe(metadata_df, "SampleID")

    nmds_df <- left_join(
      nmds_scores,
      metadata_df,
      by = "SampleID"
    )

  } else {

    message("No sample metadata were found. NMDS scores will be exported alone.")
    nmds_df <- nmds_scores
  }

  nmds_df$distance_from_center <- sqrt(
    nmds_df$NMDS1^2 + nmds_df$NMDS2^2
  )

  saveRDS(ps_hellinger, PS_HELLINGER_FILE)
  saveRDS(nmds, NMDS_OBJECT_FILE)
  saveRDS(nmds_df, NMDS_TABLE_FILE)

  write.csv(
    nmds_df,
    NMDS_CSV_FILE,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("\nNMDS heavy step completed.\n")
}


# =============================================================================
# 9. SELECT COLOR VARIABLE
# =============================================================================

color_var <- select_color_variable(
  nmds_df,
  candidate_color_variables
)

cat("\nAvailable NMDS/metadata columns:\n")
print(names(nmds_df))

if (is.null(color_var)) {
  message("No informative preferred metadata variable was found. Plotting in gray.")
} else {
  message("Using '", color_var, "' as the NMDS color variable.")
}


# =============================================================================
# 10. GENERATE NMDS PLOT WITHOUT DUPLICATED POINT LAYERS
# =============================================================================

stress_value <- round(nmds$stress, 4)

if (!is.null(color_var)) {

  color_values <- nmds_df[[color_var]]

  p <- ggplot(
    nmds_df,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = .data[[color_var]]
    )
  ) +
    geom_point(size = 3, alpha = 0.8)

  if (is.numeric(color_values)) {

    p <- p +
      scale_color_viridis_c(
        option = "plasma",
        na.value = "gray50"
      )

  } else {

    nmds_df[[color_var]] <- as.factor(nmds_df[[color_var]])

    p <- ggplot(
      nmds_df,
      aes(
        x = NMDS1,
        y = NMDS2,
        color = .data[[color_var]]
      )
    ) +
      geom_point(size = 3, alpha = 0.8) +
      scale_color_viridis_d(
        option = "plasma",
        na.value = "gray50"
      )
  }

  plot_suffix <- paste0(
    base_suffix,
    "_colored_by_",
    color_var
  )

} else {

  p <- ggplot(
    nmds_df,
    aes(
      x = NMDS1,
      y = NMDS2
    )
  ) +
    geom_point(
      size = 3,
      alpha = 0.8,
      color = "gray40"
    )

  plot_suffix <- paste0(
    base_suffix,
    "_gray_no_metadata"
  )
}

p <- p +
  labs(
    title = paste0("NMDS (Bray-Curtis) | tax_glom: ", tax_rank),
    subtitle = paste0(
      "Transformation: Hellinger | Stress = ",
      stress_value,
      if (!is.null(color_var)) paste0(" | Color: ", color_var) else ""
    ),
    x = "NMDS1",
    y = "NMDS2",
    color = color_var
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = if (is.null(color_var)) "none" else "right"
  )

print(p)


# =============================================================================
# 11. EXPORT PLOT, OUTLIERS AND SUMMARY
# =============================================================================

outliers_top10 <- nmds_df %>%
  arrange(desc(distance_from_center)) %>%
  head(10)

write.csv(
  outliers_top10,
  file.path(
    NMDS_DIR,
    paste0("NMDS_top10_outliers_", plot_suffix, ".csv")
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

ggsave(
  filename = file.path(
    NMDS_DIR,
    paste0("NMDS_", plot_suffix, ".png")
  ),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(
    NMDS_DIR,
    paste0("NMDS_", plot_suffix, ".pdf")
  ),
  plot = p,
  width = 10,
  height = 7
)

nmds_summary <- data.frame(
  Analysis = "NMDS",
  Taxonomic_Collapse = tax_rank,
  Transformation = transformation_method,
  Distance = distance_method,
  Stress = stress_value,
  Taxa = ntaxa(ps_hellinger),
  Samples = nsamples(ps_hellinger),
  Color_Variable = if (is.null(color_var)) "None" else color_var,
  Color_Variable_Found = !is.null(color_var),
  Cached_Result_Used = cached_files_exist && !FORCE_RERUN_NMDS,
  Seed = seed_value,
  Trymax = trymax_value,
  stringsAsFactors = FALSE
)

write.csv(
  nmds_summary,
  file.path(
    NMDS_DIR,
    paste0("NMDS_summary_", base_suffix, ".csv")
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 12. SESSION INFO AND FINAL MESSAGE
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    NMDS_DIR,
    "R_sessionInfo_module02_NMDS_Hellinger.txt"
  )
)

cat("\nNMDS ordination module completed successfully.\n")
cat("Outputs saved in:\n", NMDS_DIR, "\n")
cat("Stress value:", stress_value, "\n")
cat(
  "Color variable used:",
  if (is.null(color_var)) "None" else color_var,
  "\n"
)
