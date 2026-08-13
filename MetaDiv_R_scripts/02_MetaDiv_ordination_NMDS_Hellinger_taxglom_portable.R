# =============================================================================
# 02_MetaDiv_ordination_NMDS_Hellinger_taxglom_portable.R
# =============================================================================
# MetaDiv Explorer
#
# Portable ordination module using:
#   - tax_glom
#   - Hellinger transformation
#   - NMDS
#   - Bray-Curtis distance
#
# Default NMDS plots:
#   - sequencing_platform
#   - country
#   - Biome_detail
#
# Author: Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================


# =============================================================================
# 1. PORTABLE STARTUP
# =============================================================================

get_script_dir <- function() {
  
  args <- commandArgs(
    trailingOnly = FALSE
  )
  
  file_arg <- grep(
    "^--file=",
    args,
    value = TRUE
  )
  
  if (length(file_arg) > 0) {
    
    return(
      dirname(
        normalizePath(
          sub(
            "^--file=",
            "",
            file_arg[1]
          ),
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }
  
  
  source_file <- tryCatch(
    sys.frames()[[1]]$ofile,
    error = function(e) NULL
  )
  
  
  if (!is.null(source_file)) {
    
    return(
      dirname(
        normalizePath(
          source_file,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }
  
  
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}


SCRIPT_DIR <- get_script_dir()


CHECK_ENV_FILE <- file.path(
  SCRIPT_DIR,
  "check_environment.R"
)


if (!file.exists(CHECK_ENV_FILE)) {
  
  stop(
    "The environment checker was not found:\n",
    CHECK_ENV_FILE,
    "\nPlace check_environment.R in the same directory as this script."
  )
}


source(
  CHECK_ENV_FILE
)

check_metadiv_environment()


# =============================================================================
# WINDOWS / WSL PATH SUPPORT
# =============================================================================

is_wsl <- function() {
  
  if (.Platform$OS.type == "windows") {
    return(FALSE)
  }
  
  
  version_text <- tryCatch(
    
    paste(
      readLines(
        "/proc/version",
        warn = FALSE
      ),
      collapse = " "
    ),
    
    error = function(e) ""
  )
  
  
  grepl(
    "microsoft|wsl",
    version_text,
    ignore.case = TRUE
  )
}


windows_path_to_wsl <- function(path) {
  
  is_windows_path <- grepl(
    "^[A-Za-z]:[/\\\\]",
    path
  )
  
  
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
  
  
  drive <- tolower(
    substr(
      path,
      1,
      1
    )
  )
  
  
  remainder <- substring(
    path,
    3
  )
  
  
  remainder <- gsub(
    "\\\\",
    "/",
    remainder
  )
  
  
  paste0(
    "/mnt/",
    drive,
    remainder
  )
}


resolve_for_r_dir <- function() {
  
  args <- commandArgs(
    trailingOnly = TRUE
  )
  
  
  if (
    length(args) >= 1 &&
    nzchar(args[1])
  ) {
    
    selected_dir <- args[1]
    
    
  } else if (
    exists(
      "FOR_R_DIR",
      envir = .GlobalEnv
    ) &&
    nzchar(
      get(
        "FOR_R_DIR",
        envir = .GlobalEnv
      )
    )
  ) {
    
    selected_dir <- get(
      "FOR_R_DIR",
      envir = .GlobalEnv
    )
    
    
  } else if (
    .Platform$OS.type == "windows" &&
    interactive()
  ) {
    
    selected_dir <- choose.dir(
      caption = "Select the MetaDiv For_R dataset directory"
    )
    
    
    if (
      is.na(selected_dir) ||
      !nzchar(selected_dir)
    ) {
      
      stop(
        "No dataset directory was selected."
      )
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
  
  
  selected_dir <- windows_path_to_wsl(
    selected_dir
  )
  
  
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


BASE_OUTPUT_DIR <- file.path(
  FOR_R_DIR,
  "R_analysis_output"
)


EXPLORER_DIR <- BASE_OUTPUT_DIR


PS_FILE <- file.path(
  EXPLORER_DIR,
  "ps_metadiv_builder.rds"
)


NMDS_DIR <- file.path(
  EXPLORER_DIR,
  "02_Ordination_NMDS_Hellinger"
)


dir.create(
  NMDS_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# ANALYSIS PARAMETERS
# -----------------------------------------------------------------------------

tax_rank <- "genus"

distance_method <- "bray"

transformation_method <- "hellinger"

trymax_value <- 300

seed_value <- 123


# TRUE  = calculate NMDS again
# FALSE = reuse NMDS if it already exists
FORCE_RERUN_NMDS <- FALSE


# -----------------------------------------------------------------------------
# DEFAULT COLOR VARIABLES
# -----------------------------------------------------------------------------

default_color_variables <- c(
  "sequencing_platform",
  "country",
  "Biome_detail"
)


cat(
  "\nMetaDiv NMDS configuration loaded.\n"
)

cat(
  "Input For_R directory:\n",
  FOR_R_DIR,
  "\n"
)

cat(
  "Explorer directory:\n",
  EXPLORER_DIR,
  "\n"
)

cat(
  "Phyloseq object:\n",
  PS_FILE,
  "\n"
)

cat(
  "NMDS output directory:\n",
  NMDS_DIR,
  "\n"
)


# =============================================================================
# 3. LOAD PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  
  library(
    phyloseq
  )
  
  library(
    vegan
  )
  
  library(
    ggplot2
  )
  
  library(
    dplyr
  )
  
  library(
    tibble
  )
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


ps <- readRDS(
  PS_FILE
)


if (!inherits(
  ps,
  "phyloseq"
)) {
  
  stop(
    "The loaded RDS object is not a phyloseq object: ",
    PS_FILE
  )
}


ps <- prune_taxa(
  taxa_sums(ps) > 0,
  ps
)


ps <- prune_samples(
  sample_sums(ps) > 0,
  ps
)


if (ntaxa(ps) == 0) {
  
  stop(
    "No taxa with positive abundance remain after pruning."
  )
}


if (nsamples(ps) < 2) {
  
  stop(
    "NMDS requires at least two non-empty samples."
  )
}


cat(
  "\nPhyloseq object loaded.\n"
)

cat(
  "Taxa:",
  ntaxa(ps),
  "\n"
)

cat(
  "Samples:",
  nsamples(ps),
  "\n"
)

cat(
  "Total reads:",
  sum(
    sample_sums(ps)
  ),
  "\n"
)


# =============================================================================
# 5. VALIDATE TAXONOMIC RANK
# =============================================================================

available_ranks <- rank_names(
  ps
)


cat(
  "\nAvailable taxonomic ranks:\n"
)

print(
  available_ranks
)


if (!tax_rank %in% available_ranks) {
  
  stop(
    "The taxonomic rank '",
    tax_rank,
    "' does not exist in tax_table(ps).\n",
    "Available ranks: ",
    paste(
      available_ranks,
      collapse = ", "
    )
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
  paste0(
    "ps_",
    base_suffix,
    ".rds"
  )
)


NMDS_OBJECT_FILE <- file.path(
  NMDS_DIR,
  paste0(
    "nmds_object_",
    base_suffix,
    ".rds"
  )
)


NMDS_TABLE_FILE <- file.path(
  NMDS_DIR,
  paste0(
    "nmds_scores_metadata_",
    base_suffix,
    ".rds"
  )
)


NMDS_CSV_FILE <- file.path(
  NMDS_DIR,
  paste0(
    "NMDS_scores_metadata_",
    base_suffix,
    ".csv"
  )
)


# =============================================================================
# 7. HELPER FUNCTIONS
# =============================================================================

add_rownames_safe <- function(
    df,
    colname = "SampleID"
) {
  
  df <- as.data.frame(
    df,
    check.names = FALSE
  )
  
  
  if (colname %in% names(df)) {
    
    replacement <- paste0(
      colname,
      "_original"
    )
    
    
    while (
      replacement %in% names(df)
    ) {
      
      replacement <- paste0(
        replacement,
        "_copy"
      )
    }
    
    
    names(df)[
      names(df) == colname
    ] <- replacement
    
    
    message(
      "Existing metadata column '",
      colname,
      "' was renamed to '",
      replacement,
      "' to avoid duplication."
    )
  }
  
  
  tibble::rownames_to_column(
    df,
    var = colname
  )
}


# =============================================================================
# 8. HEAVY STEP: TAX_GLOM + HELLINGER + NMDS
# =============================================================================

cached_files_exist <- all(
  file.exists(
    c(
      PS_HELLINGER_FILE,
      NMDS_OBJECT_FILE
    )
  )
)


if (
  cached_files_exist &&
  !FORCE_RERUN_NMDS
) {
  
  message(
    "Existing NMDS objects found. Loading cached ordination."
  )
  
  
  ps_hellinger <- readRDS(
    PS_HELLINGER_FILE
  )
  
  
  nmds <- readRDS(
    NMDS_OBJECT_FILE
  )
  
  
} else {
  
  message(
    "Running NMDS heavy step."
  )
  
  
  # ---------------------------------------------------------------------------
  # TAXONOMIC AGGLOMERATION
  # ---------------------------------------------------------------------------
  
  ps_glom <- tax_glom(
    ps,
    taxrank = tax_rank,
    NArm = FALSE
  )
  
  
  ps_glom <- prune_taxa(
    taxa_sums(ps_glom) > 0,
    ps_glom
  )
  
  
  ps_glom <- prune_samples(
    sample_sums(ps_glom) > 0,
    ps_glom
  )
  
  
  if (ntaxa(ps_glom) < 2) {
    
    stop(
      "Fewer than two taxa remain after tax_glom at rank '",
      tax_rank,
      "'. NMDS cannot be calculated reliably."
    )
  }
  
  
  cat(
    "\nAfter tax_glom at",
    tax_rank,
    "level:\n"
  )
  
  cat(
    "Taxa:",
    ntaxa(ps_glom),
    "\n"
  )
  
  cat(
    "Samples:",
    nsamples(ps_glom),
    "\n"
  )
  
  cat(
    "Total reads:",
    sum(
      sample_sums(ps_glom)
    ),
    "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # HELLINGER TRANSFORMATION
  # ---------------------------------------------------------------------------
  
  ps_hellinger <- transform_sample_counts(
    ps_glom,
    function(x) {
      
      total <- sum(x)
      
      if (total == 0) {
        return(x)
      }
      
      sqrt(
        x / total
      )
    }
  )
  
  
  # ---------------------------------------------------------------------------
  # NMDS
  # ---------------------------------------------------------------------------
  
  set.seed(
    seed_value
  )
  
  
  nmds <- ordinate(
    ps_hellinger,
    method = "NMDS",
    distance = distance_method,
    trymax = trymax_value,
    trace = FALSE
  )
  
  
  saveRDS(
    ps_hellinger,
    PS_HELLINGER_FILE
  )
  
  
  saveRDS(
    nmds,
    NMDS_OBJECT_FILE
  )
  
  
  cat(
    "\nNMDS heavy step completed.\n"
  )
}


# =============================================================================
# 8B. BUILD NMDS TABLE USING CURRENT METADATA
# =============================================================================
#
# This part is intentionally reconstructed even when the NMDS itself
# comes from cache. This allows updated sample metadata to be used
# without recalculating the ordination.
# =============================================================================

stress_value <- round(
  nmds$stress,
  4
)


cat(
  "\nNMDS stress value:",
  stress_value,
  "\n"
)


nmds_scores <- scores(
  nmds,
  display = "sites"
)


nmds_scores <- add_rownames_safe(
  nmds_scores,
  "SampleID"
)


if (
  !is.null(
    sample_data(
      ps_hellinger,
      errorIfNULL = FALSE
    )
  )
) {
  
  metadata_df <- as(
    sample_data(ps_hellinger),
    "data.frame"
  )
  
  
  metadata_df <- add_rownames_safe(
    metadata_df,
    "SampleID"
  )
  
  
  nmds_df <- left_join(
    nmds_scores,
    metadata_df,
    by = "SampleID"
  )
  
  
} else {
  
  message(
    "No sample metadata were found. ",
    "NMDS scores will be exported alone."
  )
  
  
  nmds_df <- nmds_scores
}


nmds_df$distance_from_center <- sqrt(
  nmds_df$NMDS1^2 +
    nmds_df$NMDS2^2
)


saveRDS(
  nmds_df,
  NMDS_TABLE_FILE
)


write.csv(
  nmds_df,
  NMDS_CSV_FILE,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# =============================================================================
# 9. SELECT DEFAULT COLOR VARIABLES
# =============================================================================

available_color_variables <- default_color_variables[
  default_color_variables %in% names(nmds_df)
]


informative_color_variables <- available_color_variables[
  vapply(
    nmds_df[
      available_color_variables
    ],
    function(x) {
      
      values <- unique(
        x[
          !is.na(x) &
            trimws(
              as.character(x)
            ) != ""
        ]
      )
      
      
      length(values) > 1
    },
    logical(1)
  )
]


cat(
  "\nAvailable NMDS/metadata columns:\n"
)

print(
  names(nmds_df)
)


cat(
  "\nDefault NMDS color variables requested:\n"
)

print(
  default_color_variables
)


cat(
  "\nDefault variables found and informative:\n"
)

print(
  informative_color_variables
)


missing_color_variables <- setdiff(
  default_color_variables,
  available_color_variables
)


if (
  length(
    missing_color_variables
  ) > 0
) {
  
  message(
    "The following default metadata variables were not found: ",
    paste(
      missing_color_variables,
      collapse = ", "
    )
  )
}


noninformative_color_variables <- setdiff(
  available_color_variables,
  informative_color_variables
)


if (
  length(
    noninformative_color_variables
  ) > 0
) {
  
  message(
    "The following variables were found but contain fewer than ",
    "two informative values: ",
    paste(
      noninformative_color_variables,
      collapse = ", "
    )
  )
}


if (
  length(
    informative_color_variables
  ) == 0
) {
  
  stop(
    "None of the default NMDS color variables were found ",
    "with more than one informative value."
  )
}


# =============================================================================
# 10. GENERATE ONE NMDS PLOT FOR EACH DEFAULT VARIABLE
# =============================================================================

generated_plots <- list()

for (color_var in informative_color_variables) {
  
  message(
    "Generating NMDS colored by: ",
    color_var
  )
  
  plot_df <- nmds_df
  
  color_values <- plot_df[[color_var]]
  
  # ---------------------------------------------------------------------------
  # NUMERIC VARIABLE
  # ---------------------------------------------------------------------------
  
  if (is.numeric(color_values)) {
    
    p <- ggplot(
      plot_df,
      aes(
        x = NMDS1,
        y = NMDS2,
        color = .data[[color_var]]
      )
    ) +
      
      geom_point(
        size = 3,
        alpha = 0.8
      ) +
      
      scale_color_viridis_c(
        option = "plasma",
        na.value = "gray50"
      )
    
    # ---------------------------------------------------------------------------
    # CATEGORICAL VARIABLE
    # ---------------------------------------------------------------------------
    
  } else {
    
    plot_df[[color_var]] <- as.factor(
      plot_df[[color_var]]
    )
    
    p <- ggplot(
      plot_df,
      aes(
        x = NMDS1,
        y = NMDS2,
        color = .data[[color_var]]
      )
    ) +
      
      geom_point(
        size = 3,
        alpha = 0.8
      ) +
      
      scale_color_viridis_d(
        option = "plasma",
        na.value = "gray50"
      )
  }
  
  # ---------------------------------------------------------------------------
  # TITLES AND THEME
  # ---------------------------------------------------------------------------
  
  p <- p +
    
    labs(
      title = paste0(
        "NMDS (Bray-Curtis) | tax_glom: ",
        tax_rank
      ),
      
      subtitle = paste0(
        "Transformation: Hellinger",
        " | Stress = ",
        stress_value,
        " | Color: ",
        color_var
      ),
      
      x = "NMDS1",
      y = "NMDS2",
      color = color_var
    ) +
    
    theme_classic(
      base_size = 14
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right"
    )
  
  print(p)
  
  # ---------------------------------------------------------------------------
  # FILE NAME
  # ---------------------------------------------------------------------------
  
  plot_suffix <- paste0(
    base_suffix,
    "_colored_by_",
    color_var
  )
  
  # ---------------------------------------------------------------------------
  # EXPORT PNG
  # ---------------------------------------------------------------------------
  
  ggsave(
    filename = file.path(
      NMDS_DIR,
      paste0(
        "NMDS_",
        plot_suffix,
        ".png"
      )
    ),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  # ---------------------------------------------------------------------------
  # EXPORT PDF
  # ---------------------------------------------------------------------------
  
  ggsave(
    filename = file.path(
      NMDS_DIR,
      paste0(
        "NMDS_",
        plot_suffix,
        ".pdf"
      )
    ),
    plot = p,
    width = 10,
    height = 7
  )
  
  # ---------------------------------------------------------------------------
  # SAVE PLOT IN MEMORY
  # ---------------------------------------------------------------------------
  
  generated_plots[[color_var]] <- p
}