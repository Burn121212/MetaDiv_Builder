# =============================================================================
# 05_MetaDiv_ordination_NMDS_Hellinger_taxglom.R
# =============================================================================
# MetaDiv Explorer
#
# Ordination module using tax_glom, Hellinger transformation and NMDS.
#
# This module assumes that the phyloseq object was already created by:
# 00_01_02_03_MetaDiv_ps_explorer_alpha_resume.R
#
# Default workflow:
# 1. Load ps_metadiv_builder.rds
# 2. Collapse taxonomy using tax_glom
# 3. Apply Hellinger transformation
# 4. Run NMDS using Bray-Curtis distance
# 5. Export reusable NMDS objects, scores and plots
#
# Author: Bernardo Águila (modified to avoid duplicate SampleID column)
# Instituto de Biología, UNAM
# =============================================================================


# =============================================================================
# 1. CONFIGURATION
# =============================================================================
# change this line for each dataset (16S, ITS, CO1)
FOR_R_DIR <- "C:/Users/Montserrat/Downloads/MetaDiv Builder v1.7.3-20260626T191510Z-3-001/MetaDiv Builder v1.7.3/output/CO1/For_R/all_eukaryotes/species_only_p1_sppn08"

DATASET_OUTPUT_ROOT <- dirname(dirname(FOR_R_DIR))
RUN_NAME <- basename(FOR_R_DIR)
BASE_OUTPUT_DIR <- file.path(DATASET_OUTPUT_ROOT, "R_analysis_output")
EXPLORER_DIR <- file.path(BASE_OUTPUT_DIR, RUN_NAME)

PS_FILE <- file.path(EXPLORER_DIR, "ps_metadiv_builder.rds")
NMDS_DIR <- file.path(EXPLORER_DIR, "05_Ordination_NMDS_Hellinger")
dir.create(NMDS_DIR, recursive = TRUE, showWarnings = FALSE)

tax_rank <- "genus"
distance_method <- "bray"
transformation_method <- "hellinger"
trymax_value <- 300
seed_value <- 123
FORCE_RERUN_NMDS <- FALSE

# Preferred metadata columns for automatic coloring.
# The script will use the first available column in this list.
# If none are available, the NMDS will be plotted in gray.
candidate_color_variables <- c(
  "sequencing_platform",
  "Group",
  "Site",
  "biome",
  "country",
  "wwf_ecoregions_REALM"
)

cat("\nMetaDiv NMDS configuration loaded.\n")
cat("\nInput For_R directory:\n", FOR_R_DIR, "\n")
cat("\nExplorer directory:\n", EXPLORER_DIR, "\n")
cat("\nPhyloseq object:\n", PS_FILE, "\n")
cat("\nNMDS output directory:\n", NMDS_DIR, "\n")


# =============================================================================
# 2. LOAD PACKAGES
# =============================================================================
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)


# =============================================================================
# 3. LOAD PHYLOSEQ OBJECT
# =============================================================================
if (!file.exists(PS_FILE)) {
  stop(paste("Phyloseq object not found:", PS_FILE))
}

ps <- readRDS(PS_FILE)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

cat("\nPhyloseq object loaded.\n")
cat("Taxa:", ntaxa(ps), "\n")
cat("Samples:", nsamples(ps), "\n")
cat("Total reads:", sum(sample_sums(ps)), "\n")


# =============================================================================
# 4. VALIDATE TAXONOMIC RANK
# =============================================================================
available_ranks <- rank_names(ps)
cat("\nAvailable taxonomic ranks:\n")
print(available_ranks)

if (!tax_rank %in% available_ranks) {
  stop(paste("The taxonomic rank", tax_rank, "does not exist in tax_table(ps)."))
}


# =============================================================================
# 5. DEFINE OUTPUT SUFFIXES AND FILE PATHS
# =============================================================================
base_suffix <- paste0(
  "taxglom_",
  tax_rank,
  "_",
  transformation_method,
  "_",
  distance_method
)

PS_HELLINGER_FILE <- file.path(NMDS_DIR, paste0("ps_", base_suffix, ".rds"))
NMDS_OBJECT_FILE <- file.path(NMDS_DIR, paste0("nmds_object_", base_suffix, ".rds"))
NMDS_TABLE_FILE <- file.path(NMDS_DIR, paste0("nmds_scores_metadata_", base_suffix, ".rds"))
NMDS_CSV_FILE <- file.path(NMDS_DIR, paste0("NMDS_scores_metadata_", base_suffix, ".csv"))


# =============================================================================
# 6. HEAVY STEP: TAX_GLOM + HELLINGER + NMDS
# =============================================================================

# ---- Helper function to add rownames as a column without duplicates ----
add_rownames_safe <- function(df, colname = "SampleID") {
  if (colname %in% colnames(df)) {
    # Rename the existing column to avoid duplication
    colnames(df)[colnames(df) == colname] <- paste0(colname, "_original")
    message("NOTE: Existing column '", colname, "' was renamed to '", 
            colname, "_original' to avoid duplication.")
  }
  rownames_to_column(df, var = colname)
}

# ---- Main NMDS computation ----
if (
  file.exists(PS_HELLINGER_FILE) &&
  file.exists(NMDS_OBJECT_FILE) &&
  file.exists(NMDS_TABLE_FILE) &&
  !FORCE_RERUN_NMDS
) {
  
  message("Existing NMDS objects found. Loading previous results.")
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
  
  cat("\nAfter tax_glom at", tax_rank, "level:\n")
  cat("Taxa:", ntaxa(ps_glom), "\n")
  cat("Samples:", nsamples(ps_glom), "\n")
  cat("Total reads:", sum(sample_sums(ps_glom)), "\n")
  
  ps_hellinger <- transform_sample_counts(
    ps_glom,
    function(x) {
      if (sum(x) == 0) return(x)
      sqrt(x / sum(x))
    }
  )
  
  set.seed(seed_value)
  
  nmds <- ordinate(
    ps_hellinger,
    method = "NMDS",
    distance = distance_method,
    trymax = trymax_value
  )
  
  stress_value <- round(nmds$stress, 4)
  cat("\nNMDS stress value:", stress_value, "\n")
  
  # Extract NMDS scores (sites)
  nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
  # Add SampleID column from rownames, avoiding duplicates
  nmds_scores <- add_rownames_safe(nmds_scores, "SampleID")
  
  # Extract metadata
  metadata_df <- as(sample_data(ps_hellinger), "data.frame")
  # Add SampleID column from rownames, avoiding duplicates
  metadata_df <- add_rownames_safe(metadata_df, "SampleID")
  
  # Join scores with metadata
  nmds_df <- left_join(
    nmds_scores,
    metadata_df,
    by = "SampleID"
  )
  
  # Add distance from center
  nmds_df$distance_from_center <- sqrt(
    nmds_df$NMDS1^2 + nmds_df$NMDS2^2
  )
  
  # Save objects
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
# 7. SELECT COLOR VARIABLE FOR PLOTTING
# =============================================================================
# Find the first available column in the preferred list
color_var <- NULL
for (var in candidate_color_variables) {
  if (var %in% colnames(nmds_df)) {
    color_var <- var
    break
  }
}

if (is.null(color_var)) {
  message("No preferred color variable found. Samples will be plotted in gray.")
  color_var <- NULL
} else {
  message("Using '", color_var, "' as color variable for NMDS plot.")
}


# =============================================================================
# 8. GENERATE NMDS PLOTS
# =============================================================================

# Base plot
p <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = paste("NMDS -", tax_rank, "level"),
    subtitle = paste("Stress =", round(nmds$stress, 4)),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Add color if a variable was found
if (!is.null(color_var)) {
  p <- p +
    aes(color = .data[[color_var]]) +
    labs(color = color_var) +
    scale_color_viridis_d(option = "plasma", na.value = "gray50")
} else {
  p <- p + 
    geom_point(color = "gray40") +
    theme(legend.position = "none")
}

# Print plot to console
print(p)

# Save plots
ggsave(
  filename = file.path(NMDS_DIR, paste0("NMDS_plot_", base_suffix, ".png")),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(NMDS_DIR, paste0("NMDS_plot_", base_suffix, ".pdf")),
  plot = p,
  width = 8,
  height = 6
)

cat("\nNMDS plots saved in:\n", NMDS_DIR, "\n")

# Optional: also save a version with sample labels (if desired)
# p_labels <- p + geom_text(aes(label = SampleID), size = 2, vjust = -0.5)
# ggsave(..., file = "NMDS_plot_with_labels.png", ...)

cat("\n=== Script finished successfully ===\n")