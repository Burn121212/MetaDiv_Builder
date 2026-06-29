# =============================================================================
# 00_01_02_03_MetaDiv_ps_explorer_alpha_resume.R
# =============================================================================
# MetaDiv Explorer
#
# Combined script:
# 00 - Global configuration
# 01 - Create or load phyloseq object
# 02 - Quick exploratory analysis
# 03 - Alpha diversity analysis using vegan
#
# Author:
# Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================


# =============================================================================
# 00. GLOBAL CONFIGURATION
# =============================================================================
# change this line for each dataset 16S ITS CO1
FOR_R_DIR <- "C:/Users/berna/Desktop/MetaDiv_Builder_V1_7_3/output/ITS/For_R/all_eukaryotes/species_only_p1_sppn08"

BASE_OUTPUT_DIR <- "C:/Users/berna/Desktop/MetaDiv_Builder_V1_7_3/output/ITS/For_R/all_eukaryotes/species_only_p1_sppn08"

DATASET_OUTPUT_ROOT <- dirname(dirname(FOR_R_DIR))

RUN_NAME <- basename(FOR_R_DIR)

BASE_OUTPUT_DIR <- file.path(DATASET_OUTPUT_ROOT, "R_analysis_output")

EXPLORER_DIR <- file.path(BASE_OUTPUT_DIR, RUN_NAME)

PS_OUTPUT_NAME <- "ps_metadiv_builder.rds"

PS_FILE <- file.path(EXPLORER_DIR, PS_OUTPUT_NAME)

FORCE_REBUILD_PS <- FALSE
FORCE_RERUN_EXPLORER <- FALSE
FORCE_RERUN_ALPHA <- FALSE

dir.create(EXPLORER_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\nMetaDiv configuration loaded.\n")
cat("\nInput directory:\n", FOR_R_DIR, "\n")
cat("\nOutput directory:\n", EXPLORER_DIR, "\n")
cat("\nPhyloseq object:\n", PS_FILE, "\n")


# =============================================================================
# 00.1 LOAD PACKAGES
# =============================================================================

library(phyloseq)
library(Biostrings)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(readr)
library(vegan)


# =============================================================================
# 00.2 INPUT PATHS
# =============================================================================

abundance_path <- file.path(FOR_R_DIR, "abundance_table.csv")
taxonomy_path  <- file.path(FOR_R_DIR, "taxonomy_table.csv")
metadata_path  <- file.path(FOR_R_DIR, "sample_metadata.csv")
fasta_path     <- file.path(FOR_R_DIR, "sequences.fasta")


# =============================================================================
# 00.3 HELPER FUNCTIONS
# =============================================================================

make_minimal_metadata <- function(sample_ids) {
  
  sample_ids <- trimws(as.character(sample_ids))
  
  metadata <- data.frame(
    SampleID = sample_ids,
    Group = "All_samples",
    Site = sample_ids,
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- metadata$SampleID
  
  return(metadata)
}


load_metadata_safely <- function(metadata_path, sample_ids) {
  
  sample_ids <- trimws(as.character(sample_ids))
  
  if (!file.exists(metadata_path)) {
    message("Metadata file was not found. Minimal metadata will be created.")
    return(make_minimal_metadata(sample_ids))
  }
  
  metadata <- fread(
    metadata_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  if (nrow(metadata) == 0 || ncol(metadata) == 0) {
    message("Metadata file has zero dimensions. Minimal metadata will be created.")
    return(make_minimal_metadata(sample_ids))
  }
  
  if (!"SampleID" %in% colnames(metadata)) {
    message("SampleID column was not found. Using first column as SampleID.")
    colnames(metadata)[1] <- "SampleID"
  }
  
  metadata$SampleID <- trimws(as.character(metadata$SampleID))
  
  metadata <- metadata[
    metadata$SampleID != "" & !is.na(metadata$SampleID),
    ,
    drop = FALSE
  ]
  
  metadata <- metadata[
    !duplicated(metadata$SampleID),
    ,
    drop = FALSE
  ]
  
  if (nrow(metadata) == 0) {
    message("No valid SampleID values found. Minimal metadata will be created.")
    return(make_minimal_metadata(sample_ids))
  }
  
  rownames(metadata) <- metadata$SampleID
  
  missing_samples <- setdiff(sample_ids, rownames(metadata))
  
  if (length(missing_samples) > 0) {
    
    message(paste(
      "Adding minimal metadata for",
      length(missing_samples),
      "missing samples."
    ))
    
    missing_metadata <- make_minimal_metadata(missing_samples)
    
    for (col in setdiff(colnames(metadata), colnames(missing_metadata))) {
      missing_metadata[[col]] <- ""
    }
    
    for (col in setdiff(colnames(missing_metadata), colnames(metadata))) {
      metadata[[col]] <- ""
    }
    
    missing_metadata <- missing_metadata[, colnames(metadata), drop = FALSE]
    
    metadata <- rbind(metadata, missing_metadata)
  }
  
  if (!"Group" %in% colnames(metadata)) {
    metadata$Group <- "All_samples"
  }
  
  if (!"Site" %in% colnames(metadata)) {
    metadata$Site <- metadata$SampleID
  }
  
  rownames(metadata) <- metadata$SampleID
  
  metadata <- metadata[sample_ids, , drop = FALSE]
  
  return(metadata)
}


extract_sample_metadata <- function(ps_object) {
  
  metadata_current <- data.frame(
    sample_data(ps_object),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!"SampleID" %in% colnames(metadata_current)) {
    metadata_current$SampleID <- rownames(metadata_current)
  }
  
  return(metadata_current)
}


clean_taxonomy_table <- function(tax_df) {
  
  tax_df[] <- lapply(tax_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    return(x)
  })
  
  return(tax_df)
}


load_sequences_safely <- function(fasta_path, taxa_ids) {
  
  if (!file.exists(fasta_path)) {
    message("sequences.fasta was not found. The phyloseq object will be created without refseq.")
    return(NULL)
  }
  
  seqs <- readDNAStringSet(fasta_path)
  
  seqs <- seqs[names(seqs) %in% taxa_ids]
  
  if (length(seqs) == 0) {
    message("sequences.fasta does not contain compatible IDs. The phyloseq object will be created without refseq.")
    return(NULL)
  }
  
  return(seqs)
}


# =============================================================================
# 01. CREATE OR LOAD PHYLOSEQ OBJECT
# =============================================================================

if (file.exists(PS_FILE) && !FORCE_REBUILD_PS) {
  
  message("Existing phyloseq object found. Loading PS_FILE instead of rebuilding.")
  
  ps <- readRDS(PS_FILE)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  metadata_current <- extract_sample_metadata(ps)
  
  cat("\nLoaded phyloseq object:\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")
  
} else {
  
  message("No existing phyloseq object found, or FORCE_REBUILD_PS = TRUE.")
  message("Building phyloseq object from input files.")
  
  if (!file.exists(abundance_path)) {
    stop(paste("Required file not found:", abundance_path))
  }
  
  if (!file.exists(taxonomy_path)) {
    stop(paste("Required file not found:", taxonomy_path))
  }
  
  message("Loading abundance_table.csv with data.table::fread...")
  
  otu_dt <- fread(
    abundance_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  otu_ids <- otu_dt[[1]]
  otu_df <- otu_dt[, -1, drop = FALSE]
  rownames(otu_df) <- otu_ids
  
  rm(otu_dt)
  gc()
  
  otu_mat <- as.matrix(otu_df)
  
  rm(otu_df)
  gc()
  
  storage.mode(otu_mat) <- "numeric"
  otu_mat[is.na(otu_mat)] <- 0
  
  sample_ids <- trimws(as.character(colnames(otu_mat)))
  colnames(otu_mat) <- sample_ids
  
  cat("\nAbundance table diagnostics:\n")
  cat("Taxa:", nrow(otu_mat), "\n")
  cat("Samples:", ncol(otu_mat), "\n")
  cat("Total reads:", sum(otu_mat), "\n")
  cat("Samples with reads > 0:", sum(colSums(otu_mat) > 0), "\n")
  
  if (sum(otu_mat) == 0) {
    stop("The abundance table has a total read count of 0 after numeric conversion.")
  }
  
  message("Loading taxonomy_table.csv...")
  
  tax_dt <- fread(
    taxonomy_path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  tax_ids <- tax_dt[[1]]
  tax_df <- tax_dt[, -1, drop = FALSE]
  rownames(tax_df) <- tax_ids
  
  rm(tax_dt)
  gc()
  
  tax_df <- clean_taxonomy_table(tax_df)
  tax_mat <- as.matrix(tax_df)
  
  rm(tax_df)
  gc()
  
  message("Loading metadata...")
  
  metadata_current <- load_metadata_safely(
    metadata_path = metadata_path,
    sample_ids = sample_ids
  )
  
  common_samples <- intersect(
    colnames(otu_mat),
    rownames(metadata_current)
  )
  
  cat("\nSample matching diagnostics:\n")
  cat("Samples in abundance table:", ncol(otu_mat), "\n")
  cat("Samples in metadata:", nrow(metadata_current), "\n")
  cat("Shared samples:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    
    message("No shared sample IDs found. Rebuilding minimal metadata from abundance table.")
    
    metadata_current <- make_minimal_metadata(colnames(otu_mat))
    
    common_samples <- colnames(otu_mat)
  }
  
  otu_mat <- otu_mat[, common_samples, drop = FALSE]
  
  metadata_current <- metadata_current[common_samples, , drop = FALSE]
  
  if (nrow(metadata_current) == 0 || ncol(metadata_current) == 0) {
    stop("Metadata still has zero dimensions before creating the phyloseq object.")
  }
  
  common_taxa <- intersect(
    rownames(otu_mat),
    rownames(tax_mat)
  )
  
  if (length(common_taxa) == 0) {
    stop("No shared taxa IDs were found between abundance_table.csv and taxonomy_table.csv.")
  }
  
  otu_mat <- otu_mat[common_taxa, , drop = FALSE]
  
  tax_mat <- tax_mat[common_taxa, , drop = FALSE]
  
  cat("\nDiagnostics after ID matching:\n")
  cat("Taxa:", nrow(otu_mat), "\n")
  cat("Samples:", ncol(otu_mat), "\n")
  cat("Metadata rows:", nrow(metadata_current), "\n")
  cat("Metadata columns:", ncol(metadata_current), "\n")
  cat("Taxonomy rows:", nrow(tax_mat), "\n")
  cat("Total reads:", sum(otu_mat), "\n")
  
  seqs <- load_sequences_safely(
    fasta_path = fasta_path,
    taxa_ids = common_taxa
  )
  
  message("Creating phyloseq object...")
  
  if (!is.null(seqs)) {
    
    ps <- phyloseq(
      otu_table(otu_mat, taxa_are_rows = TRUE),
      tax_table(tax_mat),
      sample_data(metadata_current),
      refseq(seqs)
    )
    
  } else {
    
    ps <- phyloseq(
      otu_table(otu_mat, taxa_are_rows = TRUE),
      tax_table(tax_mat),
      sample_data(metadata_current)
    )
  }
  
  cat("\nBefore pruning:\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  metadata_current <- extract_sample_metadata(ps)
  
  cat("\nAfter pruning:\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")
  
  saveRDS(ps, PS_FILE)
  
  sample_depth <- data.frame(
    SampleID = sample_names(ps),
    Reads = sample_sums(ps),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    sample_depth,
    file.path(EXPLORER_DIR, "sample_read_depth.csv"),
    row.names = FALSE
  )
  
  write.csv(
    metadata_current,
    file.path(EXPLORER_DIR, "clean_metadata_used.csv"),
    row.names = FALSE
  )
  
  sink(file.path(EXPLORER_DIR, "ps_summary.txt"))
  
  cat("MetaDiv phyloseq summary\n")
  cat("========================\n\n")
  
  cat("Input directory:\n")
  cat(FOR_R_DIR, "\n\n")
  
  cat("Output directory:\n")
  cat(EXPLORER_DIR, "\n\n")
  
  cat("Saved phyloseq object:\n")
  cat(PS_FILE, "\n\n")
  
  cat("Number of taxa:", ntaxa(ps), "\n")
  cat("Number of samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n\n")
  
  cat("Taxonomic ranks available:\n")
  print(rank_names(ps))
  
  cat("\nSample variables available:\n")
  print(sample_variables(ps))
  
  cat("\nphyloseq object:\n")
  print(ps)
  
  sink()
  
  capture.output(
    sessionInfo(),
    file = file.path(EXPLORER_DIR, "R_sessionInfo_script01.txt")
  )
  
  cat("\nScript 01 completed successfully.\n")
}


# =============================================================================
# 02. QUICK EXPLORATORY ANALYSIS
# =============================================================================

OUT_DIR <- file.path(EXPLORER_DIR, "02_Explorer")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPLORER_DONE_FILE <- file.path(OUT_DIR, "MetaDiv_Explorer_DONE.txt")

if (file.exists(EXPLORER_DONE_FILE) && !FORCE_RERUN_EXPLORER) {
  
  message("Explorer outputs already exist. Skipping Script 02.")
  
} else {
  
  message("Running Script 02: Quick exploratory analysis.")
  
  TOP_N_TAXA <- 20
  
  DEFAULT_RANKS <- c(
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  
  ps <- readRDS(PS_FILE)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  metadata_current <- extract_sample_metadata(ps)
  
  cat("\nPhyloseq object loaded successfully for Script 02.\n")
  cat("Taxa:", ntaxa(ps), "\n")
  cat("Samples:", nsamples(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n")
  
  summary_table <- data.frame(
    Metric = c(
      "Number_of_taxa",
      "Number_of_samples",
      "Total_reads",
      "Minimum_reads_per_sample",
      "Maximum_reads_per_sample",
      "Mean_reads_per_sample",
      "Median_reads_per_sample"
    ),
    Value = c(
      ntaxa(ps),
      nsamples(ps),
      sum(sample_sums(ps)),
      min(sample_sums(ps)),
      max(sample_sums(ps)),
      mean(sample_sums(ps)),
      median(sample_sums(ps))
    )
  )
  
  write.csv(
    summary_table,
    file.path(OUT_DIR, "MetaDiv_Explorer_Summary.csv"),
    row.names = FALSE
  )
  
  capture.output(
    ps,
    file = file.path(OUT_DIR, "phyloseq_object_print.txt")
  )
  
  sample_depth <- data.frame(
    SampleID = sample_names(ps),
    Reads = sample_sums(ps),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    sample_depth,
    file.path(OUT_DIR, "sample_read_depth.csv"),
    row.names = FALSE
  )
  
  p_depth_bar <- ggplot(
    sample_depth,
    aes(x = reorder(SampleID, Reads), y = Reads)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Sequencing depth per sample",
      x = "Sample",
      y = "Reads"
    )
  
  ggsave(
    file.path(OUT_DIR, "sample_read_depth_barplot.png"),
    p_depth_bar,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  p_depth_hist <- ggplot(
    sample_depth,
    aes(x = Reads)
  ) +
    geom_histogram(bins = 30) +
    theme_bw() +
    labs(
      title = "Read depth distribution",
      x = "Reads",
      y = "Number of samples"
    )
  
  ggsave(
    file.path(OUT_DIR, "sample_read_depth_histogram.png"),
    p_depth_hist,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  taxa_abundance <- data.frame(
    TaxonID = taxa_names(ps),
    Total_Abundance = taxa_sums(ps),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(Total_Abundance))
  
  write.csv(
    taxa_abundance,
    file.path(OUT_DIR, "taxa_total_abundance.csv"),
    row.names = FALSE
  )
  
  top_taxa_table <- taxa_abundance %>%
    slice_head(n = TOP_N_TAXA)
  
  write.csv(
    top_taxa_table,
    file.path(OUT_DIR, paste0("top_", TOP_N_TAXA, "_taxa_by_abundance.csv")),
    row.names = FALSE
  )
  
  otu_mat_explorer <- as(otu_table(ps), "matrix")
  
  if (!taxa_are_rows(ps)) {
    otu_mat_explorer <- t(otu_mat_explorer)
  }
  
  taxa_prevalence <- data.frame(
    TaxonID = rownames(otu_mat_explorer),
    Prevalence = rowSums(otu_mat_explorer > 0),
    Total_Abundance = rowSums(otu_mat_explorer),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Prevalence_Percent = Prevalence / nsamples(ps) * 100
    ) %>%
    arrange(desc(Prevalence), desc(Total_Abundance))
  
  write.csv(
    taxa_prevalence,
    file.path(OUT_DIR, "taxa_prevalence.csv"),
    row.names = FALSE
  )
  
  p_prev <- ggplot(
    taxa_prevalence,
    aes(x = Prevalence_Percent)
  ) +
    geom_histogram(bins = 30) +
    theme_bw() +
    labs(
      title = "Taxa prevalence distribution",
      x = "Prevalence (%)",
      y = "Number of taxa"
    )
  
  ggsave(
    file.path(OUT_DIR, "taxa_prevalence_histogram.png"),
    p_prev,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  available_ranks <- rank_names(ps)
  
  ranks_to_plot <- DEFAULT_RANKS[
    DEFAULT_RANKS %in% available_ranks
  ]
  
  plot_taxonomic_barplot <- function(ps_object, rank_name, top_n = 20, out_dir = OUT_DIR) {
    
    cat("Generating barplot for rank:", rank_name, "\n")
    
    ps_glom <- tax_glom(
      ps_object,
      taxrank = rank_name,
      NArm = FALSE
    )
    
    ps_rel <- transform_sample_counts(
      ps_glom,
      function(x) {
        if (sum(x) == 0) return(x)
        x / sum(x)
      }
    )
    
    tax_df <- as.data.frame(tax_table(ps_rel))
    
    tax_df$TaxonID <- rownames(tax_df)
    
    abund_df <- data.frame(
      TaxonID = taxa_names(ps_rel),
      Total_Abundance = taxa_sums(ps_rel),
      stringsAsFactors = FALSE
    )
    
    tax_abund <- left_join(
      abund_df,
      tax_df,
      by = "TaxonID"
    )
    
    top_labels <- tax_abund %>%
      arrange(desc(Total_Abundance)) %>%
      pull(.data[[rank_name]]) %>%
      unique()
    
    top_labels <- top_labels[
      !is.na(top_labels) & top_labels != ""
    ]
    
    top_labels <- head(top_labels, top_n)
    
    ps_melt <- psmelt(ps_rel)
    
    ps_melt[[rank_name]] <- as.character(ps_melt[[rank_name]])
    
    ps_melt[[rank_name]][
      is.na(ps_melt[[rank_name]]) | ps_melt[[rank_name]] == ""
    ] <- "Unclassified"
    
    ps_melt$TaxonGroup <- ifelse(
      ps_melt[[rank_name]] %in% top_labels,
      ps_melt[[rank_name]],
      "Other"
    )
    
    plot_df <- ps_melt %>%
      group_by(Sample, TaxonGroup) %>%
      summarise(
        Abundance = sum(Abundance),
        .groups = "drop"
      )
    
    p <- ggplot(
      plot_df,
      aes(x = Sample, y = Abundance, fill = TaxonGroup)
    ) +
      geom_col() +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "right"
      ) +
      labs(
        title = paste("Relative composition -", rank_name),
        x = "Sample",
        y = "Relative abundance",
        fill = rank_name
      )
    
    ggsave(
      file.path(out_dir, paste0("barplot_top_", top_n, "_", rank_name, ".png")),
      p,
      width = 14,
      height = 7,
      dpi = 300
    )
    
    write.csv(
      plot_df,
      file.path(out_dir, paste0("barplot_data_", rank_name, ".csv")),
      row.names = FALSE
    )
    
    return(p)
  }
  
  for (rank_name in ranks_to_plot) {
    
    plot_taxonomic_barplot(
      ps_object = ps,
      rank_name = rank_name,
      top_n = TOP_N_TAXA,
      out_dir = OUT_DIR
    )
  }
  
  for (rank_name in ranks_to_plot) {
    
    ps_glom <- tax_glom(
      ps,
      taxrank = rank_name,
      NArm = FALSE
    )
    
    tax_df <- as.data.frame(tax_table(ps_glom))
    
    tax_df$TaxonID <- rownames(tax_df)
    
    abund_df <- data.frame(
      TaxonID = taxa_names(ps_glom),
      Total_Abundance = taxa_sums(ps_glom),
      stringsAsFactors = FALSE
    )
    
    rank_summary <- left_join(
      abund_df,
      tax_df,
      by = "TaxonID"
    ) %>%
      mutate(
        Rank = rank_name,
        Taxon = .data[[rank_name]]
      ) %>%
      group_by(Rank, Taxon) %>%
      summarise(
        Total_Abundance = sum(Total_Abundance),
        .groups = "drop"
      ) %>%
      arrange(desc(Total_Abundance))
    
    write.csv(
      rank_summary,
      file.path(OUT_DIR, paste0("taxonomic_summary_", rank_name, ".csv")),
      row.names = FALSE
    )
  }
  
  if (taxa_are_rows(ps)) {
    otu_for_vegan <- t(as(otu_table(ps), "matrix"))
  } else {
    otu_for_vegan <- as(otu_table(ps), "matrix")
  }
  
  spec_curve <- specaccum(
    otu_for_vegan,
    method = "random"
  )
  
  png(
    filename = file.path(OUT_DIR, "species_accumulation_curve.png"),
    width = 1400,
    height = 1000,
    res = 150
  )
  
  plot(
    spec_curve,
    xlab = "Number of samples",
    ylab = "Accumulated richness",
    main = "Taxa accumulation curve"
  )
  
  dev.off()
  
  spec_curve_df <- data.frame(
    Sites = spec_curve$sites,
    Richness = spec_curve$richness,
    SD = spec_curve$sd
  )
  
  write.csv(
    spec_curve_df,
    file.path(OUT_DIR, "species_accumulation_curve_data.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ps,
    file.path(OUT_DIR, "ps_explorer_clean.rds")
  )
  
  write.csv(
    metadata_current,
    file.path(OUT_DIR, "metadata_used_in_explorer.csv"),
    row.names = FALSE
  )
  
  write.csv(
    as.data.frame(tax_table(ps)),
    file.path(OUT_DIR, "taxonomy_used_in_explorer.csv"),
    row.names = TRUE
  )
  
  capture.output(
    sessionInfo(),
    file = file.path(OUT_DIR, "R_sessionInfo_script02.txt")
  )
  
  writeLines(
    paste("Explorer completed:", Sys.time()),
    EXPLORER_DONE_FILE
  )
  
  cat("\nMetaDiv Explorer completed successfully.\n")
  cat("\nScript 02 outputs saved in:\n")
  cat(OUT_DIR, "\n")
}


# =============================================================================
# 03. ALPHA DIVERSITY ANALYSIS
# =============================================================================

ALPHA_DIR <- file.path(EXPLORER_DIR, "03_Alpha_Diversity")

dir.create(
  ALPHA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

ALPHA_DONE_FILE <- file.path(ALPHA_DIR, "MetaDiv_Alpha_Diversity_DONE.txt")

ALPHA_MAIN_FILE <- file.path(ALPHA_DIR, "alpha_diversity_indices.csv")

if (file.exists(ALPHA_DONE_FILE) && file.exists(ALPHA_MAIN_FILE) && !FORCE_RERUN_ALPHA) {
  
  message("Alpha diversity outputs already exist. Skipping Script 03.")
  
} else {
  
  message("Running Script 03: Alpha diversity analysis.")
  
  ps <- readRDS(PS_FILE)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  metadata_current <- extract_sample_metadata(ps)
  
  otu_alpha <- as(otu_table(ps), "matrix")
  
  if (taxa_are_rows(ps)) {
    otu_alpha <- t(otu_alpha)
  }
  
  otu_alpha[is.na(otu_alpha)] <- 0
  
  otu_alpha <- otu_alpha[
    rowSums(otu_alpha) > 0,
    colSums(otu_alpha) > 0,
    drop = FALSE
  ]
  
  cat("\nAlpha diversity input diagnostics:\n")
  cat("Samples:", nrow(otu_alpha), "\n")
  cat("Taxa:", ncol(otu_alpha), "\n")
  cat("Total reads:", sum(otu_alpha), "\n")
  
  observed_richness <- specnumber(otu_alpha)
  
  shannon <- diversity(
    otu_alpha,
    index = "shannon"
  )
  
  simpson <- diversity(
    otu_alpha,
    index = "simpson"
  )
  
  inverse_simpson <- diversity(
    otu_alpha,
    index = "invsimpson"
  )
  
  pielou_evenness <- shannon / log(observed_richness)
  
  pielou_evenness[
    is.nan(pielou_evenness) | is.infinite(pielou_evenness)
  ] <- NA
  
  estimate_r <- estimateR(otu_alpha)
  
  cat("\nestimateR output dimensions:\n")
  print(dim(estimate_r))
  
  cat("\nestimateR row names:\n")
  print(rownames(estimate_r))
  
  chao1 <- estimate_r["S.chao1", ]
  ace   <- estimate_r["S.ACE", ]
  
  alpha_df <- data.frame(
    SampleID = rownames(otu_alpha),
    Reads = rowSums(otu_alpha),
    Observed_Richness = as.numeric(observed_richness),
    Shannon = as.numeric(shannon),
    Simpson = as.numeric(simpson),
    Inverse_Simpson = as.numeric(inverse_simpson),
    Pielou_Evenness = as.numeric(pielou_evenness),
    Chao1 = as.numeric(chao1),
    ACE = as.numeric(ace),
    stringsAsFactors = FALSE
  )
  
  alpha_df <- alpha_df %>%
    left_join(
      metadata_current,
      by = "SampleID"
    )
  
  write.csv(
    alpha_df,
    ALPHA_MAIN_FILE,
    row.names = FALSE
  )
  
  alpha_summary <- alpha_df %>%
    summarise(
      Samples = n(),
      Min_Reads = min(Reads, na.rm = TRUE),
      Max_Reads = max(Reads, na.rm = TRUE),
      Mean_Reads = mean(Reads, na.rm = TRUE),
      Median_Reads = median(Reads, na.rm = TRUE),
      Mean_Observed_Richness = mean(Observed_Richness, na.rm = TRUE),
      Median_Observed_Richness = median(Observed_Richness, na.rm = TRUE),
      Mean_Shannon = mean(Shannon, na.rm = TRUE),
      Median_Shannon = median(Shannon, na.rm = TRUE),
      Mean_Simpson = mean(Simpson, na.rm = TRUE),
      Median_Simpson = median(Simpson, na.rm = TRUE),
      Mean_Inverse_Simpson = mean(Inverse_Simpson, na.rm = TRUE),
      Median_Inverse_Simpson = median(Inverse_Simpson, na.rm = TRUE),
      Mean_Pielou_Evenness = mean(Pielou_Evenness, na.rm = TRUE),
      Median_Pielou_Evenness = median(Pielou_Evenness, na.rm = TRUE),
      Mean_Chao1 = mean(Chao1, na.rm = TRUE),
      Median_Chao1 = median(Chao1, na.rm = TRUE),
      Mean_ACE = mean(ACE, na.rm = TRUE),
      Median_ACE = median(ACE, na.rm = TRUE)
    )
  
  write.csv(
    alpha_summary,
    file.path(ALPHA_DIR, "alpha_diversity_summary.csv"),
    row.names = FALSE
  )
  
  alpha_long <- alpha_df %>%
    select(
      SampleID,
      Observed_Richness,
      Shannon,
      Simpson,
      Inverse_Simpson,
      Pielou_Evenness,
      Chao1,
      ACE
    ) %>%
    pivot_longer(
      cols = -SampleID,
      names_to = "Index",
      values_to = "Value"
    )
  
  p_alpha_hist <- ggplot(
    alpha_long,
    aes(x = Value)
  ) +
    geom_histogram(bins = 30) +
    facet_wrap(
      ~ Index,
      scales = "free"
    ) +
    theme_bw() +
    labs(
      title = "Alpha diversity index distributions",
      x = "Value",
      y = "Number of samples"
    )
  
  ggsave(
    file.path(ALPHA_DIR, "alpha_diversity_histograms.png"),
    p_alpha_hist,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  if ("Group" %in% colnames(alpha_df)) {
    
    p_alpha_box <- alpha_df %>%
      select(
        SampleID,
        Group,
        Observed_Richness,
        Shannon,
        Simpson,
        Inverse_Simpson,
        Pielou_Evenness,
        Chao1,
        ACE
      ) %>%
      pivot_longer(
        cols = c(
          Observed_Richness,
          Shannon,
          Simpson,
          Inverse_Simpson,
          Pielou_Evenness,
          Chao1,
          ACE
        ),
        names_to = "Index",
        values_to = "Value"
      ) %>%
      ggplot(
        aes(x = Group, y = Value)
      ) +
      geom_boxplot(
        outlier.size = 0.5
      ) +
      facet_wrap(
        ~ Index,
        scales = "free_y"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      ) +
      labs(
        title = "Alpha diversity by group",
        x = "Group",
        y = "Value"
      )
    
    ggsave(
      file.path(ALPHA_DIR, "alpha_diversity_boxplots_by_group.png"),
      p_alpha_box,
      width = 14,
      height = 9,
      dpi = 300
    )
    
  } else {
    
    message("Column 'Group' was not found in metadata. Skipping alpha boxplots by group.")
  }
  
  capture.output(
    sessionInfo(),
    file = file.path(ALPHA_DIR, "R_sessionInfo_script03_alpha_diversity.txt")
  )
  
  writeLines(
    paste("Alpha diversity completed:", Sys.time()),
    ALPHA_DONE_FILE
  )
  
  cat("\nAlpha diversity analysis completed successfully.\n")
  cat("\nScript 03 outputs saved in:\n")
  cat(ALPHA_DIR, "\n")
}


# =============================================================================
# FINAL PIPELINE MESSAGE
# =============================================================================

cat("\nPipeline finished.\n")
cat("\nResume options:\n")
cat("\nFORCE_REBUILD_PS =", FORCE_REBUILD_PS, "\n")
cat("FORCE_RERUN_EXPLORER =", FORCE_RERUN_EXPLORER, "\n")
cat("FORCE_RERUN_ALPHA =", FORCE_RERUN_ALPHA, "\n")