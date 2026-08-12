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
# 01_MetaDiv_explorer_alpha.R
# =============================================================================
# Loads the phyloseq object created by 00_MetaDiv_create_phyloseq.R and runs:
#   - exploratory summaries
#   - taxonomic barplots through genus
#   - taxa accumulation
#   - alpha-diversity indices
#
# Author: Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================

FOR_R_DIR <- resolve_for_r_dir()

OUTPUT_ROOT <- file.path(FOR_R_DIR, "R_analysis_output")
PS_FILE <- file.path(OUTPUT_ROOT, "ps_metadiv_builder.rds")

FORCE_RERUN_EXPLORER <- FALSE
FORCE_RERUN_ALPHA <- FALSE

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(readr)
  library(vegan)
})

if (!file.exists(PS_FILE)) {
  stop(
    "Phyloseq object not found:\n",
    PS_FILE,
    "\nRun 00_MetaDiv_create_phyloseq.R first.",
    call. = FALSE
  )
}

ps <- readRDS(PS_FILE)

if (!inherits(ps, "phyloseq")) {
  stop("The saved RDS object is not a phyloseq object.", call. = FALSE)
}

ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

extract_sample_metadata <- function(ps_object) {
  metadata_current <- data.frame(
    sample_data(ps_object),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!"SampleID" %in% colnames(metadata_current)) {
    metadata_current$SampleID <- rownames(metadata_current)
  }
  
  metadata_current
}

metadata_current <- extract_sample_metadata(ps)

cat("\nMetaDiv Module 01 configuration loaded.\n")
cat("Input directory:\n", FOR_R_DIR, "\n")
cat("Phyloseq object:\n", PS_FILE, "\n")
cat("Taxa:", ntaxa(ps), "\n")
cat("Samples:", nsamples(ps), "\n")
cat("Total reads:", sum(sample_sums(ps)), "\n")

# =============================================================================
# 1. QUICK EXPLORATORY ANALYSIS
# =============================================================================

OUT_DIR <- file.path(OUTPUT_ROOT, "01_Explorer")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPLORER_DONE_FILE <- file.path(OUT_DIR, "MetaDiv_Explorer_DONE.txt")

if (file.exists(EXPLORER_DONE_FILE) && !FORCE_RERUN_EXPLORER) {
  
  message("Explorer outputs already exist. Skipping exploratory analysis.")
  
} else {
  
  message("Running Module 01 exploratory analysis.")
  
  TOP_N_TAXA <- 20
  
  # Taxonomic barplots are intentionally limited to genus.
  # Species-level plots can become extremely large and slow in metabarcoding
  # datasets and are not suitable for the initial exploratory report.
  MAX_BARPLOT_RANK <- "genus"
  
  DEFAULT_RANKS <- c(
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus"
  )
  
  # Sample names are hidden when there are too many samples to display
  # legibly. Complete sample-level values remain available in CSV files.
  MAX_SAMPLE_LABELS <- 60
  
  ps <- readRDS(PS_FILE)
  
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  metadata_current <- extract_sample_metadata(ps)
  
  cat("\nPhyloseq object loaded successfully for Module 01.\n")
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
  
  show_sample_labels <- nrow(sample_depth) <= MAX_SAMPLE_LABELS
  
  p_depth_bar <- ggplot(
    sample_depth,
    aes(x = reorder(SampleID, Reads), y = Reads)
  ) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    theme(
      axis.text.y = if (show_sample_labels) {
        element_text(size = 7)
      } else {
        element_blank()
      },
      axis.ticks.y = if (show_sample_labels) {
        element_line()
      } else {
        element_blank()
      }
    ) +
    labs(
      title = "Sequencing depth per sample",
      subtitle = if (show_sample_labels) {
        NULL
      } else {
        paste0(
          nrow(sample_depth),
          " samples: labels hidden for readability; see sample_read_depth.csv"
        )
      },
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
  
  cat(
    "
Taxonomic barplot limit:",
    MAX_BARPLOT_RANK,
    "
Ranks to plot:",
    paste(ranks_to_plot, collapse = ", "),
    "
"
  )
  
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
    
    sample_count <- dplyr::n_distinct(plot_df$Sample)
    show_sample_labels <- sample_count <= MAX_SAMPLE_LABELS
    
    p <- ggplot(
      plot_df,
      aes(x = Sample, y = Abundance, fill = TaxonGroup)
    ) +
      geom_col(width = 1) +
      theme_bw() +
      theme(
        axis.text.x = if (show_sample_labels) {
          element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 6
          )
        } else {
          element_blank()
        },
        axis.ticks.x = if (show_sample_labels) {
          element_line()
        } else {
          element_blank()
        },
        legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8)
      ) +
      labs(
        title = paste("Relative composition -", rank_name),
        subtitle = if (show_sample_labels) {
          paste0("Top ", top_n, " taxa; remaining taxa grouped as Other")
        } else {
          paste0(
            "Top ", top_n, " taxa; ", sample_count,
            " sample labels hidden for readability; see barplot_data_",
            rank_name, ".csv"
          )
        },
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
    file = file.path(OUT_DIR, "R_sessionInfo_module01_explorer.txt")
  )
  
  writeLines(
    paste("Explorer completed:", Sys.time()),
    EXPLORER_DONE_FILE
  )
  
  cat("\nMetaDiv Explorer completed successfully.\n")
  cat("\nModule 01 Explorer outputs saved in:\n")
  cat(OUT_DIR, "\n")
}


# =============================================================================
# 2. ALPHA DIVERSITY ANALYSIS
# =============================================================================

ALPHA_DIR <- file.path(OUTPUT_ROOT, "01_Alpha_Diversity")

dir.create(
  ALPHA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

ALPHA_DONE_FILE <- file.path(ALPHA_DIR, "MetaDiv_Alpha_Diversity_DONE.txt")

ALPHA_MAIN_FILE <- file.path(ALPHA_DIR, "alpha_diversity_indices.csv")

if (file.exists(ALPHA_DONE_FILE) && file.exists(ALPHA_MAIN_FILE) && !FORCE_RERUN_ALPHA) {
  
  message("Alpha diversity outputs already exist. Skipping alpha-diversity analysis.")
  
} else {
  
  message("Running Module 01 alpha diversity analysis.")
  
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
    file = file.path(ALPHA_DIR, "R_sessionInfo_module01_alpha_diversity.txt")
  )
  
  writeLines(
    paste("Alpha diversity completed:", Sys.time()),
    ALPHA_DONE_FILE
  )
  
  cat("\nAlpha diversity analysis completed successfully.\n")
  cat("\nModule 01 alpha-diversity outputs saved in:\n")
  cat(ALPHA_DIR, "\n")
}



# =============================================================================
# FINAL MESSAGE
# =============================================================================

cat("\nModule 01 finished.\n")
cat("FORCE_RERUN_EXPLORER =", FORCE_RERUN_EXPLORER, "\n")
cat("FORCE_RERUN_ALPHA =", FORCE_RERUN_ALPHA, "\n")