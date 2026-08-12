# =============================================================================
# 03_MetaDiv_iNEXT_Hill_Numbers_portable.R
# =============================================================================
# MetaDiv Explorer
#
# Hill numbers and iNEXT diversity analysis.
#
# q = 0 : Species richness / Observed richness
# q = 1 : Exponential Shannon diversity
# q = 2 : Inverse Simpson diversity
#
# Author:
# Bernardo Águila
# Instituto de Biología, UNAM
# =============================================================================


# =============================================================================
# 1. PORTABLE STARTUP AND CONFIGURATION
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
    
    normalizePath(
      getwd(),
      winslash = "/",
      mustWork = TRUE
    )
    
  }
  
})

source(file.path(SCRIPT_DIR, "check_environment.R"))

check_metadiv_environment()


windows_path_to_wsl <- function(path){
  
  if(
    .Platform$OS.type != "windows" &&
    grepl("^[A-Za-z]:[/\\\\]", path)
  ){
    
    drive <- tolower(substr(path,1,1))
    
    remainder <- substring(path,3)
    
    remainder <- gsub("\\\\","/",remainder)
    
    path <- paste0("/mnt/",drive,remainder)
    
  }
  
  path
  
}


resolve_for_r_dir <- function(){
  
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args)>=1 && nzchar(args[1])){
    
    selected_dir <- args[1]
    
  }else if(
    exists("FOR_R_DIR",envir=.GlobalEnv) &&
    nzchar(get("FOR_R_DIR",envir=.GlobalEnv))
  ){
    
    selected_dir <- get("FOR_R_DIR",envir=.GlobalEnv)
    
  }else if(
    .Platform$OS.type=="windows" &&
    interactive()
  ){
    
    selected_dir <- choose.dir(
      caption="Select MetaDiv For_R dataset"
    )
    
  }else if(interactive()){
    
    selected_dir <- readline(
      "Enter MetaDiv For_R dataset directory: "
    )
    
  }else{
    
    stop(
      "Dataset directory not provided.\n",
      "Usage:\n",
      "Rscript 03_MetaDiv_iNEXT_Hill_Numbers_portable.R /path/to/For_R/run"
    )
    
  }
  
  selected_dir <- windows_path_to_wsl(selected_dir)
  
  normalizePath(
    selected_dir,
    winslash="/",
    mustWork=TRUE
  )
  
}


FOR_R_DIR <- resolve_for_r_dir()

# All R modules use the same output directory created by Module 01.
BASE_OUTPUT_DIR <- file.path(
  FOR_R_DIR,
  "R_analysis_output"
)

EXPLORER_DIR <- BASE_OUTPUT_DIR

PS_FILE <- file.path(
  EXPLORER_DIR,
  "ps_metadiv_builder.rds"
)

INEXT_DIR <- file.path(
  EXPLORER_DIR,
  "03_iNEXT_Hill_Numbers"
)

dir.create(
  INEXT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\nMetaDiv iNEXT configuration loaded.\n")
cat("Input directory:\n", FOR_R_DIR, "\n")
cat("Output directory:\n", INEXT_DIR, "\n")

# =============================================================================
# 2. LOAD PACKAGES
# =============================================================================

library(phyloseq)
library(iNEXT)
library(dplyr)
library(tidyr)
library(ggplot2)


# =============================================================================
# 3. LOAD PHYLOSEQ OBJECT
# =============================================================================

if (!file.exists(PS_FILE)) {
  stop(paste("Phyloseq object not found:", PS_FILE))
}

ps <- readRDS(PS_FILE)

ps <- prune_taxa(
  taxa_sums(ps) > 0,
  ps
)

ps <- prune_samples(
  sample_sums(ps) > 0,
  ps
)

cat("\nPhyloseq object loaded.\n")
cat("Taxa:", ntaxa(ps), "\n")
cat("Samples:", nsamples(ps), "\n")
cat("Total reads:", sum(sample_sums(ps)), "\n")


# =============================================================================
# 4. PREPARE ABUNDANCE DATA FOR iNEXT
# =============================================================================

otu_mat <- as(
  otu_table(ps),
  "matrix"
)

# iNEXT expects each community/sample as a vector of abundances.
# Therefore, we need samples as rows and taxa as columns.

if (taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

otu_mat[is.na(otu_mat)] <- 0

otu_mat <- otu_mat[
  rowSums(otu_mat) > 0,
  colSums(otu_mat) > 0,
  drop = FALSE
]

cat("\niNEXT input diagnostics:\n")
cat("Samples:", nrow(otu_mat), "\n")
cat("Taxa:", ncol(otu_mat), "\n")


# =============================================================================
# 5. FAST HILL NUMBERS DIRECTLY FROM ABUNDANCE TABLE
# =============================================================================
# These are the same diversity orders used by iNEXT:
# q = 0 observed richness
# q = 1 exp(Shannon)
# q = 2 inverse Simpson

hill_df <- data.frame(
  SampleID = rownames(otu_mat),
  Reads = rowSums(otu_mat),
  Hill_q0_Richness = rowSums(otu_mat > 0),
  stringsAsFactors = FALSE
)

relative_abundance <- otu_mat / rowSums(otu_mat)

relative_abundance[is.na(relative_abundance)] <- 0

shannon_entropy <- apply(
  relative_abundance,
  1,
  function(x) {
    x <- x[x > 0]
    -sum(x * log(x))
  }
)

simpson_concentration <- rowSums(
  relative_abundance ^ 2
)

hill_df$Hill_q1_ExpShannon <- exp(shannon_entropy)

hill_df$Hill_q2_InvSimpson <- 1 / simpson_concentration

write.csv(
  hill_df,
  file.path(INEXT_DIR, "hill_numbers_observed_per_sample.csv"),
  row.names = FALSE
)


# =============================================================================
# 6. PREPARE LIST FOR iNEXT
# =============================================================================
# For thousands of sites this can be heavy.
# With 326 samples it should run, but it may take time.
# Each list element is one sample/community.

inext_list <- lapply(
  seq_len(nrow(otu_mat)),
  function(i) {
    x <- otu_mat[i, ]
    x <- x[x > 0]
    as.numeric(x)
  }
)

names(inext_list) <- rownames(otu_mat)


# =============================================================================
# 7. RUN iNEXT
# =============================================================================

q_orders <- c(0, 1, 2)

cat("\nRunning iNEXT. This may take some time...\n")

inext_result <- iNEXT(
  inext_list,
  q = q_orders,
  datatype = "abundance",
  se = TRUE,
  nboot = 50
)

saveRDS(
  inext_result,
  file.path(INEXT_DIR, "iNEXT_result.rds")
)


# =============================================================================
# 8. EXPORT iNEXT TABLES
# =============================================================================

if (!is.null(inext_result$iNextEst)) {
  
  inext_estimates <- inext_result$iNextEst
  
  if (is.list(inext_estimates)) {
    
    for (name in names(inext_estimates)) {
      
      write.csv(
        inext_estimates[[name]],
        file.path(INEXT_DIR, paste0("iNEXT_estimates_", name, ".csv")),
        row.names = FALSE
      )
    }
  }
}

if (!is.null(inext_result$AsyEst)) {
  
  write.csv(
    inext_result$AsyEst,
    file.path(INEXT_DIR, "iNEXT_asymptotic_estimates.csv"),
    row.names = FALSE
  )
}


# =============================================================================
# 9. EXPORT OBSERVED iNEXT ASYMPTOTIC SUMMARY
# =============================================================================

if (!is.null(inext_result$AsyEst)) {
  
  asy <- inext_result$AsyEst
  
  write.csv(
    asy,
    file.path(INEXT_DIR, "hill_numbers_iNEXT_asymptotic_summary.csv"),
    row.names = FALSE
  )
}


# =============================================================================
# 10. PLOTS
# =============================================================================

p_sample_size <- ggiNEXT(
  inext_result,
  type = 1,
  facet.var = "Order.q"
) +
  theme_bw() +
  labs(
    title = "iNEXT rarefaction / extrapolation by sample size"
  )

ggsave(
  file.path(INEXT_DIR, "iNEXT_sample_size_curve.png"),
  p_sample_size,
  width = 12,
  height = 8,
  dpi = 300
)

p_coverage <- ggiNEXT(
  inext_result,
  type = 3,
  facet.var = "Order.q"
) +
  theme_bw() +
  labs(
    title = "iNEXT diversity by sample coverage"
  )

ggsave(
  file.path(INEXT_DIR, "iNEXT_sample_coverage_curve.png"),
  p_coverage,
  width = 12,
  height = 8,
  dpi = 300
)


# =============================================================================
# 11. SESSION INFO
# =============================================================================

capture.output(
  sessionInfo(),
  file = file.path(INEXT_DIR, "R_sessionInfo_module03_iNEXT.txt")
)


# =============================================================================
# 12. FINAL MESSAGE
# =============================================================================

cat("\niNEXT Hill number analysis completed successfully.\n")
cat("\nOutputs saved in:\n")
cat(INEXT_DIR, "\n")