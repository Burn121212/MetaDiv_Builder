# =============================================================================
# check_environment.R
# MetaDiv Builder
# Verifies and installs required R packages if they are missing.
# =============================================================================

check_metadiv_environment <- function() {
  
  options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    timeout = max(600, getOption("timeout"))
  )
  
  cran_packages <- c(
    "data.table",
    "dplyr",
    "tidyr",
    "tibble",
    "readr",
    "ggplot2",
    "vegan",
    "iNEXT"
  )
  
  bioc_packages <- c(
    "phyloseq",
    "Biostrings"
  )
  
  cat("\n========================================\n")
  cat("MetaDiv Builder - Environment Check\n")
  cat("========================================\n")
  
  # -----------------------------
  # Install missing CRAN packages
  # -----------------------------
  missing_cran <- cran_packages[
    !sapply(cran_packages, requireNamespace, quietly = TRUE)
  ]
  
  if(length(missing_cran) > 0){
    
    cat("\nInstalling CRAN packages...\n")
    
    install.packages(
      missing_cran,
      dependencies = TRUE
    )
    
  }
  
  # -----------------------------
  # Install BiocManager
  # -----------------------------
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    
    install.packages("BiocManager")
    
  }
  
  # -----------------------------
  # Install missing Bioconductor packages
  # -----------------------------
  missing_bioc <- bioc_packages[
    !sapply(bioc_packages, requireNamespace, quietly = TRUE)
  ]
  
  if(length(missing_bioc) > 0){
    
    cat("\nInstalling Bioconductor packages...\n")
    
    BiocManager::install(
      missing_bioc,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
  # -----------------------------
  # Final verification
  # -----------------------------
  all_packages <- c(
    cran_packages,
    bioc_packages
  )
  
  cat("\nPackage status:\n\n")
  
  for(pkg in all_packages){
    
    if(requireNamespace(pkg, quietly = TRUE)){
      
      cat("[OK] ", pkg, "\n", sep="")
      
    }else{
      
      stop(
        paste("Package could not be installed:", pkg),
        call.=FALSE
      )
      
    }
    
  }
  
  cat("\nEnvironment ready.\n")
  cat("========================================\n\n")
  
}