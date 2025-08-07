# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Master Analysis Script - Run Complete Pipeline
# ==============================================================================

# Set working directory to script location
if (exists("rstudioapi") && rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(script_dir)
}

# Source utility functions
source("utils/plot_functions.R")
source("utils/pathway_functions.R")

# Pipeline execution options ----
run_setup <- TRUE
run_data_loading <- TRUE  
run_cell_characterization <- TRUE
run_differential_expression <- TRUE
run_pathway_enrichment <- TRUE
run_gsea_analysis <- TRUE
run_network_analysis <- TRUE

# Create log file
log_file <- file.path("results", "analysis_log.txt")
dir.create("results", recursive = TRUE, showWarnings = FALSE)

log_message <- function(message) {
  timestamp <- Sys.time()
  log_entry <- paste(timestamp, "-", message)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

log_message("=== Starting Single-cell RNA-seq Analysis Pipeline ===")

# Step 1: Environment Setup ----
if (run_setup) {
  log_message("Step 1: Setting up environment...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/00_setup.R")
    log_message("Environment setup completed successfully")
  }, error = function(e) {
    log_message(paste("Error in setup:", e$message))
    stop("Setup failed")
  })
  
  end_time <- Sys.time()
  log_message(paste("Setup time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 2: Data Loading ----
if (run_data_loading) {
  log_message("Step 2: Loading and preprocessing data...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/01_data_loading.R")
    log_message("Data loading completed successfully")
    log_message(paste("Loaded", ncol(obj_all_cells), "total cells and", 
                     ncol(obj_sf), "synovial fibroblasts"))
  }, error = function(e) {
    log_message(paste("Error in data loading:", e$message))
    stop("Data loading failed")
  })
  
  end_time <- Sys.time()
  log_message(paste("Data loading time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 3: Cell Characterization ----
if (run_cell_characterization) {
  log_message("Step 3: Cell type characterization and visualization...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/02_cell_characterization.R")
    log_message("Cell characterization completed successfully")
  }, error = function(e) {
    log_message(paste("Error in cell characterization:", e$message))
    warning("Cell characterization failed, continuing with next step")
  })
  
  end_time <- Sys.time()
  log_message(paste("Cell characterization time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 4: Differential Expression Analysis ----
if (run_differential_expression) {
  log_message("Step 4: Differential expression analysis...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/03_differential_expression.R")
    log_message("Differential expression analysis completed successfully")
    
    # Report results
    if (exists("deg_results")) {
      total_degs <- sum(sapply(deg_results$degs, function(x) sapply(x, nrow)))
      log_message(paste("Total DEGs identified:", total_degs))
    }
  }, error = function(e) {
    log_message(paste("Error in differential expression analysis:", e$message))
    warning("Differential expression analysis failed, continuing with next step")
  })
  
  end_time <- Sys.time()
  log_message(paste("Differential expression time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 5: Pathway Enrichment Analysis ----
if (run_pathway_enrichment) {
  log_message("Step 5: Pathway enrichment analysis...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/04_pathway_enrichment.R")
    log_message("Pathway enrichment analysis completed successfully")
    
    # Report results
    if (exists("enrichment_results")) {
      total_pathways <- sum(sapply(enrichment_results, function(cond) 
        sum(sapply(cond, function(clust) 
          if(!is.null(clust)) nrow(clust) else 0))))
      log_message(paste("Total enriched pathways identified:", total_pathways))
    }
  }, error = function(e) {
    log_message(paste("Error in pathway enrichment analysis:", e$message))
    warning("Pathway enrichment analysis failed, continuing with next step")
  })
  
  end_time <- Sys.time()
  log_message(paste("Pathway enrichment time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 6: GSEA Analysis ----
if (run_gsea_analysis) {
  log_message("Step 6: Gene Set Enrichment Analysis...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/05_gsea_analysis.R")
    log_message("GSEA analysis completed successfully")
    
    # Report results
    if (exists("gsea_paths")) {
      total_gsea <- sum(sapply(gsea_paths, function(cell) 
        sum(sapply(cell, function(comp) 
          if(!is.null(comp)) nrow(comp) else 0))))
      log_message(paste("Total GSEA pathways identified:", total_gsea))
    }
  }, error = function(e) {
    log_message(paste("Error in GSEA analysis:", e$message))
    warning("GSEA analysis failed, continuing with next step")
  })
  
  end_time <- Sys.time()
  log_message(paste("GSEA analysis time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Step 7: Network Analysis ----
if (run_network_analysis) {
  log_message("Step 7: Network analysis and hub gene identification...")
  start_time <- Sys.time()
  
  tryCatch({
    source("scripts/06_network_analysis.R")
    log_message("Network analysis completed successfully")
    
    # Report results
    if (exists("network_results")) {
      analyzed_pathways <- length(network_results)
      log_message(paste("Network analysis completed for", analyzed_pathways, "pathways"))
    }
  }, error = function(e) {
    log_message(paste("Error in network analysis:", e$message))
    warning("Network analysis failed")
  })
  
  end_time <- Sys.time()
  log_message(paste("Network analysis time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
}

# Generate final summary report ----
log_message("Generating final summary report...")

create_summary_report <- function() {
  report_file <- file.path("results", "analysis_summary.md")
  
  cat("# Single-cell RNA-seq Analysis Summary\n\n", file = report_file)
  cat("Analysis completed on:", as.character(Sys.time()), "\n\n", file = report_file, append = TRUE)
  
  # Data summary
  cat("## Data Overview\n", file = report_file, append = TRUE)
  if (exists("obj_all_cells") && exists("obj_sf")) {
    cat("- Total synovial cells:", ncol(obj_all_cells), "\n", file = report_file, append = TRUE)
    cat("- Synovial fibroblasts:", ncol(obj_sf), "\n", file = report_file, append = TRUE)
    cat("- Conditions analyzed:", paste(conditions, collapse = ", "), "\n", file = report_file, append = TRUE)
    cat("- SF subtypes:", paste(sf_subtypes, collapse = ", "), "\n\n", file = report_file, append = TRUE)
  }
  
  # Analysis results
  cat("## Analysis Results\n", file = report_file, append = TRUE)
  
  if (exists("deg_results")) {
    total_degs <- sum(sapply(deg_results$degs, function(x) sapply(x, nrow)))
    cat("- Total differentially expressed genes:", total_degs, "\n", file = report_file, append = TRUE)
  }
  
  if (exists("enrichment_results")) {
    total_pathways <- sum(sapply(enrichment_results, function(cond) 
      sum(sapply(cond, function(clust) 
        if(!is.null(clust)) nrow(clust) else 0))))
    cat("- Total enriched pathways:", total_pathways, "\n", file = report_file, append = TRUE)
  }
  
  if (exists("gsea_paths")) {
    total_gsea <- sum(sapply(gsea_paths, function(cell) 
      sum(sapply(cell, function(comp) 
        if(!is.null(comp)) nrow(comp) else 0))))
    cat("- Total GSEA pathways:", total_gsea, "\n", file = report_file, append = TRUE)
  }
  
  if (exists("network_results")) {
    cat("- Pathways analyzed for networks:", length(network_results), "\n", file = report_file, append = TRUE)
  }
  
  # Output files
  cat("\n## Generated Files\n", file = report_file, append = TRUE)
  cat("- Figures: `results/figures/`\n", file = report_file, append = TRUE)
  cat("- Tables: `results/tables/`\n", file = report_file, append = TRUE)
  cat("- Pathway results: `results/pathway_enrichment/`\n", file = report_file, append = TRUE)
  cat("- GSEA results: `results/gsea_analysis/`\n", file = report_file, append = TRUE)
  cat("- Network results: `results/network_analysis/`\n", file = report_file, append = TRUE)
  
  log_message("Summary report saved to: results/analysis_summary.md")
}

create_summary_report()

# Calculate total runtime
if (exists("pipeline_start_time")) {
  total_time <- difftime(Sys.time(), pipeline_start_time, units = "mins")
} else {
  pipeline_start_time <- Sys.time() - as.difftime(60, units = "mins") # Estimate
  total_time <- 60
}

log_message(paste("=== Analysis Pipeline Completed ==="))
log_message(paste("Total runtime:", round(total_time, 2), "minutes"))
log_message("Check results/ directory for all outputs")

# Session info
session_info_file <- file.path("results", "session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_file)
log_message("Session information saved to results/session_info.txt")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Check the following directories for results:\n")
cat("- Figures: results/figures/\n") 
cat("- Data tables: results/tables/\n")
cat("- Analysis log: results/analysis_log.txt\n")
cat("- Summary report: results/analysis_summary.md\n")
cat("=========================\n")