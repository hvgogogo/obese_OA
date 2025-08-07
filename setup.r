# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 00: Environment Setup and Configuration
# ==============================================================================

# Install and load required packages ----
packages <- c(
  # Core analysis
  "Seurat", "dplyr", "tidyverse", "data.table",
  
  # Biomart and gene annotation
  "biomaRt", "org.Hs.eg.db", "stringr",
  
  # Visualization
  "ggplot2", "patchwork", "cowplot", "ggrepel", "pheatmap", "Cairo",
  "Nebulosa", "ggpubr", "forcats", "aplot",
  
  # Pathway analysis
  "clusterProfiler", "msigdbr", "fgsea", "enrichplot",
  
  # Network analysis
  "igraph", "STRINGdb", "ggraph", "visNetwork",
  
  # Statistical analysis
  "reshape2"
)

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# Set paths ----
base_path <- "~/Library/CloudStorage/OneDrive-CardiffUniversity"
fig_path <- "/Users/hangwei/Library/CloudStorage/OneDrive-CardiffUniversity/Documents/0.2_reserch_progress_cu/8_thesis2025/Thesis_WeiHang/C5-Figs/R"

# Create output directories
output_dir <- fig_path
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "enrichment"), recursive = TRUE, showWarnings = FALSE)

# Load gene annotation data
load(file.path(base_path, "localR/bm5.RData"))

# Define global plotting theme ----
plot_theme <- theme(
  text = element_text(size = 10, face = 'bold'),
  title = element_blank(),
  plot.title = element_blank(),
  axis.text = element_blank(),      
  axis.ticks = element_blank(),      
  axis.line = element_blank(),       
  legend.text = element_text(size = 10, face = "bold"),
  legend.title = element_text(size = 12, face = "bold"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.background = element_blank(),
  legend.key.size = unit(0.5, "cm")
)

# Define global variables ----
conditions <- c("lean_healthy", "lean_OA", "obese_OA")
cell_types_all <- c("SFs", "Endos", "Adipocytes", "SMCs", 
                    "Macrophages", "Masts", "NK-T", "B")
sf_subtypes <- c("SLSFs_CD34pos", "SLSFs_CD34neg", "LLSFs", "SLSFs_IM")

# Define color palettes ----
colors_general <- c("#9C6BA3", "#FE9C9D", "#98C897", "#9DBAD2", "#F8BC7E", "#CC976B")
color_used_sf <- c(
  "SLSFs_CD34pos" = "#6ABA88",
  "LLSFs" = "#C576F6",
  "SLSFs_CD34neg" = "#366BA1",
  "SLSFs_IM" = "#D7E6F5"
)

# Analysis parameters ----
avg_log2FC_cutoff <- 1
p_val_adj_cutoff <- 0.05
string_score_threshold <- 400

cat("Environment setup complete!\n")
cat("Base path:", base_path, "\n")
cat("Output directory:", output_dir, "\n")
cat("Loaded", length(packages), "packages successfully.\n")