# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 01: Data Loading and Preprocessing
# ==============================================================================

source("scripts/00_setup.R")

# Load all synovial cells data ----
cat("Loading all synovial cells data...\n")
load(file.path(base_path, "Documents/1_synovial_tissue/obj_no815.RData"))

# Check data structure
cat("Cohort distribution:\n")
print(table(obj_no815$cohort))
cat("\nSample distribution:\n") 
print(table(obj_no815$sample))
cat("\nCell type annotation:\n")
print(table(obj_no815$harmony_annotation))

# Store all cells object
obj_all_cells <- obj_no815
rm(obj_no815)
gc() # Garbage collection

# Load synovial fibroblasts subset ----
cat("\nLoading synovial fibroblasts subset...\n")
load(file.path(base_path, "Documents/1_synovial_tissue/0.cell_type_SFsall/dis_nor/SFs_4cohorts.RData"))

cat("SF cohort distribution:\n")
print(table(obj_SFs$cohort))
cat("\nSF sample distribution:\n")
print(table(obj_SFs$sample))  
cat("\nSF subtype annotation:\n")
print(table(obj_SFs$harmony_annotation))

# Store SF object
obj_sf <- obj_SFs
rm(obj_SFs)
gc()

# Data preprocessing functions ----
prepare_sf_annotations <- function(obj) {
  # Create binary annotations
  obj$harmony_annotation1 <- ifelse(obj$harmony_annotation == "SLSFs_CD34pos", 
                                   "SFs_CD34pos", "SFs_CD34neg")
  obj$harmony_annotation2 <- ifelse(obj$harmony_annotation == "LLSFs", 
                                   "LLSFs", "SLSFs")
  return(obj)
}

create_comparison_groups <- function(obj, cluster_col, condition_col) {
  obj$compare <- paste(obj@meta.data[[cluster_col]], 
                      obj@meta.data[[condition_col]], sep = "_")
  return(obj)
}

# Apply preprocessing to SF data ----
obj_sf <- prepare_sf_annotations(obj_sf)
obj_sf <- create_comparison_groups(obj_sf, "harmony_annotation", "disease_state")

# Prepare SCT for differential expression ----
cat("\nPreparing SCT assay for FindMarkers...\n")
obj_sf_presct <- PrepSCTFindMarkers(obj_sf)

# Export count matrix for cell deconvolution ----
cat("Preparing count matrix for cell deconvolution...\n")
count_matrix <- as.matrix(obj_sf@assays$SCT@counts)
cell_types <- obj_sf$harmony_annotation
colnames(count_matrix) <- cell_types
rownames(count_matrix) <- bm5$hgnc_symbol[match(rownames(count_matrix), bm5$ensembl_gene_id)]

# Save count matrix
save(count_matrix, file = file.path(output_dir, 'scSFs_matrix.RData'))

cat("\nData loading and preprocessing complete!\n")
cat("All cells object: obj_all_cells -", ncol(obj_all_cells), "cells\n")
cat("SF cells object: obj_sf -", ncol(obj_sf), "cells\n")
cat("SF subtypes:", paste(unique(obj_sf$harmony_annotation), collapse = ", "), "\n")
cat("Count matrix saved for cell deconvolution\n")