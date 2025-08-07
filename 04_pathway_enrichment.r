# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 04: Pathway Enrichment Analysis
# ==============================================================================

source("scripts/03_differential_expression.R")

# Pathway analysis functions ----
convert_gene_ids <- function(genes, from_type = "ENSEMBL", to_type = "SYMBOL", 
                            species = "Homo sapiens") {
  cat(sprintf("Converting %d genes from %s to %s\n", length(genes), from_type, to_type))
  
  org_db <- ifelse(species == "Homo sapiens", "org.Hs.eg.db", "org.Mm.eg.db")
  
  genes_df <- bitr(genes, fromType = from_type, toType = to_type, 
                  OrgDb = org_db, drop = TRUE)
  colnames(genes_df) <- c("gene", "converted_id")
  
  cat(sprintf("Successfully converted %d genes\n", nrow(genes_df)))
  
  result <- list(
    original_genes = genes,
    converted_genes = genes_df,
    conversion_rate = nrow(genes_df) / length(genes),
    failed_conversions = length(genes) - nrow(genes_df)
  )
  
  return(result)
}

get_pathway_genesets <- function(gene_type = "SYMBOL", species = "Homo sapiens") {
  cat("=== Retrieving pathway gene sets ===\n")
  
  gene_col <- ifelse(gene_type == "ENTREZID", "entrez_gene", "gene_symbol")
  
  # Hallmark pathways
  cat("Retrieving Hallmark pathways...\n")
  hallmark_df <- msigdbr(species = species, category = "H")
  enr_sets_hallmark <- hallmark_df %>% 
    dplyr::select(gs_name, !!sym(gene_col))
  
  # C2 pathways
  cat("Retrieving C2 pathways...\n")
  c2_df <- msigdbr(species = species, category = "C2")
  
  # REACTOME
  enr_sets_reactome <- c2_df %>% 
    dplyr::filter(gs_subcat == "CP:REACTOME") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  # WikiPathways
  enr_sets_wiki <- c2_df %>% 
    dplyr::filter(gs_subcat == "CP:WIKIPATHWAYS") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  # KEGG
  enr_sets_kegg <- c2_df %>% 
    dplyr::filter(gs_subcat == "CP:KEGG") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  # C5 Gene Ontology pathways
  cat("Retrieving C5 Gene Ontology pathways...\n")
  c5_df <- msigdbr(species = species, category = "C5")
  
  # GO Biological Process
  enr_sets_gobp <- c5_df %>% 
    dplyr::filter(gs_subcat == "GO:BP") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  # GO Cellular Component
  enr_sets_gocc <- c5_df %>% 
    dplyr::filter(gs_subcat == "GO:CC") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  # GO Molecular Function
  enr_sets_gomf <- c5_df %>% 
    dplyr::filter(gs_subcat == "GO:MF") %>%
    dplyr::select(gs_name, !!sym(gene_col))
  
  pathway_sets <- list(
    hallmark = enr_sets_hallmark,
    kegg = enr_sets_kegg,
    reactome = enr_sets_reactome,
    wiki = enr_sets_wiki,
    gocc = enr_sets_gocc,
    gomf = enr_sets_gomf,
    gobp = enr_sets_gobp
  )
  
  # Create summary
  pathway_summary <- data.frame(
    pathway_type = names(pathway_sets),
    gene_sets_count = sapply(pathway_sets, function(x) length(unique(x$gs_name))),
    total_genes = sapply(pathway_sets, nrow)
  )
  
  cat("Pathway sets summary:\n")
  print(pathway_summary)
  
  result <- list(
    pathway_sets = pathway_sets,
    summary = pathway_summary,
    gene_type_used = gene_type,
    species_used = species
  )
  
  cat("✓ Pathway gene sets retrieval completed\n\n")
  return(result)
}

extract_first_number <- function(ratio) {
  as.numeric(unlist(str_split(ratio, "/"))[1])
}

process_pathway_names <- function(df) {
  cat("Processing pathway names and descriptions...\n")
  
  original_rows <- nrow(df)
  setDT(df)
  
  df$pathway <- df$Description
  df[, pathway := gsub("_", " ", pathway)]
  df$pathway <- sub(" ", ":", df$pathway, fixed = TRUE)
  df$type <- df$pathway
  df$pathway <- sub("^[^:]*:", "", df$pathway)
  
  df <- df %>%
    mutate(
      pathway = tolower(pathway),
      pathway = str_replace(pathway, "^(\\w)", function(x) toupper(x))
    )
  
  df$type <- sub(":.*$", "", df$type)
  
  # Reorder columns
  df <- df %>% 
    dplyr::select(type, pathway, p.adjust, geneID, Count, 
                 GeneRatio, genesets_num, everything())
  
  cat(sprintf("Processed %d pathway names\n", original_rows))
  return(df)
}

perform_enrichment <- function(genes_df, pathway_sets, pvalue_cutoff = 0.05) {
  cat("=== Performing enrichment analysis ===\n")
  
  gc() # Garbage collection
  enr_results <- list()
  
  for (pathway_name in names(pathway_sets)) {
    cat(sprintf("Processing pathway set: %s\n", pathway_name))
    
    enr_results[[pathway_name]] <- enricher(
      genes_df$converted_id,
      TERM2GENE = pathway_sets[[pathway_name]],
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
    
    # Add gene set numbers and ratios
    if (nrow(enr_results[[pathway_name]]@result) > 0) {
      enr_results[[pathway_name]]@result$genesets_num <-
        sapply(enr_results[[pathway_name]]@result$BgRatio, extract_first_number)
      
      enr_results[[pathway_name]]@result$GeneRatio <-
        enr_results[[pathway_name]]@result$Count /
        enr_results[[pathway_name]]@result$genesets_num
      
      significant_count <- nrow(enr_results[[pathway_name]]@result)
    } else {
      significant_count <- 0
    }
    
    cat(sprintf("  - Significant terms found: %d\n", significant_count))
  }
  
  cat("✓ Enrichment analysis completed\n\n")
  return(enr_results)
}

process_enrichment_results <- function(enr_results, cluster_id, condition, 
                                      output_dir = "enrichment_results") {
  cat(sprintf("=== Processing results for %s - %s ===\n", condition, cluster_id))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }
  
  # Combine all pathway results
  all_paths <- data.frame()
  
  for (pathway_name in names(enr_results)) {
    if (nrow(enr_results[[pathway_name]]@result) > 0) {
      pathway_results <- enr_results[[pathway_name]]@result
      pathway_results$pathway_source <- pathway_name
      all_paths <- rbind(all_paths, pathway_results)
    }
  }
  
  if (nrow(all_paths) > 0) {
    # Process pathway names and descriptions
    all_paths <- process_pathway_names(all_paths)
    
    # Filter significant pathways
    significant_paths <- all_paths %>% 
      filter(p.adjust <= 0.01) %>%
      arrange(p.adjust)
    
    # Save results
    filename_base <- sprintf("%s_%s", condition, cluster_id)
    
    sig_file <- file.path(output_dir, paste0("paths_sig_", filename_base, ".csv"))
    all_file <- file.path(output_dir, paste0("paths_all_", filename_base, ".csv"))
    
    fwrite(significant_paths, sig_file, row.names = FALSE)
    fwrite(all_paths, all_file, row.names = FALSE)
    
    # Print summary
    cat(sprintf("Total pathways tested: %d\n", nrow(all_paths)))
    cat(sprintf("Significant pathways (p.adj <= 0.01): %d\n", nrow(significant_paths)))
    cat(sprintf("Results saved to:\n  - %s\n  - %s\n", sig_file, all_file))
    
    # Get top 10 pathways
    top_pathways <- head(significant_paths$pathway, 10)
    
    if (length(top_pathways) > 0) {
      cat("Top 10 significant pathways:\n")
      for (i in 1:length(top_pathways)) {
        cat(sprintf("  %d. %s\n", i, top_pathways[i]))
      }
    }
    
    result <- list(
      significant_paths = significant_paths,
      all_paths = all_paths
    )
    
    cat("✓ Results processing completed\n\n")
    return(result)
  } else {
    cat("No enrichment results found for this cluster\n\n")
    return(NULL)
  }
}

# Main pathway enrichment functions ----
process_enrdegs <- function(obj, condition) {
  cat(sprintf("=== STEP 1: Processing enrdegs: %s ===\n", condition))
  
  # Subset and prepare data for current condition
  obj_used <- subset(obj, disease_state %in% condition)
  obj_used <- PrepSCTFindMarkers(obj_used)
  
  cat(sprintf("Cells in %s: %d\n", condition, ncol(obj_used)))
  
  # Find all markers
  cat("Finding differentially expressed genes...\n")
  
  all.markers <- FindAllMarkers(object = obj_used, only.pos = TRUE)
  all.markers$gene_name <- bm5$hgnc_symbol[match(all.markers$gene, bm5$ensembl_gene_id)]
  
  cat(sprintf("Total markers found: %d\n", nrow(all.markers)))
  cat(sprintf("Unique clusters: %d\n", length(unique(all.markers$cluster))))
  cat(sprintf("✓ Condition %s processing completed\n\n", condition))
  
  alldegs <- list()
  sigdegs <- list()
  
  # Process significant DEGs
  for (cluster in levels(all.markers$cluster)) {
    cat(sprintf("Processing cluster: %s\n", cluster))
    
    cluster_data <- all.markers %>%
      filter(cluster == !!cluster) %>%
      arrange(p_val_adj)
    
    significant_genes <- cluster_data %>%
      filter(p_val_adj < 0.05)
    
    alldegs[[cluster]] <- cluster_data
    sigdegs[[cluster]] <- significant_genes
    
    cat(sprintf("  - Total genes: %d, Significant: %d\n",
                nrow(cluster_data), nrow(significant_genes)))
  }
  
  result <- list(
    condition = condition,
    all_degs = alldegs,
    sig_degs = sigdegs
  )
  cat(sprintf("✓ Cluster processing completed for %s\n\n", condition))
  return(result)
}

process_enr_paths <- function(enr_degs, condition) {
  cat(sprintf("=== STEP 2: Pathway enrichment for %s ===\n", condition))
  enr_results <- list()
  pathway_data <- get_pathway_genesets()
  
  for (cluster_name in names(enr_degs$sig_degs)) {
    genes <- enr_degs$sig_degs[[cluster_name]]$gene
    gene_count <- length(genes)
    
    cat(sprintf("\nCluster %s: %d genes for enrichment\n", cluster_name, gene_count))
    
    if (gene_count > 0) {
      # Convert gene IDs
      gene_conversion <- convert_gene_ids(genes)
      
      if (nrow(gene_conversion$converted_genes) > 0) {
        # Perform enrichment
        enrichment_data <- perform_enrichment(
          gene_conversion$converted_genes, 
          pathway_data$pathway_sets
        )
        
        # Process and save results
        processed_results <- process_enrichment_results(
          enrichment_data, 
          cluster_name,
          condition
        )
        
        enr_results[[cluster_name]] <- list(
          enrichment_res = enrichment_data,
          processed_dataframe = processed_results
        )
      } else {
        cat(sprintf("No genes could be converted for cluster %s\n", cluster_name))
        enr_results[[cluster_name]] <- NULL
      }
    } else {
      cat(sprintf("No genes found for cluster %s\n", cluster_name))
      enr_results[[cluster_name]] <- NULL
    }
  }
  
  result <- list()
  result[[condition]] <- enr_results
  
  cat(sprintf("✓ Pathway enrichment completed for %s\n\n", condition))
  return(result)
}

# Execute pathway enrichment analysis ----
run_pathway_enrichment <- function(obj, conditions, subcategories, output_dir) {
  cat("Starting pathway enrichment analysis...\n")
  
  # Set up identifiers
  Idents(obj) <- "harmony_annotation"
  
  enrpaths_combine <- list()
  
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    cat(sprintf("\n=== Processing condition: %s ===\n", condition))
    
    # Process DEGs and enrichment
    enr_degs <- process_enrdegs(obj, condition)
    enr_paths <- process_enr_paths(enr_degs, condition)
    
    # Save intermediate results
    filename <- paste0(condition, "_Deg_enr.RData")
    filepath <- file.path(output_dir, filename)
    save(enr_degs, enr_paths, file = filepath)
    cat("Saved:", filepath, "\n")
    
    # Combine results
    for (subcat in subcategories) {
      if (!is.null(enr_paths[[condition]][[subcat]]) && 
          !is.null(enr_paths[[condition]][[subcat]][["processed_dataframe"]])) {
        enrpaths_combine[[condition]][[subcat]] <- 
          enr_paths[[condition]][[subcat]][["processed_dataframe"]][["significant_paths"]]
      }
    }
  }
  
  # Save combined results
  save(enrpaths_combine, file = file.path(output_dir, "enr_dotpltdata.RData"))
  cat("Combined enrichment results saved!\n")
  
  return(enrpaths_combine)
}

# Main execution ----
enrichment_output_dir <- file.path(output_dir, "pathway_enrichment")
dir.create(enrichment_output_dir, recursive = TRUE, showWarnings = FALSE)

# Run pathway enrichment analysis
enrichment_results <- run_pathway_enrichment(
  obj_sf, 
  conditions, 
  sf_subtypes, 
  enrichment_output_dir
)

cat("Pathway enrichment analysis completed!\n")
cat("Results saved to:", enrichment_output_dir, "\n")