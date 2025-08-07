# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 03: Differential Expression Analysis
# ==============================================================================

source("scripts/02_cell_characterization.R")

# Helper functions ----
get_condition_pairs <- function(comparison_type, unique_conditions) {
  switch(comparison_type,
    "weight" = list(group1 = unique_conditions[1:2], group2 = unique_conditions[3]),
    "OA" = list(group1 = unique_conditions[1], group2 = unique_conditions[2:3]),
    "obeseOA" = list(group1 = unique_conditions[1], group2 = unique_conditions[3])
  )
}

get_cell_mapping <- function(cell_type, available_cells) {
  cell_map <- list(
    "LLSFs" = available_cells[1], 
    "SLSFs_CD34pos" = available_cells[3],
    "SLSFs_CD34neg" = available_cells[2], 
    "SLSFs_IM" = available_cells[4]
  )
  return(cell_map[[cell_type]])
}

process_degs <- function(marker_results, fc_cutoff = 1, p_cutoff = 0.05, removed_genes = c()) {
  clean_results <- marker_results[!marker_results$gene %in% removed_genes, ] %>% 
    arrange(desc(p_val_adj))
  clean_results$geneid <- rownames(clean_results)
  
  up_genes <- clean_results %>% 
    filter(avg_log2FC > fc_cutoff & p_val_adj < p_cutoff)
  down_genes <- clean_results %>% 
    filter(avg_log2FC < -fc_cutoff & p_val_adj < p_cutoff)
  all_degs <- rbind(up_genes, down_genes)
  
  # Get top features for plotting
  features <- c(
    head(up_genes$gene, 5)[!is.na(head(up_genes$gene, 5))],
    head(down_genes$gene, 5)[!is.na(head(down_genes$gene, 5))]
  )
  features <- features[!is.na(features)]
  
  return(list(
    all_genes = clean_results, 
    degs = all_degs, 
    up_genes = up_genes,
    down_genes = down_genes, 
    features = features
  ))
}

# Visualization functions ----
volcano_plot <- function(df, selectgenes, fc_cutoff = avg_log2FC_cutoff, p_cutoff = p_val_adj_cutoff) {
  up_genes_sig <- df$geneid[df$avg_log2FC > fc_cutoff & df$p_val_adj < p_cutoff]
  dw_genes_sig <- df$geneid[df$avg_log2FC < -fc_cutoff & df$p_val_adj < p_cutoff]
  genes_notsig <- df$geneid[df$p_val_adj >= p_cutoff]
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj),
                color = ifelse(geneid %in% genes_notsig, "notsig",
                              ifelse(geneid %in% up_genes_sig, "UP-Genes",
                                    ifelse(geneid %in% dw_genes_sig, "Down-Genes", "sig"))))) +
    geom_point(data = subset(df, !geneid %in% c(up_genes_sig, dw_genes_sig)), 
               size = 0.5, alpha = 1) +
    geom_point(data = subset(df, geneid %in% c(up_genes_sig, dw_genes_sig)), 
               size = 2, alpha = 0.7) +
    scale_color_manual(values = c("UP-Genes" = "#E41A1C", "Down-Genes" = "#377EB8", 
                                 "sig" = "gray40", "notsig" = "gray40"),
                      name = NULL, breaks = c("Down-Genes", "UP-Genes")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_text_repel(data = subset(df, gene %in% selectgenes), 
                    aes(label = gene), color = "black", fontface = "bold", 
                    size = 5, max.overlaps = 10) +
    labs(x = expression(paste("Log"[2], " fold change")), 
         y = expression("-Log"[10]*" p-value")) +
    theme_minimal() + 
    theme(legend.position = "top", panel.grid = element_blank())
}

waterfall_plot <- function(data, top_genes, tail_genes) {
  data <- data %>% 
    mutate(color = ifelse(rank %in% c(top_genes$rank, tail_genes$rank), "red", "grey"))
  
  ggplot(data) +
    geom_point(aes(x = rank, y = avg_log2FC, color = color, size = abs(avg_log2FC))) +
    scale_color_manual(values = c("red" = "red", "grey" = "grey")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_text_repel(data = top_genes, aes(x = rank + 10, y = avg_log2FC, label = gene),
                    box.padding = 0.5, nudge_x = 10, nudge_y = 0.2, 
                    hjust = "left", max.overlaps = 15) +
    geom_text_repel(data = tail_genes, aes(x = rank + 10, y = avg_log2FC, label = gene),
                    box.padding = 0.5, nudge_x = 10, nudge_y = -0.2, 
                    hjust = "right", max.overlaps = 15) +
    scale_size(range = c(1, 5)) + 
    guides(color = "none") +
    labs(x = "Rank", y = "log2FoldChange") + 
    theme_bw() + 
    theme(panel.grid = element_blank())
}

create_plots <- function(gene_results, comparison_name, method, output_path) {
  # Volcano plot
  p_volcano <- volcano_plot(gene_results$all_genes, gene_results$features) + 
    labs(title = method)
  cairo_pdf(file.path(output_path, paste0("volcano_", comparison_name, ".pdf")), 
            width = 5, height = 6)
  print(p_volcano)
  dev.off()
  
  # Waterfall plot
  if (nrow(gene_results$degs) > 0) {
    deg_ranked <- gene_results$degs %>% 
      mutate(rank = rank(avg_log2FC, ties.method = "max"))
    top_genes <- deg_ranked %>% 
      filter(!is.na(gene)) %>% 
      arrange(desc(avg_log2FC)) %>% 
      slice_head(n = 10)
    bottom_genes <- deg_ranked %>% 
      filter(!is.na(gene)) %>% 
      arrange(desc(avg_log2FC)) %>% 
      slice_tail(n = 10)
    
    p_waterfall <- waterfall_plot(deg_ranked, top_genes, bottom_genes) + 
      labs(title = method)
    cairo_pdf(file.path(output_path, paste0("waterfall_", comparison_name, ".pdf")), 
              width = 7.2, height = 9)
    print(p_waterfall)
    dev.off()
  }
}

save_results <- function(gene_results, comparison_name, output_path) {
  write.csv(gene_results$all_genes, 
            file.path(output_path, paste0("allGenes_", comparison_name, ".csv")), 
            row.names = FALSE)
  write.csv(gene_results$degs, 
            file.path(output_path, paste0("DEGs_", comparison_name, ".csv")), 
            row.names = FALSE)
}

# Main differential expression analysis ----
run_deg_analysis <- function(obj, output_dir) {
  cat("Starting differential expression analysis...\n")
  
  # Setup parameters
  methods_used <- "wilcox"
  based_clusters <- "harmony_annotation"
  condition_clu <- "disease_state"
  condi2compares <- c("weight", "OA", "obeseOA")
  condis <- c("LLSFs", "SLSFs_CD34pos", "SLSFs_CD34neg", "SLSFs_IM")
  removed_degs <- c()
  
  # Initialize storage
  allgenes_list <- degs_list <- list()
  
  # Get available cells and conditions
  icells <- sort(unique(obj$harmony_annotation))
  unique_conditions <- unique(obj@meta.data[[condition_clu]])
  
  # Create comparison groups
  obj$compare <- paste(obj@meta.data[[based_clusters]], 
                      obj@meta.data[[condition_clu]], sep = "_")
  Idents(obj) <- factor(obj$compare, levels = sort(levels(factor(obj$compare))))
  
  # Run analysis for each comparison type
  for (comparison_type in condi2compares) {
    cat("\nAnalyzing:", comparison_type, "\n")
    condition_pairs <- get_condition_pairs(comparison_type, unique_conditions)
    cell_condi1 <- condition_pairs$group1
    cell_condi2 <- condition_pairs$group2
    
    for (cell_type in condis) {
      cat("Processing:", cell_type, "\n")
      target_cell <- get_cell_mapping(cell_type, icells)
      
      ident1 <- paste(target_cell, cell_condi1, sep = "_")
      ident2 <- paste(target_cell, cell_condi2, sep = "_")
      
      clustering <- Idents(obj)
      cells1 <- names(clustering[clustering %in% ident1])
      cells2 <- names(clustering[clustering %in% ident2])
      
      if (length(cells1) == 0 || length(cells2) == 0) next
      
      tryCatch({
        # Find markers
        marker_results <- FindMarkers(obj, cells2, cells1, 
                                     logfc.threshold = 0, min.pct = 0, 
                                     test.use = methods_used) %>% 
          arrange(desc(avg_log2FC))
        marker_results$gene <- bm5$hgnc_symbol[match(rownames(marker_results), 
                                                    bm5$ensembl_gene_id)]
        
        # Process results
        deg_analysis <- process_degs(marker_results, avg_log2FC_cutoff, 
                                    p_val_adj_cutoff, removed_degs)
        comparison_name <- paste(cell_type, comparison_type, sep = "_")
        
        # Create plots and save results
        create_plots(deg_analysis, comparison_name, methods_used, output_dir)
        save_results(deg_analysis, comparison_name, output_dir)
        
        # Store results
        degs_list[[cell_type]][[comparison_type]] <- deg_analysis$degs
        allgenes_list[[cell_type]][[comparison_type]] <- deg_analysis$all_genes
        
        cat("DEGs:", nrow(deg_analysis$degs), "| Up:", nrow(deg_analysis$up_genes),
            "| Down:", nrow(deg_analysis$down_genes), "\n")
            
      }, error = function(e) cat("Error:", e$message, "\n"))
    }
  }
  
  # Save all results
  save(degs_list, allgenes_list, 
       file = file.path(output_dir, "differential_expression_results.RData"))
  
  cat("Differential expression analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(list(degs = degs_list, all_genes = allgenes_list))
}

# Execute analysis ----
deg_output_dir <- file.path(output_dir, "differential_expression")
dir.create(deg_output_dir, recursive = TRUE, showWarnings = FALSE)

# Run the analysis
deg_results <- run_deg_analysis(obj_sf_presct, deg_output_dir)

cat("Differential expression analysis pipeline completed!\n")