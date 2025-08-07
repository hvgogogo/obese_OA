# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 02: Cell Type Characterization and Visualization
# ==============================================================================

source("scripts/01_data_loading.R")

# Figure 1A: All Cells UMAP ----
create_all_cells_umap <- function(obj, reduction = "harmony_redu_sct_mt_tech") {
  cat("Creating UMAP plots for all synovial cells...\n")
  
  Idents(obj) <- "harmony_annotation"
  
  # Cohort UMAP
  p1 <- DimPlot(obj, reduction = reduction, group.by = "cohort", 
                cols = colors_general[1:3]) +
    guides(color = guide_legend(title = "Cohorts"))
  
  # Disease state UMAP  
  p2 <- DimPlot(obj, reduction = reduction, group.by = "disease_state", 
                cols = colors_general[4:6]) +
    guides(color = guide_legend(title = "Disease State"))
  
  combined_plot <- (p2 / p1) & plot_theme
  ggsave(file.path(output_dir, "figures/allcells_umap.pdf"), 
         combined_plot, width = 4.5, height = 6)
  
  # Annotated UMAP
  p_annotated <- DimPlot(obj, reduction = reduction, group.by = "harmony_annotation",
                        label = TRUE, label.size = 8) + NoLegend()
  p_annotated <- p_annotated & plot_theme
  ggsave(file.path(output_dir, "figures/allcells_umap_annotation.pdf"), 
         p_annotated, width = 6, height = 6)
  
  return(list(cohort = p1, disease = p2, annotated = p_annotated))
}

# Create marker heatmap for all cells ----
create_marker_heatmap <- function(obj, validation_markers, cluster_col = "harmony_annotation") {
  cat("Creating marker gene heatmap...\n")
  
  desired_order <- c("SFs", "Endos", "Adipocytes", "SMCs",
                     "Macrophages", "Masts", "NK-T", "B")
  
  if (length(validation_markers) > 0) {
    tryCatch({
      # Calculate average expression
      expression_data <- GetAssayData(obj, layer = "data")
      rownames(expression_data) <- bm5$hgnc_symbol[match(
        rownames(expression_data), bm5$ensembl_gene_id)]
      expression_data <- expression_data[nzchar(rownames(expression_data), keepNA = FALSE), ]
      
      cluster_info <- obj@meta.data[[cluster_col]]
      names(cluster_info) <- colnames(obj)
      unique_clusters <- sort(unique(cluster_info))
      
      avg_exp_manual <- matrix(0, nrow = length(validation_markers),
                              ncol = length(unique_clusters))
      rownames(avg_exp_manual) <- validation_markers
      colnames(avg_exp_manual) <- unique_clusters
      
      for (marker in validation_markers) {
        if (marker %in% rownames(expression_data)) {
          marker_exp <- expression_data[marker, ]
          for (cluster in unique_clusters) {
            cells_in_cluster <- names(cluster_info)[cluster_info == cluster]
            cells_in_cluster <- cells_in_cluster[cells_in_cluster %in% colnames(expression_data)]
            if (length(cells_in_cluster) > 0) {
              avg_exp_manual[marker, as.character(cluster)] <- 
                mean(marker_exp[cells_in_cluster], na.rm = TRUE)
            }
          }
        }
      }
      
      avg_exp_manual <- avg_exp_manual[validation_markers, desired_order]
      
      # Create heatmap
      pdf(file.path(output_dir, "figures/allcells_markerheatmap.pdf"), 
          width = 6, height = 6)
      marker_pheatmap <- pheatmap(avg_exp_manual,
                                 cluster_cols = FALSE, cluster_rows = FALSE,
                                 scale = "row", 
                                 color = colorRampPalette(c("blue", "white", "red"))(100),
                                 fontsize = 16, cellwidth = 25, cellheight = 20,
                                 angle_col = 45, border_color = "white")
      print(marker_pheatmap)
      dev.off()
      
      cat("Average expression matrix:\n")
      print(round(avg_exp_manual, 2))
      
    }, error = function(e) {
      cat("Error creating heatmap:", e$message, "\n")
    })
  }
}

# Figure 1B: SF Cells UMAP ----
create_sf_umap <- function(obj, reduction = "harmony_redu_sct_mt_cohort") {
  cat("Creating UMAP plots for synovial fibroblasts...\n")
  
  # CD34 positive/negative
  p1 <- DimPlot(obj, reduction = reduction, group.by = "harmony_annotation1", 
                cols = c("#F8BC7E", "#98C897"))
  p1 <- p1 & plot_theme + 
    theme(legend.position = "inside", legend.position.inside = c(0.95, 0.85),
          legend.justification = c(1, 1))
  ggsave(file.path(output_dir, "figures/SFcells_CD34pos-neg.pdf"), 
         p1, width = 2, height = 2)
  
  # Lining/sublining
  p2 <- DimPlot(obj, reduction = reduction, group.by = "harmony_annotation2", 
                cols = c("#FE9C9D", "#9DBAD2"))
  p2 <- p2 & plot_theme + 
    theme(legend.position = "inside", legend.text = element_text(size = 16),
          legend.position.inside = c(1, 1.1), legend.justification = c(1, 1))
  ggsave(file.path(output_dir, "figures/SFcells_SL-LL.pdf"), 
         p2, width = 2, height = 2)
  
  # All subtypes annotated
  p3 <- DimPlot(obj, reduction = reduction, group.by = "harmony_annotation",
                cols = c("#C576F6", "#366BA1", "#6ABA88", "#D7E6F5"),
                label = TRUE, label.size = 8) + NoLegend()
  p3 <- p3 & plot_theme
  ggsave(file.path(output_dir, "figures/SFcells_umap_annotation.pdf"), 
         p3, width = 6, height = 6)
  
  return(list(cd34 = p1, lining = p2, annotated = p3))
}

# Create violin plots ----
create_violin_plots <- function(obj, genes, group_by = "harmony_annotation") {
  cat("Creating violin plots for key genes...\n")
  
  gene_ids <- bm5$ensembl_gene_id[match(genes, bm5$hgnc_symbol)]
  
  vln_plot <- function(obj, group_by) {
    p <- VlnPlot(obj, features = gene_ids, group.by = group_by, 
                assay = "SCT", layer = "data", stack = TRUE, 
                log = FALSE, flip = FALSE, pt.size = 0.5) + 
      theme_classic()
    p$data$feature <- bm5$hgnc_symbol[match(p$data$feature, bm5$ensembl_gene_id)]
    p$data$feature <- factor(p$data$feature, levels = genes)
    return(p)
  }
  
  p1 <- vln_plot(obj, group_by) & plot_theme +
    theme(axis.text.x.top = element_blank(), 
          axis.text.x.bottom = element_blank(),
          axis.text.y = element_text(size = 10), 
          axis.ticks.x.bottom = element_blank(),
          legend.position = "none")
  
  ggsave(file.path(output_dir, "figures/SFcells_violinplot.pdf"), 
         p1, width = 5, height = 5)
  
  return(p1)
}

# Create density plots ----
create_density_plots <- function(obj, features_list, 
                                reduction = "harmony_redu_sct_mt_cohort") {
  cat("Creating density plots...\n")
  
  library(Nebulosa)
  
  generate_density_plot <- function(obj, gene_symbol, joint_1 = FALSE) {
    gene_id <- bm5$ensembl_gene_id[match(gene_symbol, bm5$hgnc_symbol)]
    
    density_plots <- plot_density(
      obj = obj, features = gene_id, slot = "data",
      reduction = reduction, size = 0.4, pal = "magma", joint = joint_1
    )
    return(density_plots)
  }
  
  theme_densityplot <- theme(
    text = element_text(size = 10, face = 'bold'),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.title = element_blank(), axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "inside", legend.position.inside = c(0.9, 1),
    legend.justification = c(1, 1),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )
  
  for (feature in features_list) {
    p <- generate_density_plot(obj, feature, TRUE) &
      ggtitle(feature) & theme_densityplot
    
    filename <- paste0("SFcells_density_", feature)
    ggsave(file.path(output_dir, "figures", paste0(filename, ".pdf")),
           p, width = 3, height = 3)
  }
}

# Create pie charts ----
create_pie_charts <- function(obj) {
  cat("Creating pie charts for cell proportions...\n")
  
  cell_counts <- table(obj$harmony_annotation, obj$disease_state)
  data_long <- as.data.frame(cell_counts)
  colnames(data_long) <- c("harmony_annotation", "disease_state", "Count")
  data_long$harmony_annotation <- factor(data_long$harmony_annotation,
                                        levels = c("SLSFs_CD34pos", "SLSFs_CD34neg", 
                                                  "LLSFs", "SLSFs_IM"))
  
  create_pie_chart <- function(disease_name) {
    disease_data <- data_long %>%
      filter(disease_state == disease_name) %>%
      mutate(
        Percentage = Count / sum(Count) * 100,
        ymax = cumsum(Percentage),
        ymin = c(0, head(cumsum(Percentage), n = -1)),
        labelPosition = (ymax + ymin) / 2,
        label = paste0(round(Percentage, 1), "%")
      )
    
    ggplot(disease_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, 
                            fill = harmony_annotation)) +
      geom_rect(color = "white", size = 0.5) +
      coord_polar(theta = "y") + xlim(c(0, 4)) +
      theme_void() + scale_fill_manual(values = color_used_sf) +
      geom_text(aes(label = label, x = 3, y = labelPosition),
               size = 3.5, color = "black", fontface = "bold") +
      labs(title = gsub("_", " ", disease_name), fill = "Cell Types") +
      theme(plot.title = element_text(size = 12, hjust = 0.5, vjust = -34, face = "bold"),
            text = element_text(size = 10, face = "bold"),
            legend.position = "bottom")
  }
  
  p1 <- create_pie_chart("lean_healthy") & ggtitle("Healthy")
  p2 <- create_pie_chart("lean_OA") & ggtitle("non-obese OA")
  p3 <- create_pie_chart("obese_OA") & ggtitle("Obese OA")
  
  combined_plot <- p1 + p2 + p3 + 
    plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path(output_dir, "figures/SFs_pieplot.pdf"),
         combined_plot, width = 9, height = 3)
  
  return(combined_plot)
}

# Main execution ----
cat("Starting cell characterization analysis...\n")

# All cells analysis
validation_markers <- c("PDGFRA", "COL1A2", "VWF", "ACACB", "ACTA2", 
                        "CD14", "KIT", "CD2", "CD79A")
all_cells_plots <- create_all_cells_umap(obj_all_cells)
create_marker_heatmap(obj_all_cells, validation_markers)

# SF cells analysis  
sf_plots <- create_sf_umap(obj_sf)

# Key gene analysis
key_genes <- c("CYBB", "CLIC5", "THY1", "CD34")
violin_plots <- create_violin_plots(obj_sf, key_genes)

# Density plots
density_features <- c("CD34", "CLIC5", "CYBB", "THY1")  
create_density_plots(obj_sf, density_features)

# Pie charts
pie_plots <- create_pie_charts(obj_sf)

cat("Cell characterization analysis complete!\n")
cat("Generated plots saved in:", file.path(output_dir, "figures/"), "\n")