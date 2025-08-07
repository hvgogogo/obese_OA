# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Utility Functions: Plotting and Visualization
# ==============================================================================

# Advanced violin plot with statistical comparisons ----
violin_plot3 <- function(obj, features, ncols = 3, group_col = "Type", 
                        filter0 = FALSE, palette = c("#56B4E9", "#E69F00", "#E63232")) {
  
  # Get Ensembl IDs for the features
  features_id <- bm5$ensembl_gene_id[match(features, bm5$hgnc_symbol)]
  
  # Extract expression data for each feature
  expr_data_list <- lapply(seq_along(features_id), function(i) {
    data <- FetchData(obj, vars = c(group_col, features_id[i]))
    colnames(data)[1:2] <- c("Type", features[i])
    data
  })
  
  # Combine data
  type_column <- expr_data_list[[1]]["Type"]
  data <- do.call(cbind, lapply(expr_data_list, function(df) df[, -which(names(df) == "Type")]))
  data <- as.data.frame(data)
  colnames(data) <- features
  data <- cbind(type_column, data)
  
  # Melt data for plotting
  data <- reshape2::melt(data, id.vars = c("Type"))
  colnames(data) <- c("Type", "Gene", "Expression")
  
  # Filter zero values if requested
  if (filter0 == TRUE) {
    data <- data %>% filter(Expression != 0)
  }
  
  # Prepare comparisons
  factors <- levels(factor(data$Type))
  compare_factor <- combn(factors, 2, simplify = FALSE)
  
  # Create violin plot
  p <- ggviolin(data, x = "Type", y = "Expression", fill = "Type",
                ylab = "Gene Expression", xlab = "Gene") +
    scale_fill_manual(values = palette) +
    stat_compare_means(
      comparisons = compare_factor,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      label = "p.signif",
      vjust = -0.0
    ) +
    facet_wrap(~Gene, ncol = ncols, nrow = ceiling(length(features)/ncols)) +
    theme_bw() +
    theme(
      text = element_text(size = 10, face = 'bold'),
      title = element_text(size = 0, face = 'bold'),
      axis.line = element_line(linewidth = 0.5),
      axis.title.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      axis.title.y.left = element_text(size = 12, face = 'bold'),
      axis.text.x.bottom = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      panel.grid = element_blank()
    )
  
  return(p)
}

# Enhanced DotPlot function ----
create_enhanced_dotplot <- function(obj, features_id, features, group_by) {
  DotPlot(object = obj,
          features = features_id,
          cols = c('blue', 'red'),
          group.by = group_by) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0)) +
    scale_x_discrete(labels = features) +
    coord_flip()
}

# Create publication-ready UMAP with custom annotations ----
create_publication_umap <- function(obj, reduction, group_by, colors = NULL, 
                                   label = TRUE, label_size = 8) {
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(min(8, length(unique(obj@meta.data[[group_by]]))), "Set2")
  }
  
  # Create base UMAP
  df <- data.frame(
    UMAP_1 = obj@reductions[[reduction]]@cell.embeddings[, 1],
    UMAP_2 = obj@reductions[[reduction]]@cell.embeddings[, 2],
    cluster = obj@meta.data[[group_by]]
  ) %>%
    mutate(center_x = ave(UMAP_1, cluster, FUN = median),
           center_y = ave(UMAP_2, cluster, FUN = median))
  
  p <- ggplot(df, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.2) +
    scale_color_manual(values = colors) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_void() +
    guides(color = "none")
  
  if (label) {
    p <- p +
      geom_point(aes(center_x, center_y), size = 10, color = "grey90", alpha = 0.6) +
      geom_text(aes(center_x, center_y, label = cluster), 
                color = "black", size = label_size, fontface = "bold")
  }
  
  return(p)
}

# Enhanced heatmap function for pathway modules ----
create_pathway_heatmap <- function(obj, pathway_genes, group_by, 
                                  pathway_names = NULL, colors = NULL) {
  
  if (is.null(colors)) {
    colors <- colorRampPalette(c("#003366", "white", "#990033"))(100)
  }
  
  # Convert gene symbols to Ensembl IDs
  features_id <- bm5$ensembl_gene_id[match(pathway_genes, bm5$hgnc_symbol)]
  features_id <- features_id[!is.na(features_id)]
  
  # Create dotplot to get average expression data
  plot_data <- create_enhanced_dotplot(obj, features_id, pathway_genes, group_by)
  marker_data <- plot_data$data
  marker_data$Zscore <- marker_data$avg.exp.scaled
  marker_data$cluster <- marker_data$id
  
  # Create heatmap
  p <- ggplot(marker_data, aes(y = features.plot, x = cluster, fill = Zscore)) +
    geom_raster() +
    scale_fill_gradient2(low = colors[1], high = colors[length(colors)], 
                        mid = colors[length(colors)/2]) +
    labs(x = NULL, y = NULL, fill = "Scaled\nExpression") +
    theme_minimal() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.position = "right"
    )
  
  # Convert Ensembl IDs back to gene symbols for labels
  levels(p$data$features.plot) <- bm5$hgnc_symbol[match(levels(p$data$features.plot), 
                                                        bm5$ensembl_gene_id)]
  
  return(p)
}

# Create waterfall plot for differential expression ----
create_waterfall_plot <- function(data, top_genes, tail_genes, title = "") {
  data <- data %>% 
    mutate(color = ifelse(rank %in% c(top_genes$rank, tail_genes$rank), "red", "grey"))
  
  ggplot(data, aes(x = rank, y = avg_log2FC)) +
    geom_point(aes(color = color, size = abs(avg_log2FC)), alpha = 0.7) +
    scale_color_manual(values = c("red" = "#E41A1C", "grey" = "grey70")) +
    scale_size_continuous(range = c(0.5, 3), name = "|log2FC|") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_text_repel(data = top_genes, 
                    aes(label = gene), 
                    nudge_y = 0.5, size = 3, fontface = "bold",
                    max.overlaps = 15) +
    geom_text_repel(data = tail_genes, 
                    aes(label = gene), 
                    nudge_y = -0.5, size = 3, fontface = "bold",
                    max.overlaps = 15) +
    labs(x = "Gene Rank", y = "log2 Fold Change", title = title) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(color = "none")
}

# Enhanced volcano plot ----
create_enhanced_volcano <- function(df, selectgenes, fc_cutoff = 1, p_cutoff = 0.05, 
                                   title = "") {
  # Categorize genes
  up_genes_sig <- df$geneid[df$avg_log2FC > fc_cutoff & df$p_val_adj < p_cutoff]
  dw_genes_sig <- df$geneid[df$avg_log2FC < -fc_cutoff & df$p_val_adj < p_cutoff]
  genes_notsig <- df$geneid[df$p_val_adj >= p_cutoff]
  
  # Create color mapping
  df$gene_category <- case_when(
    df$geneid %in% up_genes_sig ~ "Upregulated",
    df$geneid %in% dw_genes_sig ~ "Downregulated", 
    df$geneid %in% genes_notsig ~ "Not significant",
    TRUE ~ "Significant"
  )
  
  # Create plot
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = gene_category)) +
    geom_point(aes(size = ifelse(gene_category %in% c("Upregulated", "Downregulated"), 1.5, 0.8)),
               alpha = 0.7) +
    scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", 
                                 "Not significant" = "grey70", "Significant" = "grey50"),
                      name = "Gene Status") +
    scale_size_identity() +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), 
               linetype = "dashed", color = "grey50", alpha = 0.8) +
    geom_hline(yintercept = -log10(p_cutoff), 
               linetype = "dashed", color = "grey50", alpha = 0.8) +
    geom_text_repel(data = subset(df, gene %in% selectgenes), 
                    aes(label = gene), color = "black", fontface = "bold",
                    size = 4, max.overlaps = 15, box.padding = 0.5) +
    labs(x = expression(paste("Log"[2], " Fold Change")), 
         y = expression("-Log"[10], " Adjusted p-value"),
         title = title) +
    theme_minimal() +
    theme(
      legend.position = "top",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# Create publication-ready pie charts ----
create_publication_pie <- function(data, group_col, fill_col, colors = NULL, title = "") {
  
  # Calculate percentages
  pie_data <- data %>%
    count(!!sym(group_col), !!sym(fill_col)) %>%
    group_by(!!sym(group_col)) %>%
    mutate(
      percentage = n / sum(n) * 100,
      ymax = cumsum(percentage),
      ymin = c(0, head(cumsum(percentage), n = -1)),
      labelPosition = (ymax + ymin) / 2,
      label = paste0(round(percentage, 1), "%")
    )
  
  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(length(unique(pie_data[[fill_col]])), "Set3")
  }
  
  ggplot(pie_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, 
                      fill = !!sym(fill_col))) +
    geom_rect(color = "white", size = 0.8) +
    coord_polar(theta = "y") +
    xlim(c(0, 4)) +
    scale_fill_manual(values = colors) +
    geom_text(aes(label = label, x = 3, y = labelPosition),
              size = 4, color = "black", fontface = "bold") +
    facet_wrap(as.formula(paste("~", group_col)), ncol = 3) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    labs(title = title, fill = str_to_title(gsub("_", " ", fill_col)))
}

# Color palettes ----
get_publication_colors <- function(n, palette = "default") {
  switch(palette,
    "default" = RColorBrewer::brewer.pal(min(n, 11), "Spectral"),
    "immune" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
    "fibroblast" = c("#6ABA88", "#366BA1", "#C576F6", "#D7E6F5"),
    "disease" = c("#56B4E9", "#E69F00", "#E63232"),
    RColorBrewer::brewer.pal(min(n, 11), palette)
  )
}