# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Utility Functions: Pathway Analysis
# ==============================================================================

# Pathway abbreviation mapping ----
pathway_abbr <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammation",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "IFN-γ Response", 
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFα via NF-κB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6-JAK-STAT3",
  "HALLMARK_ADIPOGENESIS" = "Adipogenesis",
  "HALLMARK_FATTY_ACID_METABOLISM" = "Fatty Acid Metabolism",
  "HALLMARK_GLYCOLYSIS" = "Glycolysis",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Oxidative Phosphorylation",
  "HALLMARK_HYPOXIA" = "Hypoxia",
  "HALLMARK_APOPTOSIS" = "Apoptosis",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_TGF_BETA_SIGNALING" = "TGF-β",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING" = "WNT/β-catenin",
  "HALLMARK_ANGIOGENESIS" = "Angiogenesis",
  "HALLMARK_MYOGENESIS" = "Myogenesis",
  "HALLMARK_G2M_CHECKPOINT" = "G2/M Checkpoint",
  "HALLMARK_MITOTIC_SPINDLE" = "Mitotic Spindle",
  "HALLMARK_MYC_TARGETS_V1" = "MYC Targets"
)

# Selected pathway sets for different analyses ----
get_selected_pathways <- function(analysis_type = "general") {
  switch(analysis_type,
    "obesity" = c(
      "Transcriptional regulation of white adipocyte differentiation",
      "Mitochondrial translation",
      "Signaling by vegf",
      "Adipogenesis",
      "Gpcr ligand binding",
      "Oxidative stress induced senescence",
      "Inflammatory response",
      "Tnfa signaling via nfkb",
      "Runx1 regulates transcription of genes involved in differentiation of hscs",
      "Glycolysis",
      "Hypoxia"
    ),
    "inflammation" = c(
      "Il6 jak stat3 signaling",
      "Reactive oxygen species pathway",
      "Chemokine activity",
      "Chemokine receptor binding",
      "Ccr chemokine receptor binding",
      "Cxcr chemokine receptor binding",
      "Inflammatory response",
      "Diseases of programmed cell death",
      "Oxidative stress induced senescence"
    ),
    "osteoarthritis" = c(
      "Transcriptional regulation of white adipocyte differentiation",
      "Rrna processing in the mitochondrion",
      "Excitatory synapse",
      "Synaptic membrane",
      "Potassium channels",
      "Fatty acid metabolism",
      "Collagen biosynthesis and modifying enzymes",
      "Collagen chain trimerization",
      "Oxidative phosphorylation",
      "Signaling by retinoic acid",
      "Lipid droplet",
      "Osteopontin signaling"
    ),
    "general" = c(
      "Inflammatory response",
      "Interferon gamma response",
      "Fatty acid metabolism",
      "Glycolysis",
      "Oxidative phosphorylation",
      "G2M checkpoint",
      "Mitotic spindle",
      "MYC targets v1",
      "Angiogenesis",
      "Adipogenesis",
      "Epithelial mesenchymal transition",
      "Myogenesis",
      "Hypoxia",
      "Tnfa signaling via nfkb",
      "Apoptosis",
      "WNT beta catenin signaling",
      "Tgf beta signaling",
      "Il6 jak stat3 signaling"
    )
  )
}

# Pathway color coding ----
get_pathway_colors <- function() {
  list(
    red_pathways = "Hypoxia",
    orange_pathways = c("Inflammation", "TNFα via NF-κB", "IL6-JAK-STAT3"),
    blue_pathways = c("Adipogenesis", "Myogenesis", "EMT", "TGF-β"),
    black_pathways = "Apoptosis"
  )
}

# Enhanced pathway enrichment dot plot ----
create_pathway_dotplot <- function(enr_paths, condition_name, pathway_type = "HALLMARK",
                                  clusters = NULL, cluster_labels = NULL) {
  
  if (is.null(clusters)) {
    clusters <- c("SLSFs_CD34pos", "SLSFs_CD34neg", "LLSFs", "SLSFs_IM")
  }
  if (is.null(cluster_labels)) {
    cluster_labels <- c("CD34pos SLSFs", "CD34neg SLSFs", "LLSFs", "IM SLSFs")
  }
  
  # Prepare data
  data_list <- list()
  for (cluster in clusters) {
    if (!is.null(enr_paths[[cluster]])) {
      data_list[[cluster]] <- enr_paths[[cluster]]
      data_list[[cluster]]$cluster <- cluster
    }
  }
  
  # Combine and filter data
  df_bar <- do.call(rbind, data_list) %>%
    dplyr::filter(type == pathway_type) %>%
    dplyr::mutate(
      condition = factor(cluster, levels = clusters),
      p.adjust = -log10(p.adjust)
    )
  
  # Apply pathway abbreviations if available
  if (exists("pathway_abbr") && length(pathway_abbr) > 1) {
    df_bar <- df_bar %>%
      dplyr::mutate(pathway = dplyr::recode(pathway, !!!pathway_abbr))
  }
  
  # Filter significant pathways
  df_bar <- df_bar %>%
    dplyr::filter(p.adjust > 2) %>%
    dplyr::mutate(p.adjust = pmin(ifelse(is.infinite(p.adjust), 300, p.adjust), 20))
  
  # Get pathway colors
  pathway_colors <- get_pathway_colors()
  get_pathway_color <- function(pathway_name) {
    if (pathway_name %in% pathway_colors$red_pathways) return("red")
    if (pathway_name %in% pathway_colors$orange_pathways) return("orange")
    if (pathway_name %in% pathway_colors$blue_pathways) return("blue")
    if (pathway_name %in% pathway_colors$black_pathways) return("black")
    return("black")
  }
  
  # Clean data and set pathway levels
  df_bar <- df_bar %>%
    dplyr::filter(!is.na(pathway)) %>%
    dplyr::filter(pathway != "") %>%
    droplevels() %>%
    dplyr::select(pathway, p.adjust, condition, GeneRatio)
  
  if (nrow(df_bar) == 0) {
    cat("No significant pathways found for", condition_name, "\n")
    return(NULL)
  }
  
  # Get colors for pathways actually present
  present_pathways <- levels(forcats::fct_rev(factor(df_bar$pathway)))
  present_colors <- sapply(present_pathways, get_pathway_color)
  
  # Create dot plot
  ggplot(df_bar, aes(condition, forcats::fct_rev(factor(pathway)), 
                    fill = p.adjust, size = GeneRatio)) +
    geom_point(shape = 21, stroke = 0.1) +
    geom_vline(xintercept = seq(1.5, length(clusters), 1), linewidth = 0.5) +
    geom_hline(yintercept = seq(1.5, length(present_pathways), 1), linewidth = 0.5) +
    scale_radius(range = c(2, 6)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    scale_x_discrete(labels = cluster_labels) +
    labs(
      title = condition_name,
      size = "Gene ratio",
      fill = "-log10(p_adj)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold", color = present_colors),
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10, face = "bold")
    )
}

# GSEA bar plot creation ----
create_gsea_barplot <- function(gsea_data, selected_pathways, cell_type, comparison_type,
                               title = NULL) {
  
  # Filter for selected pathways
  selected_data <- gsea_data[tolower(gsea_data$pathway) %in% tolower(selected_pathways), ]
  
  if (nrow(selected_data) == 0) {
    cat("No matching pathways found for", cell_type, "-", comparison_type, "\n")
    return(NULL)
  }
  
  # Shorten pathway names
  shorten_pathway_name <- function(pathway_name) {
    words <- unlist(strsplit(pathway_name, "[ _-]"))
    stop_words <- c("of", "by", "via", "in", "the", "a", "an", "and", "or", "for", "to", "from")
    meaningful_words <- words[!tolower(words) %in% stop_words & nchar(words) > 2]
    
    if (length(meaningful_words) >= 3) {
      return(paste(meaningful_words[1:3], collapse = " "))
    } else if (length(meaningful_words) > 0) {
      return(paste(meaningful_words, collapse = " "))
    } else {
      return(paste(words[1:min(3, length(words))], collapse = " "))
    }
  }
  
  selected_data$pathway_short <- sapply(selected_data$pathway, shorten_pathway_name)
  selected_data$NES_category <- ifelse(selected_data$NES > 0, "Up-Paths", "Down-Paths")
  
  if (is.null(title)) {
    title <- paste(cell_type, comparison_type, sep = " - ")
  }
  
  ggplot(selected_data, aes(reorder(pathway_short, NES), NES)) +
    geom_col(aes(fill = NES_category)) + 
    coord_flip() +
    scale_fill_manual(values = c("Up-Paths" = "#E41A1C", "Down-Paths" = "#1a2b76"),
                     name = NULL, breaks = c("Down-Paths", "Up-Paths")) +
    labs(x = "", y = "NES", title = title) +
    theme_minimal() +
    theme(
      text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      axis.line.x = element_line(linewidth = 1),
      panel.grid = element_blank()
    )
}

# Pathway overlap analysis ----
calculate_pathway_overlap <- function(pathway_list1, pathway_list2, database = "hallmark") {
  
  # Get pathway gene sets
  pathway_sets <- get_pathway_sets("Homo sapiens", c(database))
  db_pathways <- pathway_sets[[database]]
  
  # Calculate Jaccard similarity for each pathway pair
  overlap_results <- data.frame(
    pathway1 = character(),
    pathway2 = character(),
    overlap_genes = integer(),
    jaccard_similarity = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (p1 in pathway_list1) {
    for (p2 in pathway_list2) {
      if (p1 %in% names(db_pathways) && p2 %in% names(db_pathways)) {
        genes1 <- db_pathways[[p1]]
        genes2 <- db_pathways[[p2]]
        
        overlap <- length(intersect(genes1, genes2))
        union <- length(union(genes1, genes2))
        jaccard <- overlap / union
        
        overlap_results <- rbind(overlap_results, data.frame(
          pathway1 = p1,
          pathway2 = p2,
          overlap_genes = overlap,
          jaccard_similarity = jaccard
        ))
      }
    }
  }
  
  return(overlap_results)
}

# Pathway gene set utilities ----
get_pathway_genes_list <- function(pathway_names, database = "hallmark") {
  pathway_sets <- get_pathway_sets("Homo sapiens", c(database))
  db_pathways <- pathway_sets[[database]]
  
  result <- list()
  for (pathway in pathway_names) {
    # Handle different naming conventions
    pathway_key <- pathway
    if (database == "hallmark" && !grepl("^HALLMARK_", pathway)) {
      pathway_key <- paste0("HALLMARK_", toupper(gsub(" ", "_", pathway)))
    }
    
    if (pathway_key %in% names(db_pathways)) {
      result[[pathway]] <- db_pathways[[pathway_key]]
    } else {
      cat("Warning: Pathway", pathway, "not found in", database, "database\n")
    }
  }
  
  return(result)
}

# Pathway enrichment summary ----
summarize_pathway_enrichment <- function(enrichment_results, top_n = 10) {
  
  summary_list <- list()
  
  for (condition in names(enrichment_results)) {
    for (cluster in names(enrichment_results[[condition]])) {
      
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        pathways <- enrichment_results[[condition]][[cluster]]
        
        if (nrow(pathways) > 0) {
          top_pathways <- pathways %>%
            arrange(p.adjust) %>%
            slice_head(n = top_n) %>%
            select(pathway, p.adjust, Count, GeneRatio)
          
          summary_list[[condition]][[cluster]] <- list(
            total_pathways = nrow(pathways),
            significant_pathways = sum(pathways$p.adjust <= 0.05),
            top_pathways = top_pathways
          )
        }
      }
    }
  }
  
  return(summary_list)
}

# Create pathway comparison heatmap ----
create_pathway_comparison_heatmap <- function(enrichment_results, pathway_list = NULL,
                                            conditions = NULL, clusters = NULL) {
  
  if (is.null(conditions)) {
    conditions <- names(enrichment_results)
  }
  if (is.null(clusters)) {
    clusters <- names(enrichment_results[[1]])
  }
  
  # Create matrix for heatmap
  comparison_matrix <- matrix(NA, nrow = length(pathway_list), 
                             ncol = length(conditions) * length(clusters))
  
  rownames(comparison_matrix) <- pathway_list
  col_names <- paste(rep(conditions, each = length(clusters)), 
                    rep(clusters, length(conditions)), sep = "_")
  colnames(comparison_matrix) <- col_names
  
  # Fill matrix with -log10(p.adjust) values
  for (i in seq_along(conditions)) {
    condition <- conditions[i]
    for (j in seq_along(clusters)) {
      cluster <- clusters[j]
      col_idx <- (i-1) * length(clusters) + j
      
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        pathways_data <- enrichment_results[[condition]][[cluster]]
        
        for (pathway in pathway_list) {
          pathway_match <- pathways_data[tolower(pathways_data$pathway) == tolower(pathway), ]
          if (nrow(pathway_match) > 0) {
            comparison_matrix[pathway, col_idx] <- -log10(pathway_match$p.adjust[1])
          }
        }
      }
    }
  }
  
  # Create heatmap
  pheatmap::pheatmap(comparison_matrix,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "none",
                    color = colorRampPalette(c("white", "yellow", "red"))(100),
                    na_col = "grey90",
                    fontsize = 10,
                    angle_col = 45,
                    main = "Pathway Enrichment Comparison\n(-log10 p.adjust)")
}

# Export pathway results ----
export_pathway_results <- function(enrichment_results, output_dir, format = "csv") {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (condition in names(enrichment_results)) {
    for (cluster in names(enrichment_results[[condition]])) {
      
      if (!is.null(enrichment_results[[condition]][[cluster]])) {
        pathways_data <- enrichment_results[[condition]][[cluster]]
        filename <- paste(condition, cluster, "pathways", sep = "_")
        
        if (format == "csv") {
          write.csv(pathways_data, 
                   file.path(output_dir, paste0(filename, ".csv")), 
                   row.names = FALSE)
        } else if (format == "xlsx") {
          openxlsx::write.xlsx(pathways_data, 
                              file.path(output_dir, paste0(filename, ".xlsx")))
        }
      }
    }
  }
  
  cat("Pathway results exported to:", output_dir, "\n")
}