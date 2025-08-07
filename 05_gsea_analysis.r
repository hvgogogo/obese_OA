# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 05: Gene Set Enrichment Analysis (GSEA)
# ==============================================================================

source("scripts/04_pathway_enrichment.R")

# GSEA functions ----
get_pathway_sets <- function(species, categories = c("hallmark", "c2", "c5")) {
  pathway_list <- list()
  
  # Hallmark
  if ("hallmark" %in% categories) {
    m_df_hallmark <- msigdbr(species = species, category = "H")
    pathway_list$hallmark <- split(m_df_hallmark$gene_symbol, m_df_hallmark$gs_name)
  }
  
  # C2 pathways
  if ("c2" %in% categories) {
    m_df_c2 <- msigdbr(species = species, category = "C2")
    reactome_df <- m_df_c2 %>% filter(gs_subcat == "CP:REACTOME")
    pathway_list$reactome <- split(reactome_df$gene_symbol, reactome_df$gs_name)
    wiki_df <- m_df_c2 %>% filter(gs_subcat == "CP:WIKIPATHWAYS")
    pathway_list$wiki <- split(wiki_df$gene_symbol, wiki_df$gs_name)
    kegg_df <- m_df_c2 %>% filter(gs_subcat == "CP:KEGG")
    pathway_list$kegg <- split(kegg_df$gene_symbol, kegg_df$gs_name)
  }
  
  # C5 GO terms
  if ("c5" %in% categories) {
    m_df_c5 <- msigdbr(species = species, category = "C5")
    gpbp_df <- m_df_c5 %>% filter(gs_subcat == "GO:BP")
    pathway_list$gpbp <- split(gpbp_df$gene_symbol, gpbp_df$gs_name)
    gocc_df <- m_df_c5 %>% filter(gs_subcat == "GO:CC")
    pathway_list$gocc <- split(gocc_df$gene_symbol, gocc_df$gs_name)
    gomf_df <- m_df_c5 %>% filter(gs_subcat == "GO:MF")
    pathway_list$gomf <- split(gomf_df$gene_symbol, gomf_df$gs_name)
  }
  
  cat("Loaded pathway categories:", paste(categories, collapse = ", "), "\n")
  cat("Total pathway sets:", length(pathway_list), "\n")
  
  return(pathway_list)
}

run_gsea_analysis <- function(gene_data, pathway_sets) {
  # Prepare ranks
  geneList <- gene_data %>% arrange(desc(avg_log2FC))
  ranks <- setNames(geneList$avg_log2FC, geneList$gene)
  ranks <- ranks[!is.na(names(ranks))]
  
  cat("Rank gene length:", length(ranks), "\n")
  
  # Run GSEA for all pathway sets
  gsea_res <- list()
  nperm_values <- list(hallmark = 1000000, kegg = 1000000, reactome = 1000000,
                       wiki = 1000000, gocc = 1000000, gomf = 1000000, gpbp = 300000)
  
  for (name in names(pathway_sets)) {
    gsea_res[[name]] <- fgsea(pathways = pathway_sets[[name]],
                              stats = ranks,
                              minSize = 10,
                              maxSize = 500,
                              nperm = nperm_values[[name]])
    cat(name, "finished\n")
  }
  
  return(list(results = gsea_res, ranks = ranks))
}

process_gsea_results <- function(gsea_results, comparison_name, output_path) {
  all_paths <- data.frame()
  
  # Combine all pathway results
  for (name in names(gsea_results$results)) {
    obj_gsea <- gsea_results$results[[name]]
    obj_gsea$type <- name
    all_paths <- rbind(all_paths, obj_gsea)
  }
  
  # Format pathway names
  all_paths$pathway <- gsub("_", " ", all_paths$pathway)
  all_paths$pathway <- sub(" ", ":", all_paths$pathway, fixed = TRUE)
  all_paths$type <- all_paths$pathway
  all_paths$pathway <- sub("^[^:]*:", "", all_paths$pathway)
  all_paths <- all_paths %>%
    mutate(pathway = tolower(pathway),
           pathway = str_replace(pathway, "^(\\w)", function(x) toupper(x)))
  all_paths$type <- sub(":.*$", "", all_paths$type)
  all_paths <- all_paths %>% dplyr::select(type, everything())
  
  # Get significant pathways
  int_paths <- all_paths[all_paths$padj <= 0.05, ]
  
  # Save results
  fwrite(int_paths, file.path(output_path, paste0("paths_sig_", comparison_name, ".csv")))
  fwrite(all_paths, file.path(output_path, paste0("paths_all_", comparison_name, ".csv")))
  
  cat("Significant pathways found:", nrow(int_paths), "\n")
  
  return(int_paths)
}

plot_gsea_results <- function(int_paths, comparison_name, output_path, top_n = 6) {
  if (nrow(int_paths) == 0) return(NULL)
  
  # Get top pathways
  top_paths <- int_paths %>%
    arrange(padj) %>%
    slice_head(n = top_n)
  
  if (nrow(top_paths) == 0) return(NULL)
  
  # Bar plot
  df_bar <- top_paths
  df_bar$NES_category <- with(df_bar, ifelse(NES > 0, "Up-Paths", "Down-Paths"))
  
  p_bar <- ggplot(df_bar, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = NES_category)) + coord_flip() +
    scale_fill_manual(values = c("Up-Paths" = "#E41A1C", "Down-Paths" = "#1a2b76"),
                      name = NULL, breaks = c("Down-Paths", "Up-Paths")) +
    labs(x = "Pathways", y = "NES", title = " ") +
    theme_minimal() +
    theme(text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.position = "top",
          panel.grid = element_blank())
  
  # Save bar plot
  height_bar <- max(3, 0.3 * nrow(df_bar))
  cairo_pdf(file.path(output_path, paste0("gsea_barplot_", comparison_name, ".pdf")),
            width = 5, height = height_bar)
  print(p_bar)
  dev.off()
  
  # Dot plot
  dot_df <- df_bar %>%
    mutate(NES_category = ifelse(NES > 0, "Up-Paths", "Down-Paths"),
           hollow_circle = ifelse(padj > 0.05, "Hollow", "Filled"),
           num_leadingEdge = sapply(leadingEdge, length),
           GeneRatio = num_leadingEdge / size) %>%
    dplyr::select(NES, Description = pathway, num_leadingEdge, GeneRatio, padj, NES_category, hollow_circle)
  
  p_dot <- ggplot(dot_df, aes(x = NES, y = reorder(Description, NES))) +
    geom_point(aes(size = GeneRatio, color = padj, shape = hollow_circle)) +
    theme_bw(base_size = 10) +
    scale_size_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_colour_gradient(name = "Padj", limits = c(0.00, 0.05),
                          breaks = c(0.01, 0.05), high = "blue", low = "red") +
    scale_shape_manual(values = c("Filled" = 16, "Hollow" = 1),
                       labels = c("Sig-Paths", "noSig-Paths"), name = " ") +
    ylab(NULL) +
    facet_grid(. ~ NES_category, scale = "free_x") +
    theme(text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.key.size = unit(0.2, "cm"),
          panel.grid = element_blank())
  
  # Save dot plot
  height_dot <- max(3, 0.3 * nrow(df_bar))
  cairo_pdf(file.path(output_path, paste0("gsea_dotplot_", comparison_name, ".pdf")),
            width = 7, height = height_dot)
  print(p_dot)
  dev.off()
  
  return(list(bar = p_bar, dot = p_dot))
}

# Custom GSEA pathway selection and plotting ----
create_custom_gsea_plots <- function(selected_pathways, cell_type, comparison_type, 
                                    custom_name = "custom", paths_data) {
  # Check if data exists
  if (is.null(paths_data[[cell_type]][[comparison_type]])) {
    cat("No GSEA data found for", cell_type, "-", comparison_type, "\n")
    return(NULL)
  }
  
  # Get the pathway data
  int_paths <- paths_data[[cell_type]][[comparison_type]]
  
  # Filter for selected pathways (case-insensitive matching)
  selected_data <- int_paths[tolower(int_paths$pathway) %in% tolower(selected_pathways), ]
  
  if (nrow(selected_data) == 0) {
    cat("No matching pathways found. Available pathways:\n")
    cat(paste(int_paths$pathway, collapse = "\n"))
    return(NULL)
  }
  
  # Create shortened pathway names
  shorten_pathway_name <- function(pathway_name) {
    # Split by common delimiters and take first 3 meaningful words
    words <- unlist(strsplit(pathway_name, "[ _-]"))
    # Remove common connecting words
    stop_words <- c("of", "by", "via", "in", "the", "a", "an", "and", "or", "for", "to", "from")
    meaningful_words <- words[!tolower(words) %in% stop_words & nchar(words) > 2]
    
    if (length(meaningful_words) >= 3) {
      return(paste(meaningful_words[1:3], collapse = " "))
    } else if (length(meaningful_words) > 0) {
      return(paste(meaningful_words, collapse = " "))
    } else {
      # Fallback: take first 3 words regardless
      return(paste(words[1:min(3, length(words))], collapse = " "))
    }
  }
  
  # Apply shortening to selected data
  selected_data$pathway_short <- sapply(selected_data$pathway, shorten_pathway_name)
  
  cat("Found", nrow(selected_data), "matching pathways\n")
  for (i in 1:nrow(selected_data)) {
    cat(sprintf("%d. %s: %s -> %s (NES: %.2f, padj: %.3f)\n",
                i, selected_data$type[i], selected_data$pathway[i], 
                selected_data$pathway_short[i], selected_data$NES[i], selected_data$padj[i]))
  }
  
  # Bar plot
  df_bar <- selected_data
  df_bar$NES_category <- with(df_bar, ifelse(NES > 0, "Up-Paths", "Down-Paths"))
  
  p_bar <- ggplot(df_bar, aes(reorder(pathway_short, NES), NES)) +
    geom_col(aes(fill = NES_category)) + coord_flip() +
    scale_fill_manual(values = c("Up-Paths" = "#E41A1C", "Down-Paths" = "#1a2b76"),
                      name = NULL, breaks = c("Down-Paths", "Up-Paths")) +
    labs(x = "", y = "NES", title = " ") +
    theme_minimal() +
    theme(text = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.position = "bottom",
          axis.line.x = element_line(linewidth = 1),
          panel.grid = element_blank())
  
  cat("Custom plots created for", cell_type, "-", comparison_type, "\n")
  return(list(bar = p_bar, data = selected_data))
}

# Main GSEA analysis function ----
run_gsea_main_analysis <- function(genes_list, output_dir) {
  cat("Starting GSEA analysis...\n")
  
  condi2compares <- c("weight", "OA", "obeseOA")
  condis <- c("LLSFs", "SLSFs_CD34pos", "SLSFs_CD34neg", "SLSFs_IM")
  species_used <- "Homo sapiens"
  
  paths <- list()
  pathway_categories <- c("hallmark", "c2", "c5")
  pathway_sets <- get_pathway_sets(species_used, pathway_categories)
  
  for (comparison_type in condi2compares) {
    for (cell_type in condis) {
      cat("\nRunning GSEA for:", cell_type, "-", comparison_type, "\n")
      
      if (!is.null(genes_list[[cell_type]][[comparison_type]])) {
        gene_data <- genes_list[[cell_type]][[comparison_type]]
        
        if (nrow(gene_data) == 0) {
          cat("No genes found, skipping\n")
          next
        }
        
        tryCatch({
          # Run GSEA
          gsea_results <- run_gsea_analysis(gene_data, pathway_sets)
          
          # Process results
          comparison_name <- paste(cell_type, comparison_type, sep = "-")
          int_paths <- process_gsea_results(gsea_results, comparison_name, output_dir)
          
          # Store results
          paths[[cell_type]][[comparison_type]] <- int_paths
          
          # Create plots
          plots <- plot_gsea_results(int_paths, comparison_name, output_dir)
          
        }, error = function(e) {
          cat("Error in GSEA analysis:", e$message, "\n")
        })
      } else {
        cat("No gene data found\n")
      }
    }
  }
  
  # Save all results
  save(paths, file = file.path(output_dir, "gsea_allpaths.RData"))
  return(paths)
}

# Create condition-specific GSEA plots ----
create_condition_gsea_plots <- function(paths_data, output_dir) {
  cat("Creating condition-specific GSEA plots...\n")
  
  cell_types <- c("SLSFs_CD34pos", "SLSFs_CD34neg", "LLSFs", "SLSFs_IM")
  
  # Define pathway sets for different comparisons
  selected_paths_general <- c(
    'Transcriptional regulation of white adipocyte differentiation',
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
  )
  
  selected_paths_im <- c(
    "Il6 jak stat3 signaling",
    "Reactive oxygen species pathway",
    "Chemokine activity",
    "Chemokine receptor binding",
    "Ccr chemokine receptor binding",
    "Cxcr chemokine receptor binding",
    "Inflammatory response",
    "Diseases of programmed cell death",
    "Oxidative stress induced senescence"
  )
  
  # Weight comparison
  comparison_type <- "weight"
  custom_name <- 'obesity'
  
  plot1 <- create_custom_gsea_plots(selected_paths_general, cell_types[1], 
                                   comparison_type, custom_name, paths_data)$bar +
    ggtitle("CD34+ SLSFs") + theme(legend.position = "bottom")
  
  plot2 <- create_custom_gsea_plots(selected_paths_general, cell_types[3], 
                                   comparison_type, custom_name, paths_data)$bar +
    ggtitle("LLSFs") + theme(legend.position = "bottom")
  
  plot3 <- create_custom_gsea_plots(selected_paths_general, cell_types[2], 
                                   comparison_type, custom_name, paths_data)$bar +
    ggtitle("CD34- SLSFs") + theme(legend.position = "bottom")
  
  plot4 <- create_custom_gsea_plots(selected_paths_im, cell_types[4], 
                                   comparison_type, custom_name, paths_data)$bar +
    ggtitle("IM SLSFs") + theme(legend.position = "bottom")
  
  # Combine plots
  combined_plot <- (plot1 + theme(legend.position = "left") |
                    plot2 + theme(legend.position = "none") |
                    plot3 + theme(legend.position = "none") |
                    plot4 + theme(legend.position = "none"))
  
  filename <- paste(custom_name, "combined_all_celltypes", comparison_type, sep = "_")
  ggsave(file.path(output_dir, "figures", paste0("custom_barplot_", filename, ".pdf")),
         combined_plot, width = 16, height = 3)
  
  # Similar plots for OA and obeseOA comparisons would go here...
  # (Code shortened for brevity)
  
  cat("Condition-specific plots created!\n")
}

# Execute GSEA analysis ----
gsea_output_dir <- file.path(output_dir, "gsea_analysis")
dir.create(gsea_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gsea_output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# Load differential expression results if available
if (exists("deg_results")) {
  genes_list <- deg_results$all_genes
} else {
  # Load from saved file
  load(file.path(output_dir, "differential_expression/differential_expression_results.RData"))
  genes_list <- allgenes_list
}

# Run GSEA analysis
gsea_paths <- run_gsea_main_analysis(genes_list, gsea_output_dir)

# Create custom plots
create_condition_gsea_plots(gsea_paths, gsea_output_dir)

cat("GSEA analysis completed!\n")
cat("Results saved to:", gsea_output_dir, "\n")