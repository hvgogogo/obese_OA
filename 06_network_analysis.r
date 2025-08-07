# ==============================================================================
# Single-cell RNA-seq Analysis of Synovial Fibroblasts in Obesity and OA
# Script 06: Network Analysis and Hub Gene Identification
# ==============================================================================

source("scripts/05_gsea_analysis.R")

# Network analysis functions ----
get_pathway_genes <- function(pathway_name, pathway_type = "hallmark") {
  cat("Extracting genes from", pathway_type, "pathway:", pathway_name, "\n")
  
  # Load pathway sets
  if (pathway_type == "hallmark") {
    pathway_sets <- get_pathway_sets("Homo sapiens", c("hallmark"))
    pathway_db <- pathway_sets$hallmark
  } else if (pathway_type %in% c("reactome", "kegg", "wiki")) {
    pathway_sets <- get_pathway_sets("Homo sapiens", c("c2"))
    pathway_db <- switch(pathway_type,
                        "reactome" = pathway_sets$reactome,
                        "kegg" = pathway_sets$kegg,
                        "wiki" = pathway_sets$wiki)
  } else {
    stop("Unsupported pathway type")
  }
  
  # Try to find pathway with HALLMARK_ prefix for hallmark pathways
  if (pathway_type == "hallmark" && !grepl("^HALLMARK_", pathway_name)) {
    pathway_name <- paste0("HALLMARK_", toupper(gsub(" ", "_", pathway_name)))
  }
  
  # Extract genes
  if (pathway_name %in% names(pathway_db)) {
    genes <- pathway_db[[pathway_name]]
    cat("Found", length(genes), "genes in pathway:", pathway_name, "\n")
    return(list(genes = genes, pathway_name = pathway_name, gene_count = length(genes)))
  } else {
    cat("Pathway not found. Available pathways:\n")
    cat(paste(head(names(pathway_db), 10), collapse = "\n"), "\n")
    return(NULL)
  }
}

build_ppi_network <- function(gene_res, score_threshold = 400) {
  cat("Building PPI network for", length(gene_res$genes), "genes\n")
  
  if (is.null(gene_res) || length(gene_res$genes) == 0) {
    stop("No genes provided")
  }
  
  # Initialize STRING database
  string_db <- STRINGdb$new(version = "12.0", species = 9606,
                            score_threshold = score_threshold,
                            input_directory = "")
  
  # Map genes to STRING IDs
  gene_df <- data.frame(gene = unique(gene_res$genes), stringsAsFactors = FALSE)
  mapped_genes <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
  
  if (nrow(mapped_genes) == 0) {
    stop("No genes could be mapped to STRING database")
  }
  
  cat("Mapped", nrow(mapped_genes), "out of", length(gene_res$genes), "genes\n")
  
  # Get protein interactions
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  if (nrow(interactions) == 0) {
    stop("No interactions found")
  }
  
  cat("Found", nrow(interactions), "interactions\n")
  
  # Create network graph
  g <- graph_from_data_frame(interactions, directed = FALSE)
  
  # Add gene symbols as vertex attributes
  string_to_gene <- setNames(mapped_genes$gene, mapped_genes$STRING_id)
  V(g)$gene_symbol <- string_to_gene[V(g)$name]
  V(g)$label <- V(g)$gene_symbol
  E(g)$weight <- interactions$combined_score
  
  cat("Network created: ", vcount(g), "nodes,", ecount(g), "edges\n")
  
  return(list(graph = g, mapped_genes = mapped_genes, string_db = string_db))
}

calculate_centrality_measures <- function(network_obj) {
  cat("Calculating centrality measures...\n")
  
  g <- network_obj$graph
  mapped_genes <- network_obj$mapped_genes
  
  # Match STRING IDs in graph to mapped genes
  vertex_genes <- data.frame(
    STRING_id = V(g)$name,
    gene = V(g)$gene_symbol,
    stringsAsFactors = FALSE
  )
  
  # Calculate centrality measures
  centrality_results <- vertex_genes %>%
    mutate(
      degree = degree(g),
      betweenness = betweenness(g, normalized = TRUE),
      closeness = closeness(g, normalized = TRUE),
      eigenvector = eigen_centrality(g)$vector,
      pagerank = page_rank(g)$vector
    ) %>%
    mutate(
      degree_rank = rank(-degree, ties.method = "min"),
      betweenness_rank = rank(-betweenness, ties.method = "min"),
      closeness_rank = rank(-closeness, ties.method = "min"),
      eigenvector_rank = rank(-eigenvector, ties.method = "min"),
      pagerank_rank = rank(-pagerank, ties.method = "min"),
      combined_rank = (degree_rank + betweenness_rank + closeness_rank +
                        eigenvector_rank + pagerank_rank) / 5
    ) %>%
    arrange(combined_rank)
  
  cat("Centrality measures calculated for", nrow(centrality_results), "genes\n")
  return(centrality_results)
}

identify_hub_genes <- function(centrality_results, n_top = 10) {
  cat("Identifying top", n_top, "hub genes...\n")
  
  hub_analysis <- list(
    top_degree = head(centrality_results[order(-centrality_results$degree), ], n_top),
    top_betweenness = head(centrality_results[order(-centrality_results$betweenness), ], n_top),
    top_closeness = head(centrality_results[order(-centrality_results$closeness), ], n_top),
    top_eigenvector = head(centrality_results[order(-centrality_results$eigenvector), ], n_top),
    top_combined = head(centrality_results, n_top)
  )
  
  # Print top combined genes
  cat("Top hub genes (combined ranking):\n")
  print(hub_analysis$top_combined[, c("gene", "degree", "betweenness", "combined_rank")])
  
  return(hub_analysis)
}

plot_hub_interactions_force <- function(centrality_results, network_obj, n_top = 3) {
  library(ggraph)
  library(ggplot2)
  
  g <- network_obj$graph
  top_genes <- head(centrality_results, n_top)
  hub_string_ids <- top_genes$STRING_id
  
  # Get edges connected to hubs only
  all_edges <- as_edgelist(g)
  hub_edges_idx <- which(all_edges[,1] %in% hub_string_ids | all_edges[,2] %in% hub_string_ids)
  hub_edges <- all_edges[hub_edges_idx, ]
  
  # Create hub-only interaction graph
  g_hub_only <- graph_from_edgelist(hub_edges, directed = FALSE)
  
  # Add attributes
  string_to_gene <- setNames(V(g)$gene_symbol, V(g)$name)
  V(g_hub_only)$gene_symbol <- string_to_gene[V(g_hub_only)$name]
  V(g_hub_only)$is_hub <- V(g_hub_only)$name %in% hub_string_ids
  V(g_hub_only)$original_degree <- degree(g)[V(g_hub_only)$name] # From original network
  
  ggraph(g_hub_only, layout = "fr") + # Force-directed layout
    geom_edge_link(alpha = 0.6, color = "steelblue", width = 0.8) +
    geom_node_point(aes(size = original_degree,
                       color = is_hub,
                       stroke = ifelse(is_hub, 2, 0.5)),
                   alpha = 0.8) +
    geom_node_text(aes(label = ifelse(!is_hub, gene_symbol, "")),
                   repel = TRUE, size = 3, fontface = "bold") +
    geom_node_label(aes(label = ifelse(is_hub, gene_symbol, "")),
                   fill = "white", alpha = 0.8,
                   label.padding = unit(0.2, "lines"),
                   label.r = unit(0, "lines"), # Square corners
                   repel = TRUE, size = 3, fontface = "bold") +
    scale_color_manual(values = c("FALSE" = "#4ECDC4", "TRUE" = "#FF6B6B"),
                      name = "Hub gene") + # Changed legend title
    scale_size_continuous(range = c(2, 8), name = "Degree") + # Reduced node sizes and added legend title
    theme_graph(base_family = "sans") +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10, family = "sans")) +
    labs(subtitle = paste("Hubs:", paste(top_genes$gene[1:min(5, n_top)], collapse = ", ")))
}

generate_summary_report <- function(gene_res, network_obj, centrality_results, hub_analysis) {
  cat("\n=== PATHWAY PPI ANALYSIS SUMMARY ===\n")
  cat("Pathway:", gene_res$pathway_name, "\n")
  cat("Total genes in pathway:", gene_res$gene_count, "\n")
  cat("Genes mapped to PPI network:", vcount(network_obj$graph), "\n")
  cat("Total interactions:", ecount(network_obj$graph), "\n")
  cat("Network density:", edge_density(network_obj$graph), "\n")
  cat("Connected components:", components(network_obj$graph)$no, "\n")
  
  cat("\nTop 5 Hub Genes:\n")
  top5 <- head(centrality_results, 5)
  for(i in 1:nrow(top5)) {
    cat(sprintf("%d. %s (degree: %d, combined_rank: %.2f)\n",
                i, top5$gene[i], top5$degree[i], top5$combined_rank[i]))
  }
  
  cat("\n=== END SUMMARY ===\n")
  
  # Return summary data
  summary_data <- list(
    pathway_info = list(
      pathway_name = gene_res$pathway_name,
      total_genes = gene_res$gene_count,
      mapped_genes = vcount(network_obj$graph)
    ),
    network_stats = list(
      nodes = vcount(network_obj$graph),
      edges = ecount(network_obj$graph),
      density = edge_density(network_obj$graph),
      components = components(network_obj$graph)$no
    ),
    top_hubs = top5
  )
  
  return(summary_data)
}

# Main network analysis function ----
run_pathway_network_analysis <- function(pathway_name, pathway_type = "hallmark", 
                                        score_threshold = 400, output_dir) {
  cat(sprintf("=== Analyzing %s pathway: %s ===\n", pathway_type, pathway_name))
  
  # Extract pathway genes
  gene_res <- get_pathway_genes(pathway_name, pathway_type)
  if (is.null(gene_res)) {
    cat("Failed to extract genes for pathway:", pathway_name, "\n")
    return(NULL)
  }
  
  # Build PPI network
  tryCatch({
    network_obj <- build_ppi_network(gene_res, score_threshold)
  }, error = function(e) {
    cat("Error building PPI network:", e$message, "\n")
    return(NULL)
  })
  
  # Calculate centrality measures
  centrality_results <- calculate_centrality_measures(network_obj)
  
  # Identify hub genes
  hub_analysis <- identify_hub_genes(centrality_results, 15)
  
  # Generate summary report
  summary_report <- generate_summary_report(gene_res, network_obj, centrality_results, hub_analysis)
  
  # Create visualization
  plot_hub <- plot_hub_interactions_force(centrality_results, network_obj, n_top = 3)
  
  # Save plot
  filename <- paste(pathway_name, pathway_type, sep = "_")
  filename <- gsub(" ", "_", filename)
  ggsave(file.path(output_dir, "figures", paste0("hub_network_", filename, ".pdf")),
         plot_hub, width = 8, height = 8, device = cairo_pdf)
  
  # Save results
  results <- list(
    gene_res = gene_res,
    network_obj = network_obj,
    centrality_results = centrality_results,
    hub_analysis = hub_analysis,
    summary_report = summary_report,
    plot = plot_hub
  )
  
  save(results, file = file.path(output_dir, paste0("network_analysis_", filename, ".RData")))
  
  cat(sprintf("Network analysis completed for %s\n", pathway_name))
  cat("Results saved to:", file.path(output_dir, paste0("network_analysis_", filename, ".RData")), "\n")
  
  return(results)
}

# Execute network analysis ----
network_output_dir <- file.path(output_dir, "network_analysis")
dir.create(network_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(network_output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

cat("Starting network analysis...\n")

# Analyze key pathways
key_pathways <- c("Adipogenesis", "Hypoxia", "Inflammatory response", "TNF-alpha signaling via NF-kappaB")

network_results <- list()
for (pathway in key_pathways) {
  cat(sprintf("\n=== Analyzing pathway: %s ===\n", pathway))
  
  result <- run_pathway_network_analysis(
    pathway_name = pathway,
    pathway_type = "hallmark",
    score_threshold = string_score_threshold,
    output_dir = network_output_dir
  )
  
  if (!is.null(result)) {
    network_results[[pathway]] <- result
  }
}

# Save combined results
save(network_results, file = file.path(network_output_dir, "all_network_results.RData"))

# Create summary table of all hub genes
cat("\n=== Creating hub gene summary ===\n")
hub_summary <- data.frame()

for (pathway in names(network_results)) {
  if (!is.null(network_results[[pathway]])) {
    top_hubs <- network_results[[pathway]]$summary_report$top_hubs
    top_hubs$pathway <- pathway
    hub_summary <- rbind(hub_summary, top_hubs[1:3, ]) # Top 3 from each pathway
  }
}

if (nrow(hub_summary) > 0) {
  # Remove duplicates and rank by combined score
  hub_summary_unique <- hub_summary %>%
    group_by(gene) %>%
    slice_min(combined_rank, n = 1) %>%
    ungroup() %>%
    arrange(combined_rank)
  
  write.csv(hub_summary_unique, 
            file.path(network_output_dir, "hub_genes_summary.csv"), 
            row.names = FALSE)
  
  cat("Hub gene summary saved to:", file.path(network_output_dir, "hub_genes_summary.csv"), "\n")
  cat("\nTop 10 hub genes across all pathways:\n")
  print(head(hub_summary_unique[, c("gene", "degree", "betweenness", "combined_rank", "pathway")], 10))
}

cat("\nNetwork analysis pipeline completed!\n")
cat("Results saved to:", network_output_dir, "\n")