# DESeq2 Custom Gene Heatmap Generator
# Usage: 
# 1. Change the file paths under "CONFIGURATION" section
# 2. Set USE_CUSTOM_GENES to TRUE and add your genes to CUSTOM_GENE_LIST
# 3. Adjust the display options to control sample labels
# 4. Run the script

# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(readr)
library(viridis)

###### CONFIGURATION - MODIFY THESE VALUES ######
# Input/Output files
DESEQ_RESULTS_FILE <- "uterus_deseq2_results.csv"  # Your DESeq2 results file
OUTPUT_FILE <- "gene_heatmap.png"          # Output heatmap image file

# Gene selection options
USE_CUSTOM_GENES <- TRUE                  # Set to TRUE to use custom gene list
CUSTOM_GENE_LIST <- c("TMEM184A", "PDE1C", "GPR146", "FAM110B")  # Your genes of interest
CUSTOM_GENE_FILE <- NULL                   # Optional: path to file with genes (one per line)
TOP_GENES_COUNT <- 10                      # Number of top genes to show (if not using custom genes)

# Heatmap appearance
SCALE_METHOD <- "row"                      # "row", "column", or "none"
CLUSTER_ROWS <- TRUE                       # Cluster genes
CLUSTER_COLS <- TRUE                       # Cluster samples
SHOW_GENE_NAMES <- TRUE                    # Show gene names on heatmap
SHOW_SAMPLE_NAMES <- FALSE                 # Set to FALSE to hide sample names completely
SAMPLE_NAME_ANGLE <- 0                     # Use 0 for horizontal, 45 for angled, 90 for vertical text
COLOR_SCHEME <- "RdBu"                     # Color scheme (e.g., "RdBu", "viridis", "YlOrRd")
################################################

# Function to load genes from file
load_genes_from_file <- function(file_path) {
  if (!is.null(file_path) && file.exists(file_path)) {
    genes <- readLines(file_path)
    genes <- genes[genes != ""]  # Remove empty lines
    return(genes)
  }
  return(NULL)
}

# Main function
create_heatmap <- function() {
  cat("Loading DESeq2 results...\n")
  
  # Read DESeq2 results
  results <- read_csv(DESEQ_RESULTS_FILE)
  names(results)[1] <- "gene_name"
  
  
  # Ensure results has gene names column
  if (!"gene" %in% colnames(results)) {
    # Try to find gene name column
    potential_gene_cols <- c("gene_name", "Gene", "X", "rowname")
    for (col in potential_gene_cols) {
      if (col %in% colnames(results)) {
        results <- results %>% rename(gene = col)
        break
      }
    }
    
    # If still no gene column, use first column
    if (!"gene" %in% colnames(results)) {
      results <- results %>% rename(gene = 1)
    }
  }
  
  # Get genes to display
  selected_genes <- NULL
  
  # Option 1: Use custom genes from file
  if (!is.null(CUSTOM_GENE_FILE)) {
    selected_genes <- load_genes_from_file(CUSTOM_GENE_FILE)
    cat("Loaded", length(selected_genes), "genes from file\n")
  }
  
  # Option 2: Use custom gene list
  if (is.null(selected_genes) && USE_CUSTOM_GENES && length(CUSTOM_GENE_LIST) > 0) {
    selected_genes <- CUSTOM_GENE_LIST
    cat("Using", length(selected_genes), "genes from custom list\n")
  }
  
  # Option 3: Use top differentially expressed genes
  if (is.null(selected_genes)) {
    # Look for log2FoldChange or similar column
    log2fc_col <- NULL
    for (col in c("log2FoldChange", "log2FC", "logFC")) {
      if (col %in% colnames(results)) {
        log2fc_col <- col
        break
      }
    }
    
    # Look for adjusted p-value column
    padj_col <- NULL
    for (col in c("padj", "FDR", "adj.P.Val", "qvalue")) {
      if (col %in% colnames(results)) {
        padj_col <- col
        break
      }
    }
    
    # Filter and select top genes
    if (!is.null(log2fc_col) && !is.null(padj_col)) {
      # Filter significant genes
      sig_genes <- results %>%
        filter(!!sym(padj_col) < 0.001) %>%
        mutate(abs_fc = abs(!!sym(log2fc_col)))
      
      # Get top up/down regulated genes
      top_genes <- sig_genes %>%
        arrange(desc(abs_fc)) %>%
        head(TOP_GENES_COUNT)
      
      selected_genes <- top_genes$gene
      cat("Selected top", length(selected_genes), "differentially expressed genes\n")
    } else {
      # If no log2FC/padj columns, just take top genes by row variance
      cat("Could not find log2FoldChange or padj columns, selecting random genes\n")
      selected_genes <- results$gene[1:min(TOP_GENES_COUNT, nrow(results))]
    }
  }
  
  # Check if any of the selected genes are in the results
  genes_found <- selected_genes[selected_genes %in% results$gene]
  if (length(genes_found) == 0) {
    stop("Error: None of the specified genes were found in the results")
  }
  
  # Filter results to selected genes
  filtered_results <- results %>%
    filter(gene %in% genes_found)
  
  cat("Creating heatmap with", nrow(filtered_results), "genes\n")
  
  # Get expression data columns (all numeric columns except statistical columns)
  stat_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", 
                 "log2FC", "logFC", "FDR", "adj.P.Val", "qvalue", "gene")
  expr_cols <- colnames(results)[!colnames(results) %in% stat_cols]
 
   # Create log2FC annotation for genes
  # First, identify which log2FC column name is used in your results
  log2fc_col <- NULL
  for (col_name in c("log2FoldChange", "log2FC", "logFC")) {
    if (col_name %in% colnames(filtered_results)) {
      log2fc_col <- col_name
      break
    }
  }
  
  # Create the row annotation data frame
  row_annotation <- filtered_results %>%
    select(gene, !!sym(log2fc_col)) %>%
    rename(log2FC = !!sym(log2fc_col)) %>%
    tibble::column_to_rownames("gene")
  
  # Define custom colors for the log2FC annotation bar
  # This creates a blue-white-red gradient for negative to positive values
  annotation_colors <- list(
    log2FC = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  
  # If we have expression data columns
  if (length(expr_cols) > 0) {
    # Extract expression data
    expr_data <- filtered_results %>%
      select(gene, all_of(expr_cols)) %>%
      tibble::column_to_rownames("gene")
    
    # Create heatmap using expression data
    cat("Generating expression heatmap...\n")
    
    # Create color map
    if (COLOR_SCHEME == "viridis") {
      colors <- viridisLite::viridis(100)
    } else {
      colors <- colorRampPalette(rev(brewer.pal(11, COLOR_SCHEME)))(100)
    }
    
    # Generate heatmap
    pheatmap(
      expr_data,
      scale = SCALE_METHOD,
      cluster_rows = CLUSTER_ROWS,
      cluster_cols = CLUSTER_COLS,
      show_rownames = SHOW_GENE_NAMES,
      show_colnames = SHOW_SAMPLE_NAMES,
      annotation_row = row_annotation,
      #annotation_colors = annotation_colors,
      # Control sample name visibility
      main = paste0("Gene Expression Heatmap (", length(genes_found), " genes)"),
      color = colors,
      fontsize_row = 10,
      angle_col = SAMPLE_NAME_ANGLE,          # Control label angle
      filename = OUTPUT_FILE,
      width = 8,
      height = 6
    )
    
    # Create a second version with different settings if needed
    # Uncomment and modify as needed
    #pheatmap(
    #  expr_data,
    #  scale = SCALE_METHOD,
    #  cluster_rows = CLUSTER_ROWS,
    #  cluster_cols = CLUSTER_COLS,
    #  show_rownames = SHOW_GENE_NAMES,
    #  show_colnames = !SHOW_SAMPLE_NAMES,     # Opposite of first version
    #  main = paste0("Gene Expression Heatmap (", length(genes_found), " genes)"),
    #  color = colors,
    #  fontsize_row = 10,
    #  angle_col = SAMPLE_NAME_ANGLE,
    #  filename = gsub(".png", "_alt.png", OUTPUT_FILE),  # Different filename
    #  width = 8,
    #  height = 6
    #)
    
  } else {
    # If no expression data, create heatmap using log2FC
    cat("No expression data found, creating heatmap using log2FoldChange...\n")
    
    # Find log2FC column
    log2fc_col <- NULL
    for (col in c("log2FoldChange", "log2FC", "logFC")) {
      if (col %in% colnames(filtered_results)) {
        log2fc_col <- col
        break
      }
    }
    
    # If log2FC column exists
    if (!is.null(log2fc_col)) {
      # Create simple heatmap
      heatmap_data <- filtered_results %>%
        select(gene, !!sym(log2fc_col)) %>%
        tibble::column_to_rownames("gene")
      
      pheatmap(
        heatmap_data,
        cluster_rows = CLUSTER_ROWS,
        cluster_cols = FALSE,
        show_rownames = SHOW_GENE_NAMES,
        show_colnames = SHOW_SAMPLE_NAMES,      # Control sample name visibility
        main = paste0("Log2 Fold Change (", length(genes_found), " genes)"),
        color = colorRampPalette(c("blue", "white", "red"))(100),
        fontsize_row = 10,
        angle_col = SAMPLE_NAME_ANGLE,          # Control label angle
        filename = OUTPUT_FILE,
        width = 8,
        height = 6
      )
    } else {
      stop("Error: No expression data or log2FoldChange column found")
    }
  }
  
  cat("Heatmap saved to:", OUTPUT_FILE, "\n")
}

# Execute the main function
create_heatmap()
