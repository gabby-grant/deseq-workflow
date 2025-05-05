# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)  # For column_to_rownames function
library(pheatmap)
library(RColorBrewer)
library(viridis)

###### GENE LIST OPTIONS - ADD THIS SECTION ######
# Set to TRUE to use a custom gene list instead of top DEGs
USE_CUSTOM_GENES <- TRUE


# Option 1: Specify genes directly in the script (comma-separated)
CUSTOM_GENE_LIST <- c("TMEM184A", "PDE1C", "GPR146", "FAM110B")

# Option 2: Read genes from a file (one gene per line)
# Set to NULL to disable, or provide file path like "my_genes.txt"
GENE_LIST_FILE <- NULL

# Function to load genes from a file
load_genes_from_file <- function(file_path) {
  if (!is.null(file_path) && file.exists(file_path)) {
    genes <- readLines(file_path)
    genes <- genes[genes != ""]  # Remove empty lines
    return(genes)
  }
  return(NULL)
}

# Get custom genes based on provided options
get_custom_genes <- function() {
  # First check if there's a file to read from
  if (!is.null(GENE_LIST_FILE)) {
    genes <- load_genes_from_file(GENE_LIST_FILE)
    if (!is.null(genes) && length(genes) > 0) {
      cat("Loaded", length(genes), "genes from file:", GENE_LIST_FILE, "\n")
      return(genes)
    }
  }
  
  # If no file or file reading failed, use the direct list
  if (USE_CUSTOM_GENES && length(CUSTOM_GENE_LIST) > 0) {
    cat("Using", length(CUSTOM_GENE_LIST), "genes from custom list\n")
    return(CUSTOM_GENE_LIST)
  }
  
  # If no custom genes defined, return NULL to use automatic selection
  return(NULL)
}
################################################

# Read your DESeq2 results file
deseq_results <- read_csv("uterus_deseq2_results.csv")

# Name first column of GEM to "gene_name" for clarity
names(deseq_results)[1] <- "gene_name"

# Create a new column to categorize genes
deseq_results <- deseq_results %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      TRUE ~ "Not significant"
    )
  )

# Filter for significant genes only
sig_genes <- deseq_results %>%
  filter(significance != "Not significant")

# Print the number of significant genes
print(paste("Number of significantly differentially expressed genes:", nrow(sig_genes)))

###### MODIFY THIS SECTION TO USE CUSTOM GENES IF PROVIDED ######
# Get custom genes if specified
custom_genes <- get_custom_genes()

# Select genes for the heatmap
if (!is.null(custom_genes)) {
  # Use custom gene list
  # Filter to only include genes that exist in our dataset
  existing_genes <- custom_genes[custom_genes %in% deseq_results$gene_name]
  
  if (length(existing_genes) == 0) {
    stop("Error: None of the specified genes were found in the results.")
  }
  
  # Report on missing genes
  missing_genes <- setdiff(custom_genes, deseq_results$gene_name)
  if (length(missing_genes) > 0) {
    warning("The following genes were not found in the results: ", 
            paste(missing_genes, collapse=", "))
  }
  
  # Select these genes from the results
  top_genes <- deseq_results %>%
    filter(gene_name %in% existing_genes)
  
  print(paste("Using", nrow(top_genes), "custom genes for heatmap visualization"))
  
} else {
  # Use automatic selection of top DEGs as before
  # Modified: Only display top 10 genes
  max_genes_to_display <- 20  # Changed from 50 to 10
  
  # Get top 5 up-regulated and top 5 down-regulated genes
  top_up <- sig_genes %>%
    filter(significance == "Up-regulated") %>%
    arrange(desc(log2FoldChange)) %>%
    head(max_genes_to_display/2)
  
  top_down <- sig_genes %>%
    filter(significance == "Down-regulated") %>%
    arrange(log2FoldChange) %>%
    head(max_genes_to_display/2)
  
  # Combine the top genes
  top_genes <- bind_rows(top_up, top_down) %>%
    arrange(desc(log2FoldChange))
  
  print(paste("Selected", nrow(top_genes), "genes for heatmap visualization"))
}
################################################
# Create sample annotation dataframe
create_sample_annotation <- function(sample_names) {
  # Create data frame for sample annotation
  sample_annotation <- data.frame(
    row.names = sample_names,
    Sample_Type = rep("Unknown", length(sample_names))
  )
  
  # Determine sample type based on naming patterns
  # Pattern recognition for GTEX (normal) vs TCGA (tumor) samples
  for (i in 1:length(sample_names)) {
    sample <- sample_names[i]
    if (grepl("GTEX|gtex|[A-Z]+N$", sample)) {  # Normal samples (GTEX or ending with N)
      sample_annotation$Sample_Type[i] <- "Normal"
    } else if (grepl("TCGA|tcga|[A-Z]+T$", sample)) {  # Tumor samples (TCGA or ending with T)
      sample_annotation$Sample_Type[i] <- "Tumor"
    }
  }
  
  return(sample_annotation)
}

# Extract expression data columns for the heatmap
sample_columns <- colnames(deseq_results)[!colnames(deseq_results) %in% 
                                            c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significance")]
# Add this after you've identified your sample columns and before creating the heatmap



# If no sample columns are found, print a message
if(length(sample_columns) == 0) {
  print("No sample expression columns found in the results file.")
  print("Creating a simplified heatmap using only log2FoldChange values...")
  
  # Create a matrix for the heatmap with just gene names and fold changes
  heatmap_data <- top_genes %>%
    select(gene_name, log2FoldChange) %>%
    mutate(Direction = ifelse(log2FoldChange > 0, "Up-regulated", "Down-regulated")) %>%
    column_to_rownames("gene_name")
  
  # Create the simplified heatmap
  simple_heatmap <- pheatmap(
    heatmap_data["log2FoldChange"], 
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    annotation_row = data.frame(Direction = heatmap_data$Direction, row.names = rownames(heatmap_data)),
    annotation_colors = annotation_colors,
    main = ifelse(!is.null(custom_genes), "Custom Gene Selection (Log2 Fold Change)", "Top 20 Differentially Expressed Genes (Log2 Fold Change)"),
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = ifelse(!is.null(custom_genes), "custom_genes_heatmap.png", "top10_deg_heatmap.png"),
    width = 8,
    height = 6
  )
} else {
  # Extract expression data for the top genes
  expr_data <- top_genes %>%
    select(gene_name, all_of(sample_columns)) %>%
    column_to_rownames("gene_name")
  
  
  # Create sample annotation
  sample_annotation <- create_sample_annotation(colnames(expr_data))
  
  # Sort samples by type
  sample_order <- order(sample_annotation$Sample_Type)
  expr_data_sorted <- expr_data[, sample_order]
  sample_annotation_sorted <- sample_annotation[sample_order, , drop=FALSE]
  
  # Find positions for gaps between groups
  group_counts <- table(sample_annotation_sorted$Sample_Type)
  gaps_col <- cumsum(group_counts)[-length(group_counts)]
  
  # Z-score normalize the sorted data
  expr_data_z_sorted <- t(scale(t(expr_data_sorted)))
  
  
  # Define colors for sample types
  annotation_colors <- list(
    Sample_Type = c("Normal" = "forestgreen", "Tumor" = "gray"),
    Direction = c("Up-regulated" = "red", "Down-regulated" = "blue")
  )
  
  # Create annotation for genes
  gene_annotation <- data.frame(
    row.names = rownames(expr_data),
    Direction = ifelse(top_genes$log2FoldChange > 0, "Up-regulated", "Down-regulated"),
    Log2FC = round(top_genes$log2FoldChange, 2)
  )
  
  # Z-score normalization for better visualization
  expr_data_z <- t(scale(t(expr_data)))
  
  # Adjust the title and output filename based on whether using custom genes or not
  title_text <- ifelse(!is.null(custom_genes), 
                       paste0("Custom Gene Selection (", nrow(top_genes), " genes, Z-score normalized)"), 
                       "Top 20 Differentially Expressed Genes (Z-score normalized)")
  
  output_file <- ifelse(!is.null(custom_genes), "custom_genes_heatmap.png", "top10_deg_heatmap.png")
  
  # Create the heatmap
  heatmap_result <- pheatmap(
    expr_data_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = gene_annotation,
    annotation_colors = annotation_colors,
    main = title_text,
    color = viridis(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = output_file,
    width = 8,
    height = 6
  )
  # Create a second heatmap with different color scheme
  heatmap_alt <- pheatmap(
    expr_data_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = gene_annotation,
    annotation_colors = annotation_colors,
    main = title_text,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = gsub(".png", "_alt.png", output_file),
    width = 8,
    height = 6
  )
  # Create the heatmap with sorted and grouped samples
  heatmap_grouped <- pheatmap(
    expr_data_z_sorted,
    cluster_rows = TRUE,
    cluster_cols = FALSE,          # Turn off column clustering to preserve sorting
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_row = gene_annotation,
    annotation_col = sample_annotation_sorted,
    annotation_colors = annotation_colors,
    #gaps_col = gaps_col,           # Add visual separation between groups
    main = "Gene Expression Heatmap (Samples Grouped by Type)",
    color = viridis(100),
    fontsize_row = 10,
    angle_col = 45,
    filename = "grouped_samples_heatmap.png",
    width = 10,                    # Increased width for better visibility
    height = 6
  )
  
  
}

# Create a heatmap focusing on the log2FoldChange and adjusted p-values
stat_matrix <- top_genes %>%
  select(gene_name, log2FoldChange, padj) %>%
  mutate(neglog10padj = -log10(padj)) %>%
  select(gene_name, log2FoldChange, neglog10padj) %>%
  column_to_rownames("gene_name")

# Normalize the statistical values for better visualization
stat_matrix_z <- as.data.frame(scale(stat_matrix))

# Create the statistical heatmap
stat_title <- ifelse(!is.null(custom_genes), 
                     paste0("Statistical Significance of Custom Gene Selection (", nrow(top_genes), " genes)"), 
                     "Statistical Significance of Top 20 Differentially Expressed Genes")

stat_file <- ifelse(!is.null(custom_genes), "custom_genes_stats_heatmap.png", "top10_deg_stats_heatmap.png")

stat_heatmap <- pheatmap(
  stat_matrix_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  main = stat_title,
  color = colorRampPalette(c("white", "orange", "purple"))(100),
  fontsize_row = 10,  # Increased font size for better readability with fewer genes
  angle_col = 0,
  filename = stat_file,
  width = 8,
  height = 6
)

print("Heatmap generation complete. Check for output PNG files.")
