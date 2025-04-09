# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)  # For column_to_rownames function
library(pheatmap)
library(RColorBrewer)
library(viridis)

# Read your DESeq2 results file
deseq_results <- read_csv("cervix_deseq2_results.csv")

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

# Modified: Only display top 10 genes
max_genes_to_display <- 10  # Changed from 50 to 10

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

# Extract expression data columns for the heatmap
sample_columns <- colnames(deseq_results)[!colnames(deseq_results) %in% 
                                            c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significance")]

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
    annotation_colors = list(Direction = c("Up-regulated" = "red", "Down-regulated" = "blue")),
    main = "Top 10 Differentially Expressed Genes (Log2 Fold Change)",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = "top10_deg_heatmap.png",
    width = 8,
    height = 6
  )
} else {
  # Extract expression data for the top genes
  expr_data <- top_genes %>%
    select(gene_name, all_of(sample_columns)) %>%
    column_to_rownames("gene_name")
  
  # Create annotation for genes
  gene_annotation <- data.frame(
    row.names = rownames(expr_data),
    Direction = ifelse(top_genes$log2FoldChange > 0, "Up-regulated", "Down-regulated"),
    Log2FC = round(top_genes$log2FoldChange, 2)
  )
  
  # Z-score normalization for better visualization
  expr_data_z <- t(scale(t(expr_data)))
  
  # Create the heatmap
  heatmap_result <- pheatmap(
    expr_data_z,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_row = gene_annotation,
    annotation_colors = list(Direction = c("Up-regulated" = "red", "Down-regulated" = "blue")),
    main = "Top 10 Differentially Expressed Genes (Z-score normalized)",
    color = viridis(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = "top10_deg_heatmap.png",
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
    annotation_colors = list(Direction = c("Up-regulated" = "red", "Down-regulated" = "blue")),
    main = "Top 10 Differentially Expressed Genes (Z-score normalized)",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize_row = 10,  # Increased font size for better readability with fewer genes
    angle_col = 45,
    filename = "top10_deg_heatmap_alt.png",
    width = 8,
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
stat_heatmap <- pheatmap(
  stat_matrix_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  main = "Statistical Significance of Top 10 Differentially Expressed Genes",
  color = colorRampPalette(c("white", "orange", "purple"))(100),
  fontsize_row = 10,  # Increased font size for better readability with fewer genes
  angle_col = 0,
  filename = "top10_deg_stats_heatmap.png",
  width = 8,
  height = 6
)

print("Heatmap generation complete. Check for output PNG files.")
