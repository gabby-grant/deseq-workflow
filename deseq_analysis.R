.libPaths(c("~/R_packages", .libPaths()))
library(DESeq2)

# function to extract sub-matrix of counts for each group
subgem <- function(gem, anot, group) {
  datalist = list()
  subanot = subset(anot, Comparison == group)
  for (id in subanot$Sample) {
    ind = which(colnames(gem) == id)
    genes = gem[0]
    exp = gem[,ind]
    datalist[[id]] <- exp
  }
  subcounts = cbind(genes, datalist)
  return(subcounts)
}

# function to extract a subset of the sample annotation matrix for each group
subanot <- function(anot, group) {
  datalist = list()
  print(str(group))
  subanot = subset(anot, Comparison == group)
  print(str(subanot))
  return(subanot)
}

# Modified function to run deseq with tissue type as parameter
run_deseq <- function(counts, annotation, tissue_type = "cervix") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)
  
  # Filter and normalize genes with low total counts across all samples
  dds <- dds[rowSums(counts(dds)) >= 50,]
  dds <- DESeq(dds)
  norm <- fpm(dds)
  
  # Dynamically create condition names based on tissue_type
  condition_gtex <- paste0(tissue_type, "_GTEX")
  condition_tcga <- paste0(tissue_type, "_TCGA")
  
  # Sort the columns in the FPM data frame
  conditionA <- which(annotation$Group == condition_gtex)  
  conditionB <- which(annotation$Group == condition_tcga)
  
  # Check if both conditions exist in the data
  if (length(conditionA) == 0 || length(conditionB) == 0) {
    stop(paste("Error: One or both conditions not found in data:", 
               condition_gtex, "and/or", condition_tcga))
  }
  
  norm <- norm[, c(conditionA, conditionB)]
  print(str(norm))
  
  # Retrieve the results with dynamic contrast
  res <- results(dds, contrast=c("Group", condition_gtex, condition_tcga))
  print(summary(res))
  res <- cbind(as.data.frame(res), norm)
  resultsNames(dds)
  return(res)
}

# Modified main function to accept tissue type parameter
# Modified main function with better file output handling
main <- function(countfile, anotfile, outfile){
  # Create proper file path
  outname = outfile
  
  # Add file extension if missing
  if(!grepl("\\.csv$", outname)) {
    outname = paste0(outname, ".csv")
  }
  
  # Print file information for debugging
  cat("Input count file:", countfile, "\n")
  cat("Input annotation file:", anotfile, "\n")
  cat("Output will be written to:", outname, "\n")
  
  # Read input files
  counts = read.delim(countfile, sep=',', header=TRUE, row.names='Hugo_Symbol')
  samples = read.delim(anotfile, sep=',', row.names = NULL, check.names=FALSE)
  groups = unique(samples$Comparison)
  
  all_results = list()
  
  for (t in groups){
    cat("Processing group:", t, "\n")
    subcounts = subgem(counts, samples, t)
    subannotation = subanot(samples, t)
    
    # Use try-catch to handle potential errors
    results = tryCatch({
      run_deseq(subcounts, subannotation)
    }, error = function(e) {
      cat("Error in DESeq2 analysis:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(results)) {
      # Filter and sort results table
      f_results = subset(results, padj < 0.05)
      o_results = f_results[order(f_results$padj),]
      all_results[[t]] = o_results
      
      # Create group-specific output file if multiple groups
      if (length(groups) > 1) {
        group_outname = paste0(sub("\\.csv$", "", outname), "_", t, ".csv")
        cat("Writing results for group", t, "to", group_outname, "\n")
        tryCatch({
          write.csv(o_results, group_outname, row.names = TRUE)
          cat("Successfully wrote", nrow(o_results), "rows to", group_outname, "\n")
        }, error = function(e) {
          cat("Error writing file:", e$message, "\n")
        })
      }
    }
  }
  
  # For single group or combined results
  if (length(groups) == 1 && length(all_results) > 0) {
    tryCatch({
      cat("Writing results to", outname, "\n")
      write.csv(all_results[[1]], outname, row.names = TRUE)
      cat("Successfully wrote", nrow(all_results[[1]]), "rows to", outname, "\n")
    }, error = function(e) {
      cat("Error writing file:", e$message, "\n")
    })
  }
  
  # Return filtered results instead of unfiltered
  if (length(all_results) == 1) {
    return(all_results[[1]])
  } else {
    return(all_results)
  }
}
# Script execution code - add this at the end of your file
if (interactive()) {
  # Set your file paths here
  count_file <- "cervix-gtex-tcga-clean.csv"  # CHANGE THIS to your actual file name
  annotation_file <- "cervix-gtex-tcga.comparison.csv"  # CHANGE THIS to your actual file name
  output_file <- "cervix_deseq2_results.csv"
  
  # Print startup message
  cat("\n========================================\n")
  cat("Starting DESeq2 Analysis\n")
  cat("========================================\n\n")
  
  # Run the analysis
  results <- main(count_file, annotation_file, output_file)
  
  # Show that execution is complete
  cat("\nAnalysis complete! Check your working directory for output files.\n")
  cat("Current working directory: ", getwd(), "\n")
}
