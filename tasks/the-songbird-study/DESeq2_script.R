# Description: This script performs DESeq2 analysis on RNA-seq count data.

# Load the necessary libraries
library(DESeq2)
library(EnhancedVolcano)

# Read input arguments
args <- commandArgs(TRUE)
input_file_path <- args[2]
output_file_path <- args[3]

# Read the CSV file
tryCatch({
  data <- read.csv(input_file_path, header = TRUE, row.names = 1)
  col_data <- data[, (1:2)]  # Extract relevant columns from the input data
  col_data$study_group <- factor(col_data$study_group)
  col_data$social_settting <- factor(col_data$social_settting)
  
  count_data_transposed <- t(data[, -(1:4)])  # Transpose count data matrix
  
  # Ensure column names of count data match row names of col data
  if (!all(colnames(count_data_transposed) == rownames(col_data))) {
    stop("Column names of count data do not match row names of col data.")
  }
}, error = function(e) {
  stop("Error reading CSV file:", e$message)
})

# Check group condition
gr <- args[1]

if (gr == '0') {
  # For binary comparison between 'S' and 'L' groups
  
  # Check if social_setting has fewer than 2 levels
  if (length(levels(col_data$social_settting)) < 2 & length(levels(col_data$study_group)) < 2) {
    # raise an error if there's only one level for both study_group and social_setting
    stop("There is only one level for both study_group and social_setting. Please check the input data.")
  } else if (length(levels(col_data$social_settting)) < 2) {
    # If there's only one level, exclude social_setting from the design formula
    dds <- DESeqDataSetFromMatrix(countData = count_data_transposed, colData = col_data, design = ~ study_group)
  } else {
    # Otherwise, include social_setting in the design formula
    dds <- DESeqDataSetFromMatrix(countData = count_data_transposed, colData = col_data, design = ~ study_group + social_settting)
  }  
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 20
  dds <- dds[keep,]
  
  # Set reference level for study_group
  dds$study_group <- relevel(dds$study_group, ref = "S")
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get differential expression results
  results <- results(dds, alpha = 0.01, contrast = c("study_group", 'S', 'L'))
  
  # Save summary to a file
  summary_file <- paste0(output_file_path, ".summary")
  sink(summary_file)
  summary(results)
  sink()
  
  # Save results to a CSV file
  write.csv(as.data.frame(results), file = output_file_path)

  # Create and save volcano plot
  volcano_plot <- EnhancedVolcano(results, lab = NA, x = 'log2FoldChange', y = 'padj')
  output_file_path <- gsub(".csv$", "", output_file_path)
  volcano_plot_file <- paste0(output_file_path, "_volcano_plot.png")
  ggsave(volcano_plot_file, plot = volcano_plot, width = 8, height = 6)
} else {
  # For custom comparison between two social settings
  
  base_group <- args[4]
  second_group <- args[5]
  
  # Check if study_group has fewer than 2 levels
  if (length(levels(col_data$social_settting)) < 2 & length(levels(col_data$study_group)) < 2) {
    # raise an error if there's only one level for both study_group and social_setting
    stop("There is only one level for both study_group and social_setting. Please check the input data.")
  } else if (length(levels(col_data$study_group)) < 2) {
    # If there's only one level, exclude study_group from the design formula
    dds <- DESeqDataSetFromMatrix(countData = count_data_transposed, colData = col_data, design = ~ social_settting)
  } else {
    # Otherwise, include study_group in the design formula
    dds <- DESeqDataSetFromMatrix(countData = count_data_transposed, colData = col_data, design = ~ study_group + social_settting)
  }  
  
  # Filter low-count genes
  keep <- rowSums(counts(dds)) >= 20
  dds <- dds[keep,]
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get differential expression results
  results <- results(dds, alpha = 0.01, contrast = c("social_settting", base_group, second_group))
  
  # Save summary to a file
  summary_file <- paste0(output_file_path, ".summary")
  sink(summary_file)
  summary(results)
  sink()
  
  # Save results to a CSV file
  write.csv(as.data.frame(results), file = output_file_path)

  # Create and save volcano plot
  volcano_plot <- EnhancedVolcano(results, lab = NA, x = 'log2FoldChange', y = 'padj')
  output_file_path <- gsub(".csv$", "", output_file_path)
  volcano_plot_file <- paste0(output_file_path, "_volcano_plot.png")
  ggsave(volcano_plot_file, plot = volcano_plot, width = 8, height = 6)
}
