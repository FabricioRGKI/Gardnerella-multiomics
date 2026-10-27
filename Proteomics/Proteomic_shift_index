# PROTEOMIC SHIFT INDEX CALCULATION
# This script processes all CSV files in a folder, parses descriptions,
# adds labels, filters by p-value, and calculates the proteomic shift index
# for upregulated and downregulated proteins separately
# Includes both MEAN and MEDIAN-based calculations

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)

# Set working directory
setwd("your/path")  # <-- THE FOLDER PATH

# Get all CSV files in the folder
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Check if any CSV files exist
if (length(csv_files) == 0) {
  stop("No CSV files found in the current directory!")
}

# Function to process a single file
process_file <- function(file_path) {
  cat("\nProcessing:", basename(file_path), "\n")
  cat("----------------------------------------\n")
  
  # Read the CSV file
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Extract base name for output
  base_name <- tools::file_path_sans_ext(basename(file_path))
  
  # STEP 1: Process the Description column
  df_parsed <- df %>%
    mutate(
      # Extract each label using regex
      locus_tag = str_extract(Description, "\\[locus_tag=([^\\]]+)\\]"),
      # Protein: extract everything after "protein=" until "] [p" or end
      protein = str_extract(Description, "(?<=\\[protein=).*?(?=\\] \\[p|\\]$)"),
      protein_id = str_extract(Description, "\\[protein_id=([^\\]]+)\\]"),
      location = str_extract(Description, "\\[location=([^\\]]+)\\]"),
      gbkey = str_extract(Description, "\\[gbkey=([^\\]]+)\\]"),
      strain = str_extract(Description, "\\[strain=([^\\]]+)\\]"),
      gene = str_extract(Description, "\\[gene=([^\\]]+)\\]")
    ) %>%
    # Clean up extracted values (remove brackets and labels)
    mutate(across(c(locus_tag, protein_id, location, gbkey, strain, gene),
                  ~ str_replace(.x, "\\[.*?=", "") %>% str_remove("\\]$")),
           # Clean protein separately (already extracted without brackets)
           protein = str_trim(protein))
  
  # Write intermediate output CSV (optional - comment out if not needed)
  temp_file <- paste0(base_name, "_desc_filtered.csv")
  write.csv(df_parsed, temp_file, row.names = FALSE)
  cat("  - Created description filtered file:", temp_file, "\n")
  
  # STEP 2: Create the label column
  df_labeled <- df_parsed %>%
    mutate(
      label = case_when(
        !is.na(gene) & gene != "" ~ gene,
        (is.na(gene) | gene == "") & protein != "hypothetical protein" ~ protein,
        (is.na(gene) | gene == "") & protein == "hypothetical protein" ~ locus_tag,
        TRUE ~ NA_character_
      )
    )
  
  # Write labeled output CSV
  labeled_file <- paste0(base_name, "_labeled.csv")
  write.csv(df_labeled, labeled_file, row.names = FALSE)
  cat("  - Created labeled file:", labeled_file, "\n")
  
  # STEP 3: Filter rows with p_value <= 0.05 (remove rows with p_value > 0.05)
  # Check if p_value column exists
  if (!"p_value" %in% colnames(df_labeled)) {
    cat("  WARNING: 'p_value' column not found! Checking for alternative names...\n")
    # Try to find p-value column with alternative names
    pval_cols <- grep("pvalue|p.value|Pvalue|P.value|p_val", colnames(df_labeled), 
                      ignore.case = TRUE, value = TRUE)
    if (length(pval_cols) > 0) {
      pval_col <- pval_cols[1]
      cat("  Using column:", pval_col, "\n")
      df_filtered <- df_labeled %>% filter(!!sym(pval_col) <= 0.05)
    } else {
      cat("  ERROR: No p-value column found. Skipping filtering step.\n")
      df_filtered <- df_labeled
    }
  } else {
    df_filtered <- df_labeled %>% filter(p_value <= 0.05)
  }
  
  cat("  - Removed rows with p-value > 0.05\n")
  cat("  - Remaining rows:", nrow(df_filtered), "\n")
  
  # STEP 4: Divide data into upregulated and downregulated based on Log2FC
  # Check if Log2FC column exists
  if (!"Log2FC" %in% colnames(df_filtered)) {
    cat("  WARNING: 'Log2FC' column not found! Checking for alternative names...\n")
    # Try to find Log2FC column with alternative names
    log2fc_cols <- grep("log2fc|log2FC|Log2FC|logFC|LogFC", colnames(df_filtered), 
                        ignore.case = TRUE, value = TRUE)
    if (length(log2fc_cols) > 0) {
      log2fc_col <- log2fc_cols[1]
      cat("  Using column:", log2fc_col, "\n")
      log2fc_col_name <- log2fc_col
    } else {
      cat("  ERROR: No Log2FC column found! Skipping PSI calculation.\n")
      return(NULL)
    }
  } else {
    log2fc_col_name <- "Log2FC"
  }
  
  # Split the data using the identified Log2FC column
  upregulated <- df_filtered %>% filter(!!sym(log2fc_col_name) > 0)
  downregulated <- df_filtered %>% filter(!!sym(log2fc_col_name) < 0)
  
  cat("  - Upregulated proteins:", nrow(upregulated), "\n")
  cat("  - Downregulated proteins:", nrow(downregulated), "\n")
  
  # STEP 5: Calculate MEAN and MEDIAN of Log2FC for each group
  # MEAN calculations
  mean_up <- if(nrow(upregulated) > 0) mean(upregulated[[log2fc_col_name]], na.rm = TRUE) else 0
  mean_down <- if(nrow(downregulated) > 0) mean(downregulated[[log2fc_col_name]], na.rm = TRUE) else 0
  
  # MEDIAN calculations
  median_up <- if(nrow(upregulated) > 0) median(upregulated[[log2fc_col_name]], na.rm = TRUE) else 0
  median_down <- if(nrow(downregulated) > 0) median(downregulated[[log2fc_col_name]], na.rm = TRUE) else 0
  
  cat("  - Mean Log2FC (Up):", round(mean_up, 4), "\n")
  cat("  - Mean Log2FC (Down):", round(mean_down, 4), "\n")
  cat("  - Median Log2FC (Up):", round(median_up, 4), "\n")
  cat("  - Median Log2FC (Down):", round(median_down, 4), "\n")
  
  # STEP 6: Count rows for each group
  n_up <- nrow(upregulated)
  n_down <- nrow(downregulated)
  n_total <- n_up + n_down
  
  cat("  - Total proteins:", n_total, "\n")
  
  # STEP 7: Calculate Proteomic Shift Index based on MEAN
  # PSI_up_mean = mean Log2FC of upregulated * (number of upregulated / total significant proteins)
  PSI_up_mean <- if(n_total > 0) mean_up * (n_up / n_total) else 0
  
  # PSI_down_mean = mean Log2FC of downregulated * (number of downregulated / total significant proteins)
  PSI_down_mean <- if(n_total > 0) mean_down * (n_down / n_total) else 0
  
  # STEP 8: Calculate Proteomic Shift Index based on MEDIAN
  # PSI_up_median = median Log2FC of upregulated * (number of upregulated / total significant proteins)
  PSI_up_median <- if(n_total > 0) median_up * (n_up / n_total) else 0
  
  # PSI_down_median = median Log2FC of downregulated * (number of downregulated / total significant proteins)
  PSI_down_median <- if(n_total > 0) median_down * (n_down / n_total) else 0
  
  cat("  - PSI (Mean-based - Upregulated):", round(PSI_up_mean, 4), "\n")
  cat("  - PSI (Mean-based - Downregulated):", round(PSI_down_mean, 4), "\n")
  cat("  - PSI (Median-based - Upregulated):", round(PSI_up_median, 4), "\n")
  cat("  - PSI (Median-based - Downregulated):", round(PSI_down_median, 4), "\n")
  
  # STEP 9: Calculate combined PSI (optional - sum of both directions)
  PSI_total_mean <- PSI_up_mean + PSI_down_mean
  PSI_total_median <- PSI_up_median + PSI_down_median
  
  cat("  - PSI Combined (Mean-based):", round(PSI_total_mean, 4), "\n")
  cat("  - PSI Combined (Median-based):", round(PSI_total_median, 4), "\n")
  
  # STEP 10: Save results to a summary file
  results <- data.frame(
    File = basename(file_path),
    Total_Significant_Proteins = n_total,
    Upregulated_Count = n_up,
    Downregulated_Count = n_down,
    Prop_Up = n_up / n_total,
    Prop_Down = n_down / n_total,
    Mean_Log2FC_Up = mean_up,
    Mean_Log2FC_Down = mean_down,
    Median_Log2FC_Up = median_up,
    Median_Log2FC_Down = median_down,
    PSI_Up_Mean = PSI_up_mean,
    PSI_Down_Mean = PSI_down_mean,
    PSI_Up_Median = PSI_up_median,
    PSI_Down_Median = PSI_down_median,
    PSI_Combined_Mean = PSI_total_mean,
    PSI_Combined_Median = PSI_total_median
  )
  
  return(results)
}

# Process all CSV files and collect results
cat("\n========================================\n")
cat("PROCESSING ALL CSV FILES\n")
cat("========================================\n")

all_results <- list()

for (file in csv_files) {
  result <- process_file(file)
  if (!is.null(result)) {
    all_results[[file]] <- result
  }
}

# Combine all results into a single dataframe
if (length(all_results) > 0) {
  summary_df <- do.call(rbind, all_results)
  
  # Write summary file
  summary_file <- "proteomic_shift_index_summary.csv"
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  cat("\n========================================\n")
  cat("SUMMARY OF RESULTS\n")
  cat("========================================\n")
  print(summary_df)
  
  cat("\n\nSummary saved to:", summary_file, "\n")
