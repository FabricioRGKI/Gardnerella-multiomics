#Merge Modules and Regulated Proteins

# This R script merges a modules assignment file with a regulated proteins file
# by matching the last 5 digits of the locus_tag (since the two files may use 
# different locus tag formats). 
# It performs an inner join on the shortened tag, keeps protein, label, Log2FC, 
# adjusted p-value, module, and subpartition columns, and writes the results to 
# an Excel workbook.
# The workbook contains three sheets: Merged_Data (the joined records), Summary 
# (input row counts, matched records, unique matched tags, and matching rate), 
# and Module_Counts (number of proteins per module).
# Input: a modules CSV (locus_tag, modules, subpartition) and a regulated 
# proteins CSV (locus_tag, protein, label_edit, Log2FC, adj_p_value).
# Output: regulated_modules.xlsx.
# Requirements: R with readr, dplyr, writexl, and stringr.
# Note: update the two input file paths at the top before running.

# ============================================
# Script: Merge Modules and Regulated Proteins
# ============================================

# Clear workspace and set working directory
rm(list = ls())
setwd("your_path")
getwd()

# Load required libraries
library(readr)
library(dplyr)
library(writexl)
library(stringr)

# ============================================
# 1. SPECIFY YOUR FILE PATHS HERE
# ============================================

# Change these paths to your actual file locations
modules_file <- "Gv_modules.csv"
regulated_file <- "Gv_coculture01.csv"

# Change this to your desired output path (or leave as is for current folder)
output_file <- "regulated_modules01.xlsx"

# ============================================
# 2. LOAD FILES
# ============================================

cat("Loading modules file:", modules_file, "\n")
modules <- read_csv(modules_file, show_col_types = FALSE)

cat("Loading regulated proteins file:", regulated_file, "\n")
regulated <- read_csv(regulated_file, show_col_types = FALSE)

# ============================================
# 3. EXTRACT LAST 5 DIGITS OF LOCUS_TAG
# ============================================

modules <- modules %>%
  mutate(locus_tag_short = str_sub(locus_tag, -5))

regulated <- regulated %>%
  mutate(locus_tag_short = str_sub(locus_tag, -5))

# ============================================
# 4. MERGE DATASETS
# ============================================

merged <- inner_join(regulated, modules, 
                     by = "locus_tag_short",
                     suffix = c("_regulated", "_modules"))

# Select and rename columns for output
output <- merged %>%
  select(
    locus_tag_short,
    locus_tag_regulated,
    locus_tag_modules,
    protein,
    label_edit,
    Log2FC,
    adj_p_value,
    modules,
    subpartition
  )

# ============================================
# 5. CREATE SUMMARY
# ============================================

summary_data <- data.frame(
  Metric = c(
    "Total regulated proteins",
    "Total modules entries",
    "Matched records",
    "Unique locus_tags matched",
    "Matching rate (%)"
  ),
  Value = c(
    nrow(regulated),
    nrow(modules),
    nrow(output),
    n_distinct(output$locus_tag_short),
    round(nrow(output) / nrow(regulated) * 100, 1)
  )
)

# Additional summary: counts per module
module_summary <- output %>%
  group_by(modules) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

# ============================================
# 6. SAVE TO EXCEL
# ============================================

write_xlsx(
  list(
    "Merged_Data" = output,
    "Summary" = summary_data,
    "Module_Counts" = module_summary
  ),
  path = output_file
)

cat("\n✓ Done! File saved as:", output_file, "\n")
cat("  - Rows matched:", nrow(output), "\n")
cat("  - Unique locus tags:", n_distinct(output$locus_tag_short), "\n")
