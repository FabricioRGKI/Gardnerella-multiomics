# ----------------------------------------------------------------------
# Script: targeted_stats_matched_pairs_with_CSV.R
# Purpose: Perform specific Welch t‑tests and generate a
# CSV source file with full statistics.
# ----------------------------------------------------------------------
rm(list = ls())
setwd("your/working/path")
getwd()

library(tidyverse)

# ----------------------------------------------------------------------
# 1. Read CSV file
# ----------------------------------------------------------------------
csv_file <- "VaBaMet_matrix.csv"   # <-- CHANGE TO YOUR FILE
if (!file.exists(csv_file)) {
  stop("File not found: ", csv_file)
}
data <- read_csv(csv_file, show_col_types = FALSE)

# ----------------------------------------------------------------------
# 2. Keep only rows with non‑missing "name"
# ----------------------------------------------------------------------
data <- data %>% filter(!is.na(name))

# ----------------------------------------------------------------------
# 3. Identify technical replicate columns
# ----------------------------------------------------------------------
tech_cols <- grep("^\\d+\\.[A-Za-z]+\\.\\d+$", names(data), value = TRUE)
if (length(tech_cols) == 0) {
  stop("No columns matching pattern '##.XX.##' (technical replicates) found.")
}

# ----------------------------------------------------------------------
# 4. Reshape to long format (keep 'short name')
# ----------------------------------------------------------------------
long_data <- data %>%
  select(name, `short name`, all_of(tech_cols)) %>%
  pivot_longer(cols = all_of(tech_cols), names_to = "colname", values_to = "intensity") %>%
  separate(colname, into = c("replicate", "condition", "tech_rep"), sep = "\\.", remove = TRUE) %>%
  select(-tech_rep) %>%
  mutate(replicate = as.integer(replicate))

# ----------------------------------------------------------------------
# 5. Filter for MbC
# ----------------------------------------------------------------------
#For (e.g.) Methyl beta carboline
target_metabolite <- "MbC"   # Adjust if needed

mbc_data <- long_data %>%
  filter(`short name` == target_metabolite)

if (nrow(mbc_data) == 0) {
  stop("No data found for short name '", target_metabolite, "'.")
}

# ----------------------------------------------------------------------
# 6. Define the specific comparisons
# ----------------------------------------------------------------------
#For Pb coculture comparisons
comparisons <- list(
  c("Lc", "LcPb"),
  c("Li", "LiPb"),
  c("Gv", "GvPb"),
  c("Fv", "PbFv"),
  c("Pb", "LcPb"),
  c("Pb", "LiPb"),
  c("Pb", "GvPb"),
  c("Pb", "PbFv")
)

#For Gv coculture comparisons
comparisons <- list(
  c("Lc", "LcGv"),
  c("Li", "LiGv"),
  c("Pb", "GvPb"),
  c("Fv", "GvFv"),
  c("Gv", "LcGv"),
  c("Gv", "LiGv"),
  c("Gv", "GvPb"),
  c("Gv", "GvFv")
)
# ----------------------------------------------------------------------
# 7. Perform Welch t‑tests for each defined pair
# ----------------------------------------------------------------------
results <- map_dfr(comparisons, function(pair) {
  cond1 <- pair[1]
  cond2 <- pair[2]
  
  # Check that both conditions exist
  if (!(cond1 %in% mbc_data$condition) | !(cond2 %in% mbc_data$condition)) {
    return(tibble(
      group1 = cond1,
      group2 = cond2,
      n1 = NA_integer_,
      n2 = NA_integer_,
      mean1 = NA_real_,
      mean2 = NA_real_,
      mean_diff = NA_real_,
      p.value = NA_real_,
      t.statistic = NA_real_,
      df = NA_real_,
      note = "Condition not found"
    ))
  }
  
  x <- mbc_data %>% filter(condition == cond1) %>% pull(intensity)
  y <- mbc_data %>% filter(condition == cond2) %>% pull(intensity)
  
  test <- tryCatch(
    t.test(x, y, var.equal = FALSE, paired = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(test)) {
    tibble(
      group1 = cond1,
      group2 = cond2,
      n1 = length(x),
      n2 = length(y),
      mean1 = mean(x, na.rm = TRUE),
      mean2 = mean(y, na.rm = TRUE),
      mean_diff = mean1 - mean2,
      p.value = NA_real_,
      t.statistic = NA_real_,
      df = NA_real_,
      note = "Test failed"
    )
  } else {
    tibble(
      group1 = cond1,
      group2 = cond2,
      n1 = length(x),
      n2 = length(y),
      mean1 = mean(x, na.rm = TRUE),
      mean2 = mean(y, na.rm = TRUE),
      mean_diff = mean1 - mean2,
      p.value = test$p.value,
      t.statistic = test$statistic,
      df = test$parameter,
      note = ""
    )
  }
})

# ----------------------------------------------------------------------
# 8. Apply FDR correction
# ----------------------------------------------------------------------
results <- results %>%
  filter(!is.na(p.value)) %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),
    significance = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  arrange(p.adj)

# ----------------------------------------------------------------------
# 9. Prepare  output CSV
# ----------------------------------------------------------------------
# Add metadata columns
output_table <- results %>%
  mutate(
    Metabolite = target_metabolite,
    Filter = "Conditions containing 'Pb' or two-letter base conditions",
    Test = "Welch's two-tailed t-test (unequal variances)",
    Correction = paste0("Benjamini-Hochberg FDR across ", nrow(results), " tests"),
    Comparison = paste(group1, "vs", group2),
    .before = 1
  ) %>%
  select(
    Metabolite, Filter, Comparison,
    group1, group2,
    n1, n2,
    mean1, mean2, mean_diff,
    t.statistic, df,
    p.value, p.adj, significance,
    Test, Correction,
    note
  )

# Round numeric columns for cleaner output
output_table <- output_table %>%
  mutate(across(c(mean1, mean2, mean_diff, t.statistic, df), ~ round(.x, 3)),
         across(c(p.value, p.adj), ~ signif(.x, 3)))

# ----------------------------------------------------------------------
# 10. Save CSV
# ----------------------------------------------------------------------
output_file <- paste0(target_metabolite, "_targeted_stats_source.csv")
write_csv(output_table, output_file)
cat("\nSource data saved to:", output_file, "\n")

# ----------------------------------------------------------------------
# 11. Print summary to console
# ----------------------------------------------------------------------
cat("\n========== Targeted Welch t‑tests for", target_metabolite, "==========\n")
cat("FDR correction (Benjamini–Hochberg) applied across", nrow(results), "tests.\n\n")
print(output_table %>% select(Comparison, n1, n2, mean1, mean2, p.adj, significance))
cat("\n============================================================\n")
