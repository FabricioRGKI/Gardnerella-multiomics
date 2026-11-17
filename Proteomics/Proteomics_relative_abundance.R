rm(list = ls())
setwd("/your/path")
getwd()

# ============================================================================
# Proteomics_relative_abundance.R
# ----------------------------------------------------------------------------
# Intensity-based relative abundance between two species in a coculture.
# - Keeps only the coculture abundance columns (control/monoculture excluded
#   automatically by exact label match).
# - Matches species by [strain=Genus_species_...]
# - NA filter: at most 1 NA across coculture replicates.
# - Half-min imputation on remaining NAs.
# - Per-column normalization.
# - Intensity-based relative abundance per species per replicate.
# - Paired between-species test (Welch + exact permutation; n = 3).
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

proteomics_species_rel_abund <- function(
    file_path,
    mix_condition,
    species_names,
    output_prefix = "species_rel_abund"
) {
  
  # --------------------------------------------------------------------------
  # Logging helpers
  # --------------------------------------------------------------------------
  hr <- function(title) {
    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("  ", title, "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")
  }
  step <- function(n, title) {
    cat("\n", strrep("-", 78), "\n", sep = "")
    cat("  STEP ", n, ": ", title, "\n", sep = "")
    cat(strrep("-", 78), "\n", sep = "")
  }
  show_df <- function(x, n = 6, label = NULL) {
    if (!is.null(label)) cat(label, "\n")
    cat("  dim: ", nrow(x), " rows x ", ncol(x), " cols\n", sep = "")
    if (nrow(x) == 0) { cat("  (empty)\n"); return(invisible()) }
    print(utils::head(as.data.frame(x), n), row.names = FALSE)
    if (nrow(x) > n) cat("  ... (", nrow(x) - n, " more rows)\n", sep = "")
  }
  
  # --------------------------------------------------------------------------
  # STEP 0: Read input
  # --------------------------------------------------------------------------
  step(0, "Read input file")
  cat("  file_path       : ", file_path, "\n", sep = "")
  cat("  mix_condition   : ", mix_condition, "\n", sep = "")
  cat("  species_names   : ", paste(species_names, collapse = ", "), "\n", sep = "")
  
  df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  cat("  Loaded ", nrow(df), " rows x ", ncol(df), " columns\n", sep = "")
  cat("  Column names (first 15):\n")
  cat("    ", paste(head(names(df), 15), collapse = "\n    "), "\n", sep = "")
  if (ncol(df) > 15) cat("    ... (", ncol(df) - 15, " more)\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 1: Keep ONLY coculture columns
  # --------------------------------------------------------------------------
  step(1, "Select coculture columns (exact label match)")
  
  abund_cols_all <- grep("^Abundance:", names(df), value = TRUE)
  cat("  Total 'Abundance:' columns found: ", length(abund_cols_all), "\n", sep = "")
  
  label <- sub("_R[0-9]+.*$", "", sub("^.*:\\s*", "", abund_cols_all))
  cat("  Parsed condition labels (unique):\n")
  print(sort(unique(label)))
  
  mix_cols <- abund_cols_all[label == mix_condition]
  dropped  <- setdiff(abund_cols_all, mix_cols)
  
  cat("\n  Coculture columns KEPT (", length(mix_cols), "):\n", sep = "")
  cat("    ", paste(mix_cols, collapse = "\n    "), "\n", sep = "")
  cat("\n  Other columns DROPPED (", length(dropped), "):\n", sep = "")
  cat("    ", paste(dropped, collapse = "\n    "), "\n", sep = "")
  
  if (length(mix_cols) == 0)
    stop("No abundance columns matched mix_condition = '", mix_condition, "'.")
  
  meta_cols <- intersect(c(
    "Description", "Sum PEP Score", "Coverage [%]", "# Peptides", "# PSMs",
    "# Unique Peptides", "# AAs", "MW [kDa]", "calc. pI",
    "Score MS Amanda 2.0: MS Amanda 2.0",
    "# Peptides (by Search Engine): MS Amanda 2.0",
    "Biological Process", "Cellular Component", "Molecular Function",
    "Pfam IDs", "Entrez Gene ID", "Gene ID", "Gene Symbol"
  ), names(df))
  cat("\n  Metadata columns carried through (", length(meta_cols), "):\n", sep = "")
  cat("    ", paste(meta_cols, collapse = ", "), "\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 1b: Species matching via [strain=...] tag
  # --------------------------------------------------------------------------
  step("1b", "Species matching via [strain=...] tag")
  
  strain_tag <- sub(".*\\[strain=([^]]+)\\].*", "\\1", df$Description)
  strain_tag[!grepl("\\[strain=", df$Description)] <- NA_character_
  norm_strain <- tolower(gsub("_", " ", strain_tag))
  
  cat("  Rows with a [strain=...] tag: ", sum(!is.na(strain_tag)),
      " / ", nrow(df), "\n", sep = "")
  cat("  Distinct strain tags (up to 30):\n")
  strains_all <- sort(unique(strain_tag[!is.na(strain_tag)]))
  cat("    ", paste(head(strains_all, 30), collapse = "\n    "), "\n", sep = "")
  if (length(strains_all) > 30) cat("    ... (", length(strains_all) - 30, " more)\n", sep = "")
  
  species_prefixes <- paste0(tolower(species_names), " ")
  
  is_target <- Reduce(`|`, lapply(species_prefixes, function(p) {
    !is.na(norm_strain) & startsWith(norm_strain, p)
  }))
  
  n_matched <- sum(is_target)
  cat("\n  Rows matching requested species: ", n_matched, "\n", sep = "")
  for (sp in species_names) {
    p <- paste0(tolower(sp), " ")
    cat("    ", sp, ": ", sum(!is.na(norm_strain) & startsWith(norm_strain, p)),
        " rows\n", sep = "")
  }
  
  if (n_matched == 0) {
    strains <- sort(unique(strain_tag[!is.na(strain_tag)]))
    stop("No rows matched the species names.\n",
         "Available strain tags (first 30):\n  ",
         paste(head(strains, 30), collapse = "\n  "))
  }
  
  filtered_df <- df[is_target, c(meta_cols, mix_cols), drop = FALSE]
  cat("\n  Subset to target rows + relevant cols:\n")
  show_df(filtered_df, n = 3)
  
  # --------------------------------------------------------------------------
  # STEP 2: >= 2 unique peptides
  # --------------------------------------------------------------------------
  step(2, "Filter: # Unique Peptides >= 2")
  
  n_before <- nrow(filtered_df)
  filtered_df <- filtered_df[filtered_df[["# Unique Peptides"]] >= 2, , drop = FALSE]
  cat("  Rows before: ", n_before, "\n", sep = "")
  cat("  Rows after : ", nrow(filtered_df), "  (removed ",
      n_before - nrow(filtered_df), ")\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 3: NA filter (<=1 NA across coculture reps)
  # --------------------------------------------------------------------------
  step(3, "Filter: at most 1 NA across coculture replicates")
  
  abundance_cols <- mix_cols
  na_count <- rowSums(is.na(filtered_df[, abundance_cols, drop = FALSE]))
  
  cat("  NA count distribution per row:\n")
  print(table(na_count, useNA = "ifany"))
  
  n_before <- nrow(filtered_df)
  filtered_df <- filtered_df[na_count <= 1, , drop = FALSE]
  cat("  Rows before: ", n_before, "\n", sep = "")
  cat("  Rows after : ", nrow(filtered_df), "  (removed ",
      n_before - nrow(filtered_df), ")\n", sep = "")
  
  if (nrow(filtered_df) == 0)
    stop("No rows survived filtering.")
  
  # --------------------------------------------------------------------------
  # STEP 4: Half-min imputation
  # --------------------------------------------------------------------------
  step(4, "Half-min imputation of remaining NAs")
  
  n_na_before <- sum(is.na(filtered_df[, abundance_cols, drop = FALSE]))
  cat("  NA cells before imputation: ", n_na_before, "\n", sep = "")
  
  replace_na_with_half_min <- function(x) {
    m <- as.matrix(x)
    storage.mode(m) <- "double"
    row_min <- apply(m, 1, min, na.rm = TRUE)
    row_min[!is.finite(row_min)] <- NA_real_
    half <- row_min / 2
    na_idx <- which(is.na(m), arr.ind = TRUE)
    if (nrow(na_idx) > 0) {
      m[na_idx] <- half[na_idx[, 1]]
    }
    as.data.frame(m, check.names = FALSE)
  }
  
  abundance_data <- replace_na_with_half_min(
    filtered_df[, abundance_cols, drop = FALSE]
  )
  
  non_abundance_cols <- setdiff(names(filtered_df), abundance_cols)
  filtered_df <- cbind(
    filtered_df[, non_abundance_cols, drop = FALSE],
    abundance_data
  )
  
  n_na_after <- sum(is.na(filtered_df[, abundance_cols, drop = FALSE]))
  cat("  NA cells after imputation : ", n_na_after, "\n", sep = "")
  cat("  Example (first 3 rows, abundance cols only):\n")
  show_df(filtered_df[1:min(3, nrow(filtered_df)), abundance_cols, drop = FALSE])
  
  # ---- FULL TABLE DUMP TO CSV ----
  step4_file <- paste0(output_prefix, "_", mix_condition, "_step4_imputed.csv")
  write.csv(filtered_df, step4_file, row.names = FALSE)
  cat("  Full filtered_df (after imputation) written to: ", step4_file,
      "  (", nrow(filtered_df), " rows x ", ncol(filtered_df), " cols)\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 5: Normalize each coculture column by its own sum
  # --------------------------------------------------------------------------
  step(5, "Normalize each coculture column by its own sum")
  
  col_sums <- colSums(filtered_df[, abundance_cols, drop = FALSE], na.rm = TRUE)
  cat("  Column sums BEFORE normalization:\n")
  print(col_sums)
  
  normalized_abundance <- as.data.frame(
    sweep(
      as.matrix(filtered_df[, abundance_cols, drop = FALSE]),
      2, col_sums, "/"
    ),
    check.names = FALSE
  )
  
  cat("\n  Column sums AFTER normalization (should be 1):\n")
  print(colSums(normalized_abundance, na.rm = TRUE))
  
  normalized_df <- cbind(
    filtered_df[, non_abundance_cols, drop = FALSE],
    normalized_abundance
  )
  
  cat("\n  First 3 normalized rows:\n")
  show_df(normalized_df[1:min(3, nrow(normalized_df)), abundance_cols, drop = FALSE])
  
  # ---- FULL TABLE DUMP TO CSV ----
  step5_file <- paste0(output_prefix, "_", mix_condition, "_step5_normalized.csv")
  write.csv(normalized_df, step5_file, row.names = FALSE)
  cat("  Full normalized_df written to: ", step5_file,
      "  (", nrow(normalized_df), " rows x ", ncol(normalized_df), " cols)\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 6: Replicate check
  # --------------------------------------------------------------------------
  step(6, "Confirm coculture replicate columns")
  cat("  Found ", length(mix_cols), " coculture replicate(s)\n", sep = "")
  cat("    ", paste(mix_cols, collapse = ", "), "\n", sep = "")
  if (length(mix_cols) != 3)
    warning("Expected 3 coculture replicates, found ", length(mix_cols))
  
  # --------------------------------------------------------------------------
  # STEP 7: Per species, per replicate — intensity sum + row count
  # --------------------------------------------------------------------------
  step(7, "Per species x replicate aggregation (intensity sum + row count)")
  
  # Pre-compute strain tags once, directly on normalized_df (fixes the
  # df-vs-normalized_df logical-index misalignment).
  desc_norm        <- normalized_df$Description
  strain_norm      <- sub(".*\\[strain=([^]]+)\\].*", "\\1", desc_norm)
  strain_norm[!grepl("\\[strain=", desc_norm)] <- NA_character_
  norm_strain_norm <- tolower(gsub("_", " ", strain_norm))
  
  per_species_rep <- do.call(rbind, lapply(species_names, function(sp) {
    sp_prefix <- paste0(tolower(sp), " ")
    hit <- !is.na(norm_strain_norm) & startsWith(norm_strain_norm, sp_prefix)
    
    sub <- normalized_df[hit, mix_cols, drop = FALSE]
    
    # ---- DEBUG PRINT ----
    cat("    [debug] species=", sp,
        "  nrow(sub)=", nrow(sub),
        "  class(sub[[1]])=", class(sub[[1]]),
        "  positives=", sum(!is.na(sub[[1]]) & sub[[1]] > 0),
        "\n", sep = "")
    
    cat("    ", sp, ": ", nrow(sub), " rows contributed\n", sep = "")
    
    if (nrow(sub) == 0) {
      warning("Species '", sp, "' has 0 rows after filtering.")
      return(data.frame(
        species = sp, replicate = mix_cols,
        n_rows = NA_integer_, intensity = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    data.frame(
      species   = sp,
      replicate = mix_cols,
      n_rows    = vapply(sub, function(col) {
        sum(!is.na(col) & col > 0)
      }, integer(1)),
      intensity = vapply(sub, function(col) {
        sum(col, na.rm = TRUE)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }))
  
  cat("\n  Aggregated per species x replicate:\n")
  print(as.data.frame(per_species_rep), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 8: Intensity-based relative abundance per replicate
  # --------------------------------------------------------------------------
  step(8, "Intensity-based relative abundance per replicate")
  
  rel_per_rep <- per_species_rep %>%
    group_by(replicate) %>%
    mutate(
      total_intensity_all_species = sum(intensity, na.rm = TRUE),
      rel_intensity = ifelse(total_intensity_all_species > 0,
                             intensity / total_intensity_all_species,
                             NA_real_)
    ) %>%
    ungroup()
  
  cat("  Per-replicate intensities + relative abundances:\n")
  print(as.data.frame(rel_per_rep), row.names = FALSE)
  
  check_sums <- rel_per_rep %>%
    group_by(replicate) %>%
    summarise(total_rel = sum(rel_intensity, na.rm = TRUE), .groups = "drop")
  cat("\n  Sanity check — sum of rel_intensity per replicate (should be ~1):\n")
  print(as.data.frame(check_sums), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 9: Per-species summary across replicates
  # --------------------------------------------------------------------------
  step(9, "Per-species summary across replicates")
  
  rel_summary <- rel_per_rep %>%
    group_by(species) %>%
    summarise(
      n_replicates         = sum(!is.na(intensity)),
      mean_intensity       = mean(intensity,   na.rm = TRUE),
      median_intensity     = median(intensity, na.rm = TRUE),
      sd_intensity         = sd(intensity,     na.rm = TRUE),
      mean_rel_abundance   = mean(rel_intensity,   na.rm = TRUE),
      median_rel_abundance = median(rel_intensity, na.rm = TRUE),
      sd_rel_abundance     = sd(rel_intensity,     na.rm = TRUE),
      min_rel_abundance    = min(rel_intensity,    na.rm = TRUE),
      max_rel_abundance    = max(rel_intensity,    na.rm = TRUE),
      mean_rows            = mean(n_rows,      na.rm = TRUE),
      median_rows          = median(n_rows,    na.rm = TRUE),
      .groups = "drop"
    )
  
  n_rep <- length(mix_cols)
  if (n_rep >= 2) {
    tcrit <- qt(0.975, df = n_rep - 1)
    se_rel <- rel_summary$sd_rel_abundance / sqrt(n_rep)
    rel_summary$ci95_low  <- rel_summary$mean_rel_abundance - tcrit * se_rel
    rel_summary$ci95_high <- rel_summary$mean_rel_abundance + tcrit * se_rel
  } else {
    rel_summary$ci95_low  <- NA_real_
    rel_summary$ci95_high <- NA_real_
  }
  
  cat("  Summary table (with 95% CI on mean relative abundance):\n")
  print(as.data.frame(rel_summary), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 10: Between-species paired test
  # --------------------------------------------------------------------------
  step(10, "Between-species paired test")
  
  wide_rel <- rel_per_rep %>%
    select(species, replicate, rel_intensity) %>%
    pivot_wider(names_from = species, values_from = rel_intensity)
  
  cat("  Wide format (replicate x species):\n")
  print(as.data.frame(wide_rel), row.names = FALSE)
  
  test_results <- NULL
  if (ncol(wide_rel) == 3) {
    sp1 <- names(wide_rel)[2]; sp2 <- names(wide_rel)[3]
    x1 <- wide_rel[[sp1]]; x2 <- wide_rel[[sp2]]
    d  <- x1 - x2
    n  <- sum(!is.na(d))
    
    cat("\n  Paired differences (", sp1, " - ", sp2, "):\n", sep = "")
    print(d)
    cat("  n = ", n, ", mean(d) = ", mean(d, na.rm = TRUE),
        ", sd(d) = ", sd(d, na.rm = TRUE), "\n", sep = "")
    
    if (n >= 2) {
      tt <- tryCatch(t.test(x1, x2, paired = TRUE), error = function(e) NULL)
      
      signs <- expand.grid(rep(list(c(-1, 1)), n))
      perm_means <- apply(signs, 1, function(s) mean(s * d[!is.na(d)]))
      obs <- mean(d, na.rm = TRUE)
      p_two <- mean(abs(perm_means) >= abs(obs) - 1e-12)
      p_one <- if (obs >= 0) mean(perm_means >= obs - 1e-12)
      else          mean(perm_means <= obs + 1e-12)
      
      dz <- if (sd(d, na.rm = TRUE) > 0) obs / sd(d, na.rm = TRUE) else NA_real_
      se <- sd(d, na.rm = TRUE) / sqrt(n)
      tc <- qt(0.975, df = n - 1)
      
      test_results <- tibble(
        species_1            = sp1,
        species_2            = sp2,
        n_replicates         = n,
        mean_rel_species_1   = mean(x1, na.rm = TRUE),
        mean_rel_species_2   = mean(x2, na.rm = TRUE),
        mean_difference      = obs,
        sd_difference        = sd(d, na.rm = TRUE),
        ci95_diff_low        = obs - tc * se,
        ci95_diff_high       = obs + tc * se,
        cohens_dz            = dz,
        t_statistic          = if (!is.null(tt)) unname(tt$statistic) else NA_real_,
        df                   = if (!is.null(tt)) unname(tt$parameter) else NA_real_,
        p_paired_t_two_sided = if (!is.null(tt)) tt$p.value else NA_real_,
        p_perm_exact_two     = p_two,
        p_perm_exact_one     = p_one,
        direction            = if (obs > 0) paste(sp1, ">", sp2)
        else if (obs < 0) paste(sp2, ">", sp1)
        else "equal"
      )
      
      cat("\n  Between-species test result:\n")
      print(as.data.frame(test_results), row.names = FALSE)
    }
  } else {
    cat("  Skipping test: need exactly 2 species, found ",
        ncol(wide_rel) - 1, ".\n", sep = "")
    warning("Between-species test skipped.")
  }
  
  # --------------------------------------------------------------------------
  # STEP 11: Save outputs
  # --------------------------------------------------------------------------
  step(11, "Save output CSVs")
  
  files <- c(
    per_replicate        = paste0(output_prefix, "_", mix_condition, "_per_replicate.csv"),
    summary              = paste0(output_prefix, "_", mix_condition, "_summary.csv"),
    between_species_test = paste0(output_prefix, "_", mix_condition, "_between_species_test.csv")
  )
  
  write.csv(rel_per_rep, files[["per_replicate"]], row.names = FALSE)
  write.csv(rel_summary, files[["summary"]],       row.names = FALSE)
  if (!is.null(test_results))
    write.csv(test_results, files[["between_species_test"]], row.names = FALSE)
  
  cat("  Files written:\n")
  for (f in files) cat("    ", f, "\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 12: Final summary
  # --------------------------------------------------------------------------
  hr("FINAL SUMMARY")
  
  cat("\n=== INTENSITY-BASED RELATIVE ABUNDANCE PER REPLICATE ===\n")
  print(as.data.frame(rel_per_rep), row.names = FALSE)
  
  cat("\n=== SUMMARY ACROSS ", n_rep, " REPLICATES ===\n", sep = "")
  print(as.data.frame(rel_summary), row.names = FALSE)
  
  if (!is.null(test_results)) {
    cat("\n=== BETWEEN-SPECIES TEST ===\n")
    print(as.data.frame(test_results), row.names = FALSE)
    cat("\nNote: with n = 3, the exact permutation test has a minimum achievable\n",
        "two-sided p-value of 0.25 by construction. Rely on the effect size,\n",
        "95% CI, and paired t-test for interpretation.\n")
  }
  
  invisible(list(
    per_replicate        = rel_per_rep,
    summary              = rel_summary,
    between_species_test = test_results
  ))
}

#RUN CALL EXAMPLES

# LcGv — Lactobacillus crispatus + Gardnerella vaginalis
res_Lc2Gv <- proteomics_species_rel_abund(
  file_path     = "proteins.csv",
  mix_condition = "Lc2Gv",
  species_names = c("Lactobacillus crispatus", "Gardnerella vaginalis"),
  output_prefix = "coculture"
)

# LcFv — Lactobacillus crispatus + Fannyhessea vaginae
res_Lc2Fv <- proteomics_species_rel_abund(
  file_path     = "proteins.csv",
  mix_condition = "Lc2Fv",
  species_names = c("Lactobacillus crispatus", "Fannyhessea vaginae"),
  output_prefix = "coculture"
)

# LcLi — Lactobacillus crispatus + Lactobacillus iners
res_Lc2Li <- proteomics_species_rel_abund(
  file_path     = "proteins.csv",
  mix_condition = "Lc2Li",
  species_names = c("Lactobacillus crispatus", "Lactobacillus iners"),
  output_prefix = "coculture"
)

# LcPb — Lactobacillus crispatus + Prevotella bivia
res_Lc2Pb <- proteomics_species_rel_abund(
  file_path     = "proteins.csv",
  mix_condition = "Lc2Pb",
  species_names = c("Lactobacillus crispatus", "Prevotella bivia"),
  output_prefix = "coculture"
)

______
#FILES WITH SEPARATED SPECIES
rm(list = ls())
setwd("/your/path")
getwd()

# ----------------------------------------------------------------------------
# Reads pairs of split files, combines them, and calculates relative abundance.
# ============================================================================

library(dplyr)
library(tidyr)
library(tibble)

process_combined_coculture <- function(file_1, species_1, file_2, species_2, mix_condition) {
  
  # --------------------------------------------------------------------------
  # Logging helpers
  # --------------------------------------------------------------------------
  hr <- function(title) {
    cat("\n", strrep("=", 78), "\n", sep = "")
    cat("  ", title, "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")
  }
  step <- function(n, title) {
    cat("\n", strrep("-", 78), "\n", sep = "")
    cat("  STEP ", n, ": ", title, "\n", sep = "")
    cat(strrep("-", 78), "\n", sep = "")
  }
  show_df <- function(x, n = 6, label = NULL) {
    if (!is.null(label)) cat(label, "\n")
    cat("  dim: ", nrow(x), " rows x ", ncol(x), " cols\n", sep = "")
    if (nrow(x) == 0) { cat("  (empty)\n"); return(invisible()) }
    print(utils::head(as.data.frame(x), n), row.names = FALSE)
    if (nrow(x) > n) cat("  ... (", nrow(x) - n, " more rows)\n", sep = "")
  }
  
  hr(paste0("Processing: ", mix_condition, " (", species_1, " + ", species_2, ")"))
  
  # --------------------------------------------------------------------------
  # STEP 1: Read both files
  # --------------------------------------------------------------------------
  step(1, "Read input files")
  cat("  file_1: ", file_1, "\n", sep = "")
  cat("  file_2: ", file_2, "\n", sep = "")
  
  df1 <- read.csv(file_1, stringsAsFactors = FALSE, check.names = FALSE)
  df2 <- read.csv(file_2, stringsAsFactors = FALSE, check.names = FALSE)
  
  cat("  df1: ", nrow(df1), " rows x ", ncol(df1), " cols\n", sep = "")
  cat("  df2: ", nrow(df2), " rows x ", ncol(df2), " cols\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 2: Keep only needed columns
  # --------------------------------------------------------------------------
  step(2, "Extract only needed columns (Description, # Unique Peptides, Abundance)")
  
  cols_to_keep <- c("Description", "# Unique Peptides",
                    grep("^Abundance:", names(df1), value = TRUE))
  
  cat("  Columns kept (", length(cols_to_keep), "):\n", sep = "")
  cat("    ", paste(cols_to_keep, collapse = "\n    "), "\n", sep = "")
  
  df1 <- df1[, cols_to_keep, drop = FALSE]
  df2 <- df2[, cols_to_keep, drop = FALSE]
  
  df1$Description <- as.character(df1$Description)
  df2$Description <- as.character(df2$Description)
  
  cat("  df1 after column subset: ", nrow(df1), " rows x ", ncol(df1), " cols\n", sep = "")
  cat("  df2 after column subset: ", nrow(df2), " rows x ", ncol(df2), " cols\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 3: Combine rows
  # --------------------------------------------------------------------------
  step(3, "Combine df1 + df2 with bind_rows()")
  
  df <- bind_rows(df1, df2)
  cat("  Combined: ", nrow(df), " rows x ", ncol(df), " cols\n", sep = "")
  cat("  Source split: df1 = ", nrow(df1), ", df2 = ", nrow(df2),
      "  (sum = ", nrow(df1) + nrow(df2), ")\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 4: Identify coculture columns
  # --------------------------------------------------------------------------
  step(4, "Identify coculture columns")
  
  abund_cols_all <- grep("^Abundance:", names(df), value = TRUE)
  label <- sub("_R[0-9]+.*$", "", sub("^.*:\\s*", "", abund_cols_all))
  
  cat("  Total 'Abundance:' columns: ", length(abund_cols_all), "\n", sep = "")
  cat("  Parsed condition labels (unique):\n")
  print(sort(unique(label)))
  
  mix_cols <- abund_cols_all[label == mix_condition]
  
  if (length(mix_cols) == 0) {
    stop("No abundance columns matched mix_condition '", mix_condition, "'")
  }
  
  cat("\n  mix_cols KEPT (", length(mix_cols), "):\n", sep = "")
  cat("    ", paste(mix_cols, collapse = "\n    "), "\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 5: Match species via [strain=...] tag
  # --------------------------------------------------------------------------
  step(5, "Match species via [strain=...] tag")
  
  strain_tag <- sub(".*\\[strain=([^]]+)\\].*", "\\1", df$Description)
  strain_tag[!grepl("\\[strain=", df$Description)] <- NA_character_
  norm_strain <- tolower(gsub("_", " ", strain_tag))
  
  cat("  Rows with a [strain=...] tag: ", sum(!is.na(strain_tag)),
      " / ", nrow(df), "\n", sep = "")
  
  sp_prefix_1 <- paste0(tolower(species_1), " ")
  sp_prefix_2 <- paste0(tolower(species_2), " ")
  
  is_sp1 <- !is.na(norm_strain) & startsWith(norm_strain, sp_prefix_1)
  is_sp2 <- !is.na(norm_strain) & startsWith(norm_strain, sp_prefix_2)
  
  cat("  Rows matching '", species_1, "': ", sum(is_sp1), "\n", sep = "")
  cat("  Rows matching '", species_2, "': ", sum(is_sp2), "\n", sep = "")
  cat("  Rows matching neither           : ",
      sum(!is_sp1 & !is_sp2), "\n", sep = "")
  
  if (sum(is_sp1) == 0) warning("No rows matched species '", species_1, "'")
  if (sum(is_sp2) == 0) warning("No rows matched species '", species_2, "'")
  
  # --------------------------------------------------------------------------
  # STEP 6: Filter >= 2 unique peptides
  # --------------------------------------------------------------------------
  step(6, "Filter: # Unique Peptides >= 2")
  
  keep_peptides <- df[["# Unique Peptides"]] >= 2
  
  cat("  Rows before: ", nrow(df), "\n", sep = "")
  cat("  Rows passing (>= 2 unique peptides): ", sum(keep_peptides), "\n", sep = "")
  cat("  Rows removed: ", sum(!keep_peptides), "\n", sep = "")
  
  df <- df[keep_peptides, , drop = FALSE]
  is_sp1 <- is_sp1[keep_peptides]
  is_sp2 <- is_sp2[keep_peptides]
  
  cat("  After filter — ", species_1, ": ", sum(is_sp1),
      "  |  ", species_2, ": ", sum(is_sp2), "\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 7: Filter <= 1 NA across mix columns
  # --------------------------------------------------------------------------
  step(7, "Filter: at most 1 NA across coculture replicate columns")
  
  na_count <- rowSums(is.na(df[, mix_cols, drop = FALSE]))
  
  cat("  NA count distribution per row:\n")
  print(table(na_count, useNA = "ifany"))
  
  keep_na <- na_count <= 1
  
  cat("  Rows before: ", nrow(df), "\n", sep = "")
  cat("  Rows passing (<= 1 NA): ", sum(keep_na), "\n", sep = "")
  cat("  Rows removed: ", sum(!keep_na), "\n", sep = "")
  
  df <- df[keep_na, , drop = FALSE]
  is_sp1 <- is_sp1[keep_na]
  is_sp2 <- is_sp2[keep_na]
  
  if (nrow(df) == 0) stop("No rows survived filtering.")
  
  cat("  After filter — ", species_1, ": ", sum(is_sp1),
      "  |  ", species_2, ": ", sum(is_sp2), "\n", sep = "")
  
  # --------------------------------------------------------------------------
  # STEP 8: Half-min imputation
  # --------------------------------------------------------------------------
  step(8, "Half-min imputation of remaining NAs")
  
  m <- as.matrix(df[, mix_cols, drop = FALSE])
  storage.mode(m) <- "double"
  
  na_before <- sum(is.na(m))
  cat("  NA cells before imputation: ", na_before, "\n", sep = "")
  
  row_min <- apply(m, 1, min, na.rm = TRUE)
  row_min[!is.finite(row_min)] <- NA_real_
  half <- row_min / 2
  
  na_idx <- which(is.na(m), arr.ind = TRUE)
  if (nrow(na_idx) > 0) {
    cat("  Cells to impute: ", nrow(na_idx), "\n", sep = "")
    m[na_idx] <- half[na_idx[, 1]]
  } else {
    cat("  No cells to impute.\n")
  }
  
  na_after <- sum(is.na(m))
  cat("  NA cells after imputation : ", na_after, "\n", sep = "")
  
  cat("\n  Min / median / max per column AFTER imputation:\n")
  print(data.frame(
    col    = mix_cols,
    min    = apply(m, 2, min,    na.rm = TRUE),
    median = apply(m, 2, median, na.rm = TRUE),
    max    = apply(m, 2, max,    na.rm = TRUE),
    row.names = NULL
  ))
  
  # --------------------------------------------------------------------------
  # STEP 9: Normalize by COMBINED column sum
  # --------------------------------------------------------------------------
  step(9, "Normalize by combined column sum")
  
  col_sums <- colSums(m, na.rm = TRUE)
  cat("  Column sums BEFORE normalization:\n")
  print(col_sums)
  
  m_norm <- sweep(m, 2, col_sums, "/")
  
  cat("\n  Column sums AFTER normalization (should be 1):\n")
  print(colSums(m_norm, na.rm = TRUE))
  
  # --------------------------------------------------------------------------
  # STEP 10: Sum intensity per species per replicate
  # --------------------------------------------------------------------------
  step(10, "Sum intensity per species per replicate")
  
  intensity_sp1 <- colSums(m_norm[is_sp1, , drop = FALSE], na.rm = TRUE)
  intensity_sp2 <- colSums(m_norm[is_sp2, , drop = FALSE], na.rm = TRUE)
  
  n_rows_sp1 <- sum(is_sp1)
  n_rows_sp2 <- sum(is_sp2)
  
  cat("  ", species_1, " — rows: ", n_rows_sp1, "\n", sep = "")
  cat("    intensity per replicate:\n")
  print(intensity_sp1)
  
  cat("\n  ", species_2, " — rows: ", n_rows_sp2, "\n", sep = "")
  cat("    intensity per replicate:\n")
  print(intensity_sp2)
  
  # --------------------------------------------------------------------------
  # STEP 11: Build per-replicate table
  # --------------------------------------------------------------------------
  step(11, "Build per-replicate table")
  
  per_replicate <- bind_rows(
    data.frame(species = species_1, replicate = mix_cols,
               n_rows = n_rows_sp1, intensity = as.numeric(intensity_sp1)),
    data.frame(species = species_2, replicate = mix_cols,
               n_rows = n_rows_sp2, intensity = as.numeric(intensity_sp2))
  )
  
  cat("  Per-replicate table:\n")
  print(as.data.frame(per_replicate), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 12: Calculate relative abundance
  # --------------------------------------------------------------------------
  step(12, "Calculate relative abundance per replicate")
  
  rel_per_rep <- per_replicate %>%
    group_by(replicate) %>%
    mutate(
      total_intensity = sum(intensity),
      rel_intensity = intensity / total_intensity
    ) %>%
    ungroup()
  
  cat("  Relative abundance table:\n")
  print(as.data.frame(rel_per_rep), row.names = FALSE)
  
  check_sums <- rel_per_rep %>%
    group_by(replicate) %>%
    summarise(total_rel = sum(rel_intensity), .groups = "drop")
  cat("\n  Sanity check — sum of rel_intensity per replicate (should be 1):\n")
  print(as.data.frame(check_sums), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 13: Summary per species
  # --------------------------------------------------------------------------
  step(13, "Summary per species across replicates")
  
  rel_summary <- rel_per_rep %>%
    group_by(species) %>%
    summarise(
      n_replicates         = n(),
      mean_intensity       = mean(intensity),
      mean_rel_abundance   = mean(rel_intensity),
      median_rel_abundance = median(rel_intensity),
      sd_rel_abundance     = sd(rel_intensity),
      min_rel_abundance    = min(rel_intensity),
      max_rel_abundance    = max(rel_intensity),
      mean_rows            = mean(n_rows),
      .groups = "drop"
    )
  
  cat("  Summary table:\n")
  print(as.data.frame(rel_summary), row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 14: Between-species paired test
  # --------------------------------------------------------------------------
  step(14, "Between-species paired test")
  
  wide_rel <- rel_per_rep %>%
    select(species, replicate, rel_intensity) %>%
    pivot_wider(names_from = species, values_from = rel_intensity)
  
  cat("  Wide format (replicate x species):\n")
  print(as.data.frame(wide_rel), row.names = FALSE)
  
  x1 <- wide_rel[[species_1]]
  x2 <- wide_rel[[species_2]]
  d  <- x1 - x2
  n  <- sum(!is.na(d))
  
  cat("\n  Paired differences (", species_1, " - ", species_2, "):\n", sep = "")
  print(d)
  cat("  n = ", n, ", mean(d) = ", mean(d, na.rm = TRUE),
      ", sd(d) = ", sd(d, na.rm = TRUE), "\n", sep = "")
  
  test_results <- NULL
  if (n >= 2) {
    tt <- tryCatch(t.test(x1, x2, paired = TRUE), error = function(e) NULL)
    
    signs <- expand.grid(rep(list(c(-1, 1)), n))
    perm_means <- apply(signs, 1, function(s) mean(s * d[!is.na(d)]))
    obs <- mean(d, na.rm = TRUE)
    p_two <- mean(abs(perm_means) >= abs(obs) - 1e-12)
    
    dz <- if (sd(d, na.rm = TRUE) > 0) obs / sd(d, na.rm = TRUE) else NA_real_
    se <- sd(d, na.rm = TRUE) / sqrt(n)
    tc <- qt(0.975, df = n - 1)
    
    test_results <- tibble(
      coculture            = mix_condition,
      species_1            = species_1,
      species_2            = species_2,
      n_replicates         = n,
      mean_rel_species_1   = mean(x1, na.rm = TRUE),
      mean_rel_species_2   = mean(x2, na.rm = TRUE),
      mean_difference      = obs,
      ci95_diff_low        = obs - tc * se,
      ci95_diff_high       = obs + tc * se,
      cohens_dz            = dz,
      p_paired_t_two_sided = if (!is.null(tt)) tt$p.value else NA_real_,
      p_perm_exact_two     = p_two,
      direction            = if (obs > 0) paste(species_1, ">", species_2) else paste(species_2, ">", species_1)
    )
    
    cat("\n  Between-species test result:\n")
    print(as.data.frame(test_results), row.names = FALSE)
  } else {
    cat("  Skipping test: n < 2.\n")
  }
  
  return(list(
    per_replicate = rel_per_rep,
    summary       = rel_summary,
    test          = test_results
  ))
}

# ============================================================================
# EXECUTION: Run for all pairs
# ============================================================================

coculture_pairs <- list(
  list(cc="GvFv", s1="Gardnerella vaginalis", f1="GvFv_G_working.csv", s2="Fannyhessea vaginae", f2="GvFv_F_working.csv"),
  list(cc="GvPb", s1="Gardnerella vaginalis", f1="GvPb_G_working.csv", s2="Prevotella bivia",      f2="GvPb_P_working.csv"),
  list(cc="LiFv", s1="Lactobacillus iners",   f1="LiFv_L_working.csv", s2="Fannyhessea vaginae", f2="LiFv_F_working.csv"),
  list(cc="LiGv", s1="Lactobacillus iners",   f1="LiGv_L_working.csv", s2="Gardnerella vaginalis",f2="LiGv_G_working.csv"),
  list(cc="LiPb", s1="Lactobacillus iners",   f1="LiPb_L_working.csv", s2="Prevotella bivia",      f2="LiPb_P_working.csv"),
  list(cc="PbFv", s1="Prevotella bivia",      f1="PbFv_P2_working.csv",s2="Fannyhessea vaginae", f2="PbFv_F2_working.csv")
)

all_summaries <- list()
all_tests     <- list()
all_per_rep   <- list()

for (pair in coculture_pairs) {
  res <- process_combined_coculture(
    file_1 = pair$f1, species_1 = pair$s1,
    file_2 = pair$f2, species_2 = pair$s2,
    mix_condition = pair$cc
  )
  
  all_summaries[[pair$cc]] <- res$summary %>% mutate(coculture = pair$cc, .before = 1)
  all_tests[[pair$cc]]     <- res$test
  all_per_rep[[pair$cc]]   <- res$per_replicate %>% mutate(coculture = pair$cc, .before = 1)
}

# Combine
final_summary <- bind_rows(all_summaries)
final_tests   <- bind_rows(all_tests)

# Print to console
cat("FINAL SUMMARY ACROSS ALL COCULTURES")
print(as.data.frame(final_summary), row.names = FALSE)

cat("BETWEEN-SPECIES TESTS")
print(as.data.frame(final_tests), row.names = FALSE)

# Save to disk
write.csv(final_summary, "ALL_split_cocultures_summary.csv", row.names = FALSE)
write.csv(final_tests,   "ALL_split_cocultures_tests.csv", row.names = FALSE)
write.csv(bind_rows(all_per_rep), "ALL_split_cocultures_per_replicate.csv", row.names = FALSE)

cat("\nFiles saved to disk:\n")
cat("  - ALL_split_cocultures_summary.csv\n")
cat("  - ALL_split_cocultures_tests.csv\n")
cat("  - ALL_split_cocultures_per_replicate.csv\n")
