# =============================================================================
# R script: PCoA with Overlaid Ellipses for Lc, Gv, and LcGv
# =============================================================================
# 
# This script:
#   1. Reads the metabolite dataset.
#   2. Preprocesses (log10, Z‑score).
#   3. Performs PCoA (Manhattan distance).
#   4. Creates scatter plots with:
#        - Points colored by group, shaped by biological replicate.
#        - Ellipses around: all Lc‑containing groups, all Gv‑containing groups,
#          and the specific LcGv group.
#   5. Saves as PNG and SVG.
#
# =============================================================================

library(ggplot2)
library(svglite)

# --- 1. Define file paths ---
input_file   <- "VaBaMet_matrix.csv"   # adjust if needed
output_dir   <- "VaBaMet_PCoA"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 2. Read and prepare data ---
df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

# Changed pattern: digit, dot, letters (e.g., 1.Fv, 2.LcPb, 3.GvFv)
intensity_cols <- grep("^\\d\\.[A-Za-z]+$", colnames(df), value = TRUE)
if (length(intensity_cols) == 0) stop("No 'digit.letters' columns found (e.g., 1.Fv).")

# Parse sample information
sample_names <- intensity_cols
bio   <- sub("^(\\d+)\\..*$", "\\1", sample_names)           # digit before dot
group <- sub("^\\d+\\.(.*)$", "\\1", sample_names)          # letters after dot

sample_info <- data.frame(
  sample = sample_names,
  bio    = bio,
  group  = group,
  stringsAsFactors = FALSE
)

# Add grouping categories for ellipses (same logic)
sample_info$ellipse_group <- "Other"
sample_info$ellipse_group[grepl("Lc", sample_info$group) & sample_info$group != "LcGv"] <- "Lc-containing"
sample_info$ellipse_group[grepl("Gv", sample_info$group) & sample_info$group != "LcGv"] <- "Gv-containing"
sample_info$ellipse_group[sample_info$group == "LcGv"] <- "LcGv"
print(unique(sample_info[, c("group", "ellipse_group")]))

# Create numeric matrix: samples × metabolites
mat <- as.matrix(df[, intensity_cols])
rownames(mat) <- df$`short name`          # metabolites as row names
mat <- t(mat)                             # now rows = samples, columns = metabolites

# --- 3. Data transformation --------------------------------------------
# 3a. For each metabolite (column), replace zeros with half the column's min positive value
mat_nozero <- mat
for (j in seq_len(ncol(mat))) {
  col_vals <- mat[, j]
  pos <- col_vals[col_vals > 0]
  if (length(pos) > 0) {
    half_min <- min(pos) / 2
    col_vals[col_vals == 0] <- half_min
  }
  # if no positives, leave as zeros (will be removed later if constant)
  mat_nozero[, j] <- col_vals
}

# 3b. log2 transformation
mat_log <- log2(mat_nozero)

# 3c. Pareto scaling per metabolite (column‑wise)
pareto_scale_col <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  if (sigma == 0) return(x - mu)
  (x - mu) / sqrt(sigma)
}
mat_scaled <- apply(mat_log, 2, pareto_scale_col)

# Remove constant columns (zero variance after scaling)
constant_cols <- apply(mat_scaled, 2, function(x) sd(x, na.rm = TRUE) == 0)
if (any(constant_cols)) {
  mat_scaled <- mat_scaled[, !constant_cols]
  cat("Removed", sum(constant_cols), "constant metabolites.\n")
}

# --- 4. PCoA (Manhattan distance) ---
dist_mat <- dist(mat_scaled, method = "manhattan")
pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

pcoa_scores <- as.data.frame(pcoa_res$points)
colnames(pcoa_scores) <- c("PCo1", "PCo2")
pcoa_scores <- cbind(sample_info, pcoa_scores)

pcoa_var <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)[1:2]

# --- 5. Custom color palette ---
group_colors <- c(
  "Lc"   = "#E3377E",
  "LcLi" = "#C9319C",
  "LcGv" = "#722CB9",
  "LcPb" = "#9C7A5A",
  "LcFv" = "#936796",
  "Li"   = "#B02CBA",
  "LiGv" = "#5927D8",
  "LiPb" = "#837478",
  "LiFv" = "#7A62B4",
  "Gv"   = "#0121F5",
  "GvPb" = "#2C6E96",
  "GvFv" = "#225CD2",
  "Pb"   = "#56BC37",
  "PbFv" = "#4DA973",
  "Fv"   = "#4497AF"
)

# --- 6. Plotting function ---
create_ordination_plot <- function(scores, x_col, y_col, x_var, y_var, title) {
  ellipse_data <- scores[scores$ellipse_group != "Other", ]
  
  p <- ggplot(scores, aes_string(x = x_col, y = y_col)) +
    stat_ellipse(data = ellipse_data,
                 aes(color = ellipse_group, fill = ellipse_group),
                 geom = "polygon", alpha = 0.15, level = 0.9, size = 0.2) +
    geom_point(aes(color = group, shape = bio), size = 3, alpha = 0.85) +
    scale_color_manual(
      name = "Group / Ellipse",
      values = group_colors,
      breaks = names(group_colors),
      guide = guide_legend(override.aes = list(shape = 16))
    ) +
    scale_fill_manual(
      name = "Ellipse (95% coverage)",
      values = c("Lc-containing" = "#f19ec1",
                 "Gv-containing" = "#91a0ff",
                 "LcGv" = "#C19FDF"),
      guide = guide_legend(override.aes = list(alpha = 0.3))
    ) +
    scale_shape_manual(
      name = "Biological Replicate",
      values = c("1" = 16, "2" = 17, "3" = 15)
    ) +
    labs(
      title = title,
      x = paste0(x_col, " (", x_var, "% variance)"),
      y = paste0(y_col, " (", y_var, "% variance)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.box = "vertical"
    ) +
    coord_fixed(ratio = 1)
  return(p)
}

# --- 8. Generate plot ---
pcoa_plot <- create_ordination_plot(
  pcoa_scores, "PCo1", "PCo2", pcoa_var[1], pcoa_var[2],
  "PCoA of Bacterial Metabolite Profiles (Manhattan distance)\nwith Lc / Gv / LcGv Ellipses"
)

# --- 9. Save plot ---
ggsave(file.path(output_dir, "PCoA_manhattan.png"), pcoa_plot, width = 6, height = 6, dpi = 300)
ggsave(file.path(output_dir, "PCoA_manhattan.svg"), pcoa_plot, width = 6, height = 6, device = svglite)

cat("Plots saved to:", output_dir, "\n")
                       
