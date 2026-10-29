#Heatmap for untergeted metabolites

rm(list = ls())
setwd("your/working/path")
getwd()

library(pheatmap)
library(RColorBrewer)
library(svglite)

# --- 1. Define file paths ---
input_file     <- "VaBaMet_matrix.csv"
output_png     <- "VaBaMet_heatmap.png"
output_svg     <- "VaBaMet_heatmap.svg"

# --- 2. Read the CSV file ---
df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

intensity_cols <- grep("^\\d\\.[A-Za-z]+$", colnames(df), value = TRUE)

# --- 3. Create row names: use 'short name' if available, else 'name' ---
if ("short name" %in% colnames(df)) {
  row_names <- df$`short name`
  row_names[is.na(row_names)] <- ""
} else {
  row_names <- df$name
}

mat <- as.matrix(df[, intensity_cols])
rownames(mat) <- row_names

# --- 4. Data transformation --------------------------------------------
# 4a. Row‑wise zero replacement: replace zeros with half the row's minimum positive value,
#     then log2‑transform (no additional offset needed)
mat_nozero <- mat   # copy
for (i in seq_len(nrow(mat))) {
  x <- mat[i, ]
  pos <- x[x > 0]
  if (length(pos) > 0) {
    half_min <- min(pos) / 2
    x[x == 0] <- half_min
  } else {
    # No positive values in this row – all zeros; leave as is
    # (the row will likely be removed later due to zero variance)
  }
  mat_nozero[i, ] <- x
}
mat_log <- log2(mat_nozero)

# 4b. Pareto scaling (row‑wise): mean‑centre, divide by sqrt(SD)
pareto_scale <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  if (sigma == 0) return(x - mu)   # if SD = 0 → centre only
  (x - mu) / sqrt(sigma)
}
mat_scaled <- t(apply(mat_log, 1, pareto_scale))

# 4c. Remove rows that become constant (e.g., all zeros after log)
constant_rows <- apply(mat_scaled, 1, function(x) all(is.na(x) | sd(x, na.rm = TRUE) == 0))
if (any(constant_rows)) {
  warning("Removing ", sum(constant_rows), " rows with zero variance.")
  mat_scaled <- mat_scaled[!constant_rows, ]
}

# --- 5. Prepare column annotations with custom colors -----------------
sample_groups <- sub("^\\d\\.(.*)$", "\\1", colnames(mat_scaled))
annotation_col <- data.frame(Group = sample_groups)
rownames(annotation_col) <- colnames(mat_scaled)

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
annotation_colors <- list(Group = group_colors)

# --- 6. Symmetric colour palette centred on zero (white) ---------------
max_abs <- max(abs(mat_scaled), na.rm = TRUE)
n_colors <- 99
heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(n_colors)
breaks <- seq(-max_abs, max_abs, length.out = n_colors + 1)

# --- 7. Create the heatmap --------------------------------------------
hm <- pheatmap(mat_scaled,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "single",      # average linkage
               color = heatmap_colors,
               breaks = breaks,                    # symmetric around zero
               scale = "none",
               border_color = NA,
               fontsize_row = 6,
               fontsize_col = 8,
               annotation_col = annotation_col,
               annotation_colors = annotation_colors,
               annotation_legend = TRUE,
               na_col = "grey80",
               silent = TRUE
)

# --- 9. Save PNG -------------------------------------------------------
png(filename = output_png, width = 7, height = 9, units = "in", res = 300)
print(hm)
dev.off()
cat("PNG saved to:", output_png, "\n")

# --- 10. Save SVG with fallback ----------------------------------------
save_svg <- function(hm, file, width = 7, height = 9) {
  cairo_ok <- capabilities("cairo") && requireNamespace("Cairo", quietly = TRUE)
  if (cairo_ok) {
    svg(filename = file, width = width, height = height)
    print(hm)
    dev.off()
    cat("SVG saved (using Cairo) to:", file, "\n")
  } else {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Cairo not available and 'svglite' not installed. Install svglite: install.packages('svglite')")
    }
    svglite::svglite(file = file, width = width, height = height)
    print(hm)
    dev.off()
    cat("SVG saved (using svglite) to:", file, "\n")
  }
}

save_svg(hm, output_svg)
