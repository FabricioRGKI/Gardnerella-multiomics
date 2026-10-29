# ----------------------------------------------------------------------
# Script: Targeted boxplots
#          Input file should contain technical replicate columns (##.XX.##)
# For single metabolites
# Output: PNG and SVG files.
# ----------------------------------------------------------------------

rm(list = ls())
setwd("your/working/folder")
getwd()

library(tidyverse)
library(svglite)

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
# 3. Define colour mapping (as originally provided)
# ----------------------------------------------------------------------
condition_colors <- c(
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

# ----------------------------------------------------------------------
# 4. Identify technical replicate columns (pattern: number.letters.number)
# ----------------------------------------------------------------------
tech_cols <- grep("^\\d+\\.[A-Za-z]+\\.\\d+$", names(data), value = TRUE)
if (length(tech_cols) == 0) {
  stop("No columns matching pattern '##.XX.##' (technical replicates) found.")
}

# ----------------------------------------------------------------------
# 5. Reshape all data to long format (keep 'name' column for filtering)
# ----------------------------------------------------------------------
long_data <- data %>%
  select(name, all_of(tech_cols)) %>%   # removed 'short name'
  pivot_longer(cols = all_of(tech_cols), names_to = "colname", values_to = "intensity") %>%
  separate(colname, into = c("replicate", "condition", "tech_rep"), sep = "\\.", remove = TRUE) %>%
  select(-tech_rep) %>%
  mutate(replicate = as.integer(replicate))

# ----------------------------------------------------------------------
# 6. Filter for the metabolite of interest (using 'name' column)
# ----------------------------------------------------------------------
#For (e.g.) Methyl-beta-carboline
target_metabolite <- "Methyl-beta-carboline"   # exact match as in 'name' column

metab_data <- long_data %>%
  filter(name == target_metabolite)

if (nrow(metab_data) == 0) {
  stop("Metabolite '", target_metabolite, "' not found in 'name' column.")
}

# ----------------------------------------------------------------------
# 7. Keep only specified conditions
# ----------------------------------------------------------------------
#For Pb cocultures comparisons
keep_conditions <- c("PbFv", "LiPb", "GvPb", "LcPb", "Lc", "Li", "Gv", "Pb", "Fv")

#For Gv cocultures comparisons
keep_conditions <- c("LcGv", "LiGv", "GvPb", "GvFv", "Lc", "Li", "Gv", "Pb", "Fv")

metab_data <- metab_data %>%
  filter(condition %in% keep_conditions)

# ----------------------------------------------------------------------
# 8. Order conditions on x‑axis (custom order, only those present)
# ----------------------------------------------------------------------
#For Pb cocultures comparisons
desired_order <- c("Lc", "LcLi", "LcGv", "LcPb", "LcFv", 
                   "Li", "LiGv", "LiPb", "LiFv", 
                   "Gv", "GvPb", "GvFv", 
                   "Pb", "PbFv", "Fv")
#For Gv cocultures comparisons
desired_order <- c("Lc", "LcLi", "LcGv", "LcPb", "LcFv", 
                   "Li", "LiGv", "LiPb", "LiFv", 
                   "Gv", "GvPb", "Pb", "GvFv", "Fv")

# Determine final order based on which desired conditions are actually present
final_order <- intersect(desired_order, unique(metab_data$condition))
extra <- setdiff(unique(metab_data$condition), desired_order)
final_order <- c(final_order, sort(extra))

# Convert condition to a factor with the specified levels
metab_data$condition <- factor(metab_data$condition, levels = final_order)

# ----------------------------------------------------------------------
# 9. Create the plot with requested styling
# ----------------------------------------------------------------------
p <- ggplot(metab_data, aes(x = condition, y = intensity)) +
  
  # ---------- Black jittered points (all data visible) ----------
# alpha: transparency (0.5 = 50% opaque)
# width: horizontal jitter spread
# size: point diameter
geom_jitter(colour = "black", alpha = 0.5, width = 0.2, size = 2.2) +
  
  # ---------- Boxplot with thin black outline ----------
# colour: outline colour and linewidth: outline thickness
# fill: uses condition_colors mapping
# outlier.shape = NA because points already shown as jitter
geom_boxplot(aes(fill = condition),
             colour = "black", linewidth = 0.3,
             outlier.shape = NA, alpha = 0.9) +
  
# ---------- Median line (thinner than default) ----------
# middle.linewidth controls the thickness of the median bar
stat_summary(fun = median, geom = "crossbar", width = 0.5, colour = "black",
             middle.linewidth = 0.4, show.legend = FALSE) +
  
  # ---------- Fill colours from the original palette ----------
scale_fill_manual(values = condition_colors, na.value = "gray70") +
  
  # ---------- Labels ----------
labs(title = target_metabolite,
     x = "Condition", y = "Intensity") +
  
  # ---------- Minimal journal style ----------
theme_classic(base_size = 14) +
  theme(
    # Axis text: angled to avoid overlap
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = rel(1.65)),
    axis.text.y = element_text(colour = "black"),
    # Axis lines and ticks
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black"),
    # Remove all gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Hide legend (only one metabolite)
    legend.position = "none",
    # Centre and bold title
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ----------------------------------------------------------------------
# 10. Save outputs (PNG + SVG)
# ----------------------------------------------------------------------
#methyl_beta_carboline
ggsave("MbC_single_metabolites.png", plot = p, width = 8, height = 6, dpi = 300)
svglite("MbC_single_metabolites.svg", width = 8, height = 6)
print(p)
dev.off()
