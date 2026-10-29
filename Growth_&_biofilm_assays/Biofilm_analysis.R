#Biofilm analysis

rm(list = ls())
setwd("your/working/path")
getwd()

# =============================================================================
# Biofilm formation analysis – OD values (boxplots + SVG output)
# =============================================================================

library(tidyr)      # pivot_longer, separate
library(dplyr)      # data manipulation
library(ggplot2)    # plotting
library(ggpubr)     # stat_compare_means (optional, for easy annotation)
library(broom)      # tidy test results
library(svglite)    # SVG graphics device

# -----------------------------------------------------------------------------
# 1. Read the CSV file (adjust the file path as needed)
# -----------------------------------------------------------------------------
df <- read.csv("joined_biofilm_ODs.csv", check.names = FALSE, stringsAsFactors = FALSE)

# Add a replicate ID (row number) to keep track of individual measurements
df$Replicate <- seq_len(nrow(df))

# -----------------------------------------------------------------------------
# 2. Convert from wide to long format
# -----------------------------------------------------------------------------
df_long <- df %>%
  pivot_longer(
    cols = -Replicate,
    names_to = "Condition",
    values_to = "OD"
  )

# -----------------------------------------------------------------------------
# 3. Parse the condition column into organism, treatment, concentration and type
# -----------------------------------------------------------------------------
parse_condition <- function(cond) {
  parts <- strsplit(cond, "_")[[1]]
  org <- parts[1]                     # Gv, Gp, Pb, Lc, Ca, Ec or Sa
  
  # Control conditions: DM (for H) and DBM (for S)
  if (length(parts) == 2) {
    control_label <- parts[2]
    treatment <- ifelse(control_label == "DM", "H", "S")  # DM is H control, DBM is S control
    conc <- "Control"
    type <- control_label
  } else {
    # Treatment conditions: e.g. Ca_S_1.0
    treatment <- parts[2]              # S or H
    conc <- parts[3]                   # 0.01, 0.1, 1.0
    type <- paste(treatment, conc, sep = "_")
  }
  
  data.frame(Organism = org,
             Treatment = treatment,
             Concentration = conc,
             Type = type,
             stringsAsFactors = FALSE)
}

condition_info <- do.call(rbind, lapply(df_long$Condition, parse_condition))
df_long <- cbind(df_long, condition_info)

# Ensure correct ordering of concentrations for the x-axis
df_long$Concentration <- factor(df_long$Concentration,
                                levels = c("Control", "0.01", "0.1", "1.0"))

# The "compound" for colouring is the Treatment (S or H)
df_long$Compound <- factor(df_long$Treatment, levels = c("S", "H"))

# -----------------------------------------------------------------------------
# 4. Boxplots – coloured by compound (S/H), faceted by organism
# -----------------------------------------------------------------------------
p <- ggplot(df_long, aes(x = Concentration, y = OD, fill = Compound)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.8, outlier.shape = NA) +
  # Add individual points (jittered) for transparency
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1.5, alpha = 0.6, pch = 21, colour = "grey30") +
  scale_fill_manual(values = c("S" = "#56BC37", "H" = "#2C6E96"),
                    labels = c("S (control DBM)", "H (control DM)")) +
  labs(x = "Concentration", y = "OD (biofilm formation)",
       fill = "Compound") +
  facet_wrap(~ Organism, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom")

# Display the plot in the R session
print(p)

# Save as SVG using svglite
svglite("biofilm_boxplot.svg", width = 5, height = 8)
print(p)
dev.off()

# -----------------------------------------------------------------------------
# 5. Statistics – Wilcoxon rank-sum test with effect sizes and descriptives
# -----------------------------------------------------------------------------
library(effsize)   # for Cliff's delta (install.packages("effsize"))

treatment_data <- df_long %>% filter(Concentration != "Control")
control_data   <- df_long %>% filter(Concentration == "Control")

groups <- treatment_data %>%
  distinct(Organism, Treatment, Concentration)

stat_results <- list()

for (i in seq_len(nrow(groups))) {
  org  <- groups$Organism[i]
  trt  <- groups$Treatment[i]
  conc <- groups$Concentration[i]
  
  # --- Pull the two groups -------------------------------------------------
  trt_od <- treatment_data %>%
    filter(Organism == org, Treatment == trt, Concentration == conc) %>%
    pull(OD)
  
  ctrl_label <- ifelse(trt == "S", "DBM", "DM")
  
  ctrl_od <- control_data %>%
    filter(Organism == org, Type == ctrl_label) %>%
    pull(OD)
  
  # --- Only test if both groups have >= 2 observations ---------------------
  if (length(trt_od) >= 2 && length(ctrl_od) >= 2) {
    
    # Two-tailed Wilcoxon rank-sum (Mann–Whitney U)
    test_res <- wilcox.test(trt_od, ctrl_od, exact = FALSE,
                            conf.int = TRUE, conf.level = 0.95)
    
    n1 <- length(trt_od)
    n2 <- length(ctrl_od)
    U  <- as.numeric(test_res$statistic)
    
    # Rank-biserial correlation (non-parametric effect size, −1 to +1)
    r_rb <- 1 - (2 * U) / (n1 * n2)
    
    # Cliff's delta (equivalent interpretation, robust to ties)
    cliff_d <- tryCatch({
      cd <- effsize::cliff.delta(trt_od, ctrl_od)
      as.numeric(cd$estimate)
    }, error = function(e) NA_real_)
    
    # Hodges–Lehmann estimate of median difference + 95% CI
    hl_diff <- tryCatch(as.numeric(test_res$estimate), error = function(e) NA_real_)
    hl_ci   <- tryCatch(as.numeric(test_res$conf.int),   error = function(e) c(NA, NA))
    
    # --- Assemble the row --------------------------------------------------
    stat_results[[i]] <- data.frame(
      Organism          = org,
      Treatment         = trt,
      Concentration     = conc,
      Control           = ctrl_label,
      N_Treatment       = n1,
      N_Control         = n2,
      Median_Treatment  = median(trt_od, na.rm = TRUE),
      IQR_Treatment     = IQR(trt_od,    na.rm = TRUE),
      Median_Control    = median(ctrl_od, na.rm = TRUE),
      IQR_Control       = IQR(ctrl_od,    na.rm = TRUE),
      Delta_Median      = median(trt_od, na.rm = TRUE) - median(ctrl_od, na.rm = TRUE),
      HL_Diff           = hl_diff,
      HL_CI_Low         = hl_ci[1],
      HL_CI_High        = hl_ci[2],
      U                 = U,
      r_rb              = r_rb,
      Cliffs_delta      = cliff_d,
      p.value           = test_res$p.value,
      stringsAsFactors  = FALSE
    )
  }
}

stat_table <- bind_rows(stat_results)

# --- Significance stars (based on raw p-value) -------------------------------
stat_table$signif <- cut(stat_table$p.value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***", "**", "*", "ns"),
                         right  = FALSE)

# --- Round numeric columns for readability ----------------------------------
stat_table <- stat_table %>%
  mutate(across(
    c(Median_Treatment, IQR_Treatment, Median_Control, IQR_Control,
      Delta_Median, HL_Diff, HL_CI_Low, HL_CI_High,
      U, r_rb, Cliffs_delta, p.value),
    ~ round(.x, 4)
  ))

# --- Print a compact results table ------------------------------------------
cat("\n=== BIOFILM STATISTICAL RESULTS (Wilcoxon rank-sum, uncorrected) ===\n\n")
print(stat_table %>%
        select(Organism, Treatment, Concentration, Control,
               N_Treatment, N_Control,
               Median_Treatment, IQR_Treatment,
               Median_Control,   IQR_Control,
               Delta_Median, HL_Diff, HL_CI_Low, HL_CI_High,
               U, r_rb, Cliffs_delta, p.value, signif))

# --- Save full table --------------------------------------------------------
write.csv(stat_table, "biofilm_statistics.csv", row.names = FALSE)
cat("\nFull statistics table saved to 'biofilm_statistics.csv'.\n")

# --- Also produce a condensed table -----------------------
pub_table <- stat_table %>%
  mutate(
    `Median (IQR) Treatment` = sprintf("%.2f (%.2f)", Median_Treatment, IQR_Treatment),
    `Median (IQR) Control`   = sprintf("%.2f (%.2f)", Median_Control,   IQR_Control),
    `Δ Median [95% CI]`      = sprintf("%.2f [%.2f, %.2f]",
                                       Delta_Median, HL_CI_Low, HL_CI_High),
    `U`                      = round(U, 1),
    `r_rb`                   = round(r_rb, 2),
    `p`                      = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  select(Organism, Treatment, Concentration, Control,
         `Median (IQR) Treatment`, `Median (IQR) Control`,
         `Δ Median [95% CI]`, U, r_rb, p, signif)

write.csv(pub_table, "biofilm_statistics_condensed.csv", row.names = FALSE)
cat("Condensed table saved to 'biofilm_statistics_condensed.csv'.\n")
