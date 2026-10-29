#Growth curve and inhibition test
rm(list = ls())
setwd("your/working/path")
getwd()

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(svglite)
library(broom)

# -----------------------------------------------------------------------------
# USER CONFIGURATION: Metabolite to control mapping
# -----------------------------------------------------------------------------
# For each metabolite (S = Sphingosine, H = Harmane), specify the control code
metabolite_to_control <- c(
  "S" = "DBM",   # DBM is control for Sphingosine
  "H" = "DM"     # DM is control for Harmane
)

# -----------------------------------------------------------------------------
# Read and combine the three experimental replicates
# -----------------------------------------------------------------------------
cat("\n=== READING EXPERIMENTAL REPLICATES ===\n")

# Replace with your actual file names
data <- read.csv("xp01.csv", sep = ",", check.names = FALSE) %>%
  mutate(Experiment = "Exp1")

data <- read.csv("xp02.csv", sep = ",", check.names = FALSE) %>%
  mutate(Experiment = "Exp2")

data <- read.csv("xp03.csv", sep = ",", check.names = FALSE) %>%
  mutate(Experiment = "Exp3")

# If you have multiple files, combine them with bind_rows()

# Combine all data
data <- bind_rows(data1, data2, data3)

cat("Data loaded successfully!\n")
cat("Total rows:", nrow(data), "\n")
cat("Columns:", ncol(data), "\n")

# -----------------------------------------------------------------------------
# Helper: parse column name into components (new format)
# -----------------------------------------------------------------------------
parse_column_name <- function(colname) {
  # Expected formats:
  # Treatment:   Organism_MetaboliteReplicate_Concentration  e.g. "Ca_S1_1"
  # Control:     Organism_ControlCodeReplicate               e.g. "Ca_DM1" or "Ca_DBM2"
  # Media wells: "M_S_1", "M_H_0.1", "M_DM", "M_DBM" -> skip
  # Also skip well names like "A1", "B12", "H1", etc.
  
  # Skip well names (single letter + number)
  if (grepl("^[A-H][0-9]+$", colname)) return(NULL)
  
  # Skip columns starting with "M_" (media-only controls)
  if (grepl("^M_", colname)) return(NULL)
  
  # Pattern for treatments: Organism_S[1-3]_0.1 or 1  (or H[1-3])
  treatment_pattern <- "^([A-Za-z]+)_([SH])([1-3])_([0-9.]+)$"
  if (grepl(treatment_pattern, colname)) {
    parts <- str_match(colname, treatment_pattern)
    organism <- parts[2]
    metabolite <- parts[3]      # "S" or "H"
    replicate <- as.integer(parts[4])
    concentration <- as.numeric(parts[5])
    return(list(
      type = "treatment",
      organism = organism,
      metabolite = metabolite,
      replicate = replicate,
      concentration = concentration,
      full_code = colname
    ))
  }
  
  # Pattern for controls: Organism_DM[1-3] or DBM[1-3]
  control_pattern <- "^([A-Za-z]+)_(DM|DBM)([1-3])$"
  if (grepl(control_pattern, colname)) {
    parts <- str_match(colname, control_pattern)
    organism <- parts[2]
    control_code <- parts[3]    # "DM" or "DBM"
    replicate <- as.integer(parts[4])
    # Determine which metabolite this control corresponds to
    metabolite <- ifelse(control_code == "DM", "H", "S")
    return(list(
      type = "control",
      organism = organism,
      metabolite = metabolite,
      control_code = control_code,
      replicate = replicate,
      full_code = colname
    ))
  }
  
  # Unknown format: skip
  return(NULL)
}

# -----------------------------------------------------------------------------
# Assign compound/metabolite info
# -----------------------------------------------------------------------------
get_compound_info <- function(parsed) {
  if (parsed$type == "treatment") {
    compound_name <- paste0(parsed$metabolite, "_", parsed$concentration, "mM")
    return(list(
      type = "Treatment",
      compound = compound_name,
      metabolite = parsed$metabolite,
      concentration_numeric = parsed$concentration,
      concentration_label = paste0(parsed$concentration, " mM"),
      control_code = NA_character_
    ))
  } else if (parsed$type == "control") {
    return(list(
      type = "Control",
      compound = parsed$control_code,   # "DM" or "DBM"
      metabolite = parsed$metabolite,
      concentration_numeric = NA_real_,
      concentration_label = "Control",
      control_code = parsed$control_code
    ))
  } else {
    return(NULL)
  }
}

# -----------------------------------------------------------------------------
# Main processing function (adapted for new format)
# -----------------------------------------------------------------------------
create_growth_curves <- function(data) {
  cat("Data dimensions:", dim(data), "\n")
  cat("Processing", length(unique(data$Experiment)), "experiments together...\n")
  
  # Extract time from first column (has "s" suffix)
  time_col <- data[[1]]
  time_numeric <- as.numeric(str_remove(time_col, "s$"))
  data$Time <- time_numeric / 3600   # convert to hours
  
  # Remove the raw time label and temperature columns (first two columns)
  data <- data[, -c(1, 2)]
  
  # Identify columns to keep: those that are NOT just well names (A1, H12, etc.)
  well_pattern <- "^[A-H][0-9]+$"
  cols_to_keep <- c("Time", "Experiment")
  for (col in names(data)) {
    if (!grepl(well_pattern, col) && col != "Time" && col != "Experiment") {
      cols_to_keep <- c(cols_to_keep, col)
    }
  }
  data <- data[, cols_to_keep]
  
  processed_data <- list()
  
  # Loop over all measurement columns
  for (col_name in names(data)[!names(data) %in% c("Time", "Experiment")]) {
    parsed <- parse_column_name(col_name)
    if (is.null(parsed)) next
    
    info <- get_compound_info(parsed)
    if (is.null(info)) next
    
    # For each experiment separately
    for (exp in unique(data$Experiment)) {
      exp_data <- data[data$Experiment == exp, ]
      
      col_data <- data.frame(
        Time = exp_data$Time,
        Well = paste0(col_name, "_", exp),
        OD = exp_data[[col_name]],
        Organism = parsed$organism,
        Compound = info$compound,
        Metabolite = info$metabolite,
        Concentration_numeric = info$concentration_numeric,
        Concentration_label = info$concentration_label,
        Replicate = parsed$replicate,
        Experiment = exp,
        Type = info$type,
        Control_Code = info$control_code,
        stringsAsFactors = FALSE
      )
      
      key <- paste(parsed$organism, info$compound, exp, sep = "_")
      if (!key %in% names(processed_data)) {
        processed_data[[key]] <- col_data
      } else {
        processed_data[[key]] <- rbind(processed_data[[key]], col_data)
      }
    }
  }
  
  combined_data <- bind_rows(processed_data)
  
  cat("\nProcessed data summary:\n")
  cat("Total wells:", length(unique(combined_data$Well)), "\n")
  cat("Experiments:", paste(unique(combined_data$Experiment), collapse=", "), "\n")
  cat("Organisms:", paste(unique(combined_data$Organism), collapse=", "), "\n")
  cat("Metabolites:", paste(unique(combined_data$Metabolite), collapse=", "), "\n")
  cat("Compounds:", paste(unique(combined_data$Compound), collapse=", "), "\n")
  
  # Calculate mean and SD per time point, grouping by relevant factors
  summary_data <- combined_data %>%
    group_by(Time, Organism, Metabolite, Compound, Concentration_label, 
             Type, Control_Code, Experiment) %>%
    summarise(
      Mean_OD = mean(OD, na.rm = TRUE),
      SD_OD = sd(OD, na.rm = TRUE),
      SEM_OD = SD_OD / sqrt(n()),
      N = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      Upper = Mean_OD + SD_OD,
      Lower = pmax(Mean_OD - SD_OD, 0)
    )
  
  return(list(raw_data = combined_data, summary_data = summary_data))
}

# -----------------------------------------------------------------------------
# Plotting functions (adjusted for new groupings)
# -----------------------------------------------------------------------------
plot_growth_curves_faceted <- function(summary_data) {
  # Separate treatments and controls
  treatment_data <- summary_data %>% filter(Type == "Treatment")
  control_data   <- summary_data %>% filter(Type == "Control")
  
  # Prepare expanded control data: for each organism and treatment compound,
  # attach the corresponding control (based on metabolite)
  expanded_controls <- list()
  
  for (org in unique(treatment_data$Organism)) {
    org_treatments <- treatment_data %>% filter(Organism == org)
    for (met in unique(org_treatments$Metabolite)) {
      # Find the control code for this metabolite
      ctrl_code <- metabolite_to_control[met]
      if (is.na(ctrl_code)) next
      
      ctrl_sub <- control_data %>%
        filter(Organism == org, Control_Code == ctrl_code)
      if (nrow(ctrl_sub) == 0) next
      
      # Assign the control to all treatment compounds of this metabolite
      # (i.e., both concentrations will share the same control)
      treat_compounds <- org_treatments %>% filter(Metabolite == met) %>% pull(Compound) %>% unique()
      for (comp in treat_compounds) {
        ctrl_copy <- ctrl_sub %>% mutate(Compound = comp)
        expanded_controls[[paste(org, met, comp)]] <- ctrl_copy
      }
    }
  }
  
  expanded_controls_df <- bind_rows(expanded_controls)
  
  # Create a colour palette for compounds (metabolite + concentration)
  compounds <- unique(treatment_data$Compound)
  # Define colors: Sphingosine in greens, Harmane in teals
  n_S <- sum(grepl("^S_", compounds))
  n_H <- sum(grepl("^H_", compounds))
  color_S <- colorRampPalette(c("#56BC37", "#4D9537"))(n_S)
  color_H <- colorRampPalette(c("#2C6E96", "#2C5B77"))(n_H)
  compound_colors <- setNames(c(color_S, color_H), compounds)
  
  ggplot() +
    # Treatment ribbons and lines
    geom_ribbon(data = treatment_data,
                aes(x = Time, ymin = Lower, ymax = Upper,
                    group = interaction(Compound, Experiment),
                    fill = Compound),
                alpha = 0.2) +
    geom_line(data = treatment_data,
              aes(x = Time, y = Mean_OD,
                  group = interaction(Compound, Experiment),
                  color = Compound),
              linewidth = 0.3) +
    
    # Control ribbons and lines (black dashed)
    geom_ribbon(data = expanded_controls_df,
                aes(x = Time, ymin = Lower, ymax = Upper,
                    group = interaction(Compound, Experiment)),
                alpha = 0.15, fill = "black") +
    geom_line(data = expanded_controls_df,
              aes(x = Time, y = Mean_OD,
                  group = interaction(Compound, Experiment)),
              linewidth = 0.3, color = "black", linetype = "dashed") +
    
    facet_grid(Organism ~ Compound, scales = "free_y") +
    labs(
      title = "Growth Curves: Treatment vs. Respective Control",
      subtitle = "Solid coloured lines = treatments (Sphingosine or Harmane); black dashed lines = matched controls (DBM for S, DM for H)\nShaded areas = ± SD within each experiment",
      x = "Time (hours)", y = "Optical Density (OD)",
      color = "Treatment", fill = "Treatment"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0, hjust = 1)
    ) +
    scale_color_manual(values = compound_colors) +
    scale_fill_manual(values = compound_colors)
}

# Alternative plot (by organism, not grid by compound)
plot_growth_curves_by_organism <- function(summary_data) {
  treatment_data <- summary_data %>% filter(Type == "Treatment")
  control_data   <- summary_data %>% filter(Type == "Control")
  
  # Expand controls similarly
  expanded_controls <- list()
  for (org in unique(treatment_data$Organism)) {
    org_treatments <- treatment_data %>% filter(Organism == org)
    for (met in unique(org_treatments$Metabolite)) {
      ctrl_code <- metabolite_to_control[met]
      if (is.na(ctrl_code)) next
      ctrl_sub <- control_data %>% filter(Organism == org, Control_Code == ctrl_code)
      if (nrow(ctrl_sub) == 0) next
      expanded_controls[[paste(org, met)]] <- ctrl_sub %>% mutate(Metabolite = met)
    }
  }
  expanded_controls_df <- bind_rows(expanded_controls)
  
  compounds <- unique(treatment_data$Compound)
  n_S <- sum(grepl("^S_", compounds))
  n_H <- sum(grepl("^H_", compounds))
  color_S <- colorRampPalette(c("#56BC37", "#4D9537"))(n_S)
  color_H <- colorRampPalette(c("#2C6E96", "#2C5B77"))(n_H)
  compound_colors <- setNames(c(color_S, color_H), compounds)
  
  ggplot() +
    geom_ribbon(data = treatment_data,
                aes(x = Time, ymin = Lower, ymax = Upper,
                    group = interaction(Compound, Experiment),
                    fill = Compound), alpha = 0.2) +
    geom_line(data = treatment_data,
              aes(x = Time, y = Mean_OD,
                  group = interaction(Compound, Experiment),
                  color = Compound), linewidth = 0.3) +
    geom_ribbon(data = expanded_controls_df,
                aes(x = Time, ymin = Lower, ymax = Upper,
                    group = interaction(Metabolite, Experiment)),
                alpha = 0.15, fill = "black") +
    geom_line(data = expanded_controls_df,
              aes(x = Time, y = Mean_OD,
                  group = interaction(Metabolite, Experiment)),
              linewidth = 0.3, color = "black", linetype = "dashed") +
    facet_wrap(~ Organism, scales = "free_y", ncol = 2) +
    labs(title = "Growth Curves by Organism",
         subtitle = "Coloured lines = treatments (Sphingosine/Harmane); black dashed = controls (DBM/DM)",
         x = "Time (hours)", y = "OD") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = compound_colors) +
    scale_fill_manual(values = compound_colors)
}

# -----------------------------------------------------------------------------
# AUC calculation (individual replicates)
# -----------------------------------------------------------------------------
calculate_auc_individual <- function(raw_data) {
  cat("\n=== CALCULATING INDIVIDUAL AUC ===\n")
  auc_individual <- raw_data %>%
    group_by(Organism, Metabolite, Compound, Concentration_label, 
             Replicate, Type, Control_Code, Experiment) %>%
    arrange(Time) %>%
    summarise(
      AUC = {
        if(n() >= 2) {
          valid_times <- Time[!is.na(OD)]
          valid_ods <- OD[!is.na(OD)]
          if(length(valid_times) >= 2) {
            sum(diff(valid_times) * (valid_ods[-1] + valid_ods[-length(valid_ods)]) / 2)
          } else NA_real_
        } else NA_real_
      },
      Max_OD = max(OD, na.rm = TRUE),
      Final_OD = last(OD[!is.na(OD)]),
      .groups = 'drop'
    )
  return(auc_individual)
}

# -----------------------------------------------------------------------------
# Statistical tests using the metabolite-to-control mapping
# -----------------------------------------------------------------------------
perform_statistical_tests <- function(auc_individual) {
  cat("\n=== STATISTICAL ANALYSIS ===\n")
  
  treatments <- auc_individual %>% filter(Type == "Treatment")
  controls <- auc_individual %>% filter(Type == "Control")
  
  results <- list()
  
  for (org in unique(treatments$Organism)) {
    cat("\n--- Organism:", org, "---\n")
    org_treat <- treatments %>% filter(Organism == org)
    
    for (met in unique(org_treat$Metabolite)) {
      # Get control code for this metabolite
      ctrl_code <- metabolite_to_control[met]
      if (is.na(ctrl_code)) {
        cat("  Skipping metabolite", met, "- no control mapping.\n")
        next
      }
      
      # For each concentration of this metabolite
      for (conc in unique(org_treat$Concentration_label[org_treat$Metabolite == met])) {
        comp_data <- org_treat %>% filter(Metabolite == met, Concentration_label == conc)
        ctrl_data <- controls %>% filter(Organism == org, Control_Code == ctrl_code)
        
        if (nrow(comp_data) >= 2 && nrow(ctrl_data) >= 2) {
          t_test_result <- tryCatch({
            t.test(comp_data$AUC, ctrl_data$AUC)
          }, error = function(e) NULL)
          
          if (!is.null(t_test_result)) {
            pooled_sd <- sqrt((sd(comp_data$AUC)^2 + sd(ctrl_data$AUC)^2) / 2)
            cohens_d <- (mean(comp_data$AUC) - mean(ctrl_data$AUC)) / pooled_sd
          } else {
            cohens_d <- NA
          }
          
          compound_name <- unique(comp_data$Compound)
          results[[paste(org, met, conc, sep="_")]] <- data.frame(
            Organism = org,
            Metabolite = met,
            Concentration = conc,
            Compound = compound_name,
            Control = ctrl_code,
            Treatment_N = nrow(comp_data),
            Control_N = nrow(ctrl_data),
            Treatment_Mean_AUC = mean(comp_data$AUC, na.rm = TRUE),
            Treatment_SD_AUC = sd(comp_data$AUC, na.rm = TRUE),
            Control_Mean_AUC = mean(ctrl_data$AUC, na.rm = TRUE),
            Control_SD_AUC = sd(ctrl_data$AUC, na.rm = TRUE),
            Difference = mean(comp_data$AUC, na.rm = TRUE) - mean(ctrl_data$AUC, na.rm = TRUE),
            P_value = ifelse(!is.null(t_test_result), t_test_result$p.value, NA),
            T_statistic = ifelse(!is.null(t_test_result), t_test_result$statistic, NA),
            Cohens_d = cohens_d,
            Significance = case_when(
              is.null(t_test_result) ~ "Test failed",
              t_test_result$p.value < 0.001 ~ "*** p < 0.001",
              t_test_result$p.value < 0.01 ~ "** p < 0.01",
              t_test_result$p.value < 0.05 ~ "* p < 0.05",
              TRUE ~ "ns"
            )
          )
        }
      }
    }
  }
  
  stats_summary <- bind_rows(results)
  cat("\n=== STATISTICAL TEST RESULTS ===\n")
  if (nrow(stats_summary) > 0) {
    print(stats_summary %>%
            select(Organism, Metabolite, Concentration, Treatment_Mean_AUC, Control_Mean_AUC,
                   P_value, Significance, Cohens_d) %>%
            mutate(across(where(is.numeric), ~round(., 3))))
  } else {
    cat("No statistical tests could be performed.\n")
  }
  return(stats_summary)
}

# -----------------------------------------------------------------------------
# Inhibition calculation
# -----------------------------------------------------------------------------
calculate_inhibition_with_stats <- function(auc_individual, stats_results) {
  cat("\n=== CALCULATING INHIBITION WITH STATISTICS ===\n")
  
  # 1. Summarise control AUCs per organism and control code
  control_means <- auc_individual %>%
    filter(Type == "Control") %>%
    group_by(Organism, Control_Code) %>%
    summarise(
      Control_AUC_mean = mean(AUC, na.rm = TRUE),
      Control_AUC_sd   = sd(AUC, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # 2. For each treatment row, determine the appropriate control code via metabolite
  treatment_data <- auc_individual %>%
    filter(Type == "Treatment") %>%
    mutate(
      Control_Code = metabolite_to_control[Metabolite]
    ) %>%
    filter(!is.na(Control_Code))
  
  # 3. Join with control means
  inhibition_data <- treatment_data %>%
    left_join(control_means, by = c("Organism", "Control_Code")) %>%
    mutate(
      Inhibition_AUC = (1 - (AUC / Control_AUC_mean)) * 100
    ) %>%
    filter(!is.na(Inhibition_AUC))
  
  # 4. Summarise inhibition per organism, metabolite, concentration
  inhibition_summary <- inhibition_data %>%
    group_by(Organism, Metabolite, Concentration_label, Compound, Control_Code) %>%
    summarise(
      Mean_Inhibition = mean(Inhibition_AUC, na.rm = TRUE),
      SD_Inhibition   = sd(Inhibition_AUC, na.rm = TRUE),
      SEM_Inhibition  = SD_Inhibition / sqrt(n()),
      N = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      Inhibition_Level = case_when(
        Mean_Inhibition >= 80 ~ "Strong (≥80%)",
        Mean_Inhibition >= 50 ~ "Moderate (50-79%)",
        Mean_Inhibition >= 20 ~ "Weak (20-49%)",
        Mean_Inhibition > 0   ~ "Very Weak (1-19%)",
        TRUE                  ~ "No Inhibition/Growth"
      )
    )
  
  # 5. Attach statistical results
  if (nrow(stats_results) > 0) {
    inhibition_summary <- inhibition_summary %>%
      left_join(
        stats_results %>% select(Organism, Metabolite, Concentration, P_value, Significance, Cohens_d),
        by = c("Organism", "Metabolite", "Concentration_label" = "Concentration")
      )
  }
  
  return(list(raw_inhibition = inhibition_data, summary_inhibition = inhibition_summary))
}

# -----------------------------------------------------------------------------
# Inhibition plot
# -----------------------------------------------------------------------------
plot_inhibition_with_stats <- function(inhibition_summary) {
  if (nrow(inhibition_summary) == 0) return(NULL)
  
  # Create significance labels
  inhibition_summary <- inhibition_summary %>%
    mutate(
      sig_label = case_when(
        is.na(P_value) ~ "ns",
        P_value < 0.001 ~ "***",
        P_value < 0.01  ~ "**",
        P_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      y_position = if_else(
        !is.na(SEM_Inhibition),
        Mean_Inhibition + SEM_Inhibition + 5,
        Mean_Inhibition + 5
      ),
      # Create a nice x-axis label
      x_label = paste0(Metabolite, "\n", Concentration_label)
    )
  
  ggplot(inhibition_summary, aes(x = x_label, y = Mean_Inhibition, fill = Metabolite)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = Mean_Inhibition - SEM_Inhibition, 
                      ymax = Mean_Inhibition + SEM_Inhibition),
                  width = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_text(aes(y = y_position, label = sig_label), vjust = 0, size = 5) +
    facet_wrap(~ Organism, ncol = 1) +
    labs(
      title = "Growth Inhibition by Metabolite",
      subtitle = "Error bars = SEM; * p<0.05, ** p<0.01, *** p<0.001",
      y = "Inhibition (%)", 
      x = "Metabolite and Concentration", 
      fill = "Metabolite"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = c("S" = "#56BC37", "H" = "#2C6E96"))
}

# -----------------------------------------------------------------------------
# Complete analysis pipeline
# -----------------------------------------------------------------------------
perform_complete_analysis <- function(raw_data, summary_data) {
  cat("\n=== COMPLETE ANALYSIS WITH STATISTICS ===\n")
  
  # 1. AUC
  cat("\n1. Calculating individual AUC values...\n")
  auc_individual <- calculate_auc_individual(raw_data)
  cat("\n=== AVAILABLE CONTROLS BY ORGANISM ===\n")
  auc_individual %>%
    filter(Type == "Control") %>%
    distinct(Organism, Control_Code, Metabolite) %>%
    arrange(Organism, Control_Code) %>%
    print()
  
  # 2. Statistical tests
  cat("\n2. Performing statistical tests...\n")
  stats_results <- perform_statistical_tests(auc_individual)
  
  # 3. Inhibition
  cat("\n3. Calculating inhibition...\n")
  inhibition_results <- calculate_inhibition_with_stats(auc_individual, stats_results)
  
  # 4. Save files
  write.csv(auc_individual, "auc_individual_values.csv", row.names = FALSE)
  if (nrow(stats_results) > 0) write.csv(stats_results, "statistical_test_results.csv", row.names = FALSE)
  if (nrow(inhibition_results$raw_inhibition) > 0) write.csv(inhibition_results$raw_inhibition, "inhibition_detailed.csv", row.names = FALSE)
  if (nrow(inhibition_results$summary_inhibition) > 0) write.csv(inhibition_results$summary_inhibition, "inhibition_summary.csv", row.names = FALSE)
  
  # 5. Plot
  inhibition_plot <- plot_inhibition_with_stats(inhibition_results$summary_inhibition)
  
  return(list(
    auc_individual = auc_individual,
    stats_results = stats_results,
    inhibition_results = inhibition_results,
    inhibition_plot = inhibition_plot
  ))
}

# =============================================================================
# RUN THE ANALYSIS
# =============================================================================
cat("\n=== STARTING COMPLETE ANALYSIS ===\n\n")

processed <- create_growth_curves(data)

if (nrow(processed$summary_data) > 0) {
  
  # Growth curve plots
  cat("\n=== GENERATING GROWTH CURVES ===\n")
  tryCatch({
    gc_plot <- plot_growth_curves_faceted(processed$summary_data)
    print(gc_plot)
    ggsave("growth_curves_faceted.svg", gc_plot, width = 6, height = 9)
    cat("Growth curves (faceted by organism and compound) saved.\n")
  }, error = function(e) cat("Error in growth curves plot:", e$message, "\n"))
  
  # Full analysis with statistics
  analysis_results <- perform_complete_analysis(processed$raw_data, processed$summary_data)
  
  if (!is.null(analysis_results$inhibition_plot)) {
    tryCatch({
      print(analysis_results$inhibition_plot)
      ggsave("inhibition_plot.svg", analysis_results$inhibition_plot, width = 6, height = 8)
      cat("Inhibition plot saved as 'inhibition_plot.svg'\n")
    }, error = function(e) cat("Error in inhibition plot:", e$message, "\n"))
  }
  
  cat("\n=== FINAL SUMMARY ===\n")
  print(analysis_results$inhibition_results$summary_inhibition %>%
          select(Organism, Metabolite, Concentration_label, Mean_Inhibition, 
                 SEM_Inhibition, Inhibition_Level, Significance))
  
} else {
  cat("No data was processed. Please check your data structure.\n")
}
             
