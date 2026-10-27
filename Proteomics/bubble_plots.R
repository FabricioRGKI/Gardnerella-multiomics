#Bubble plots for Differential Protein Abundance

#Example using the coculture mix LcFv with control Fv (Regulated proteins of Fannyhessea vaginae in the coculture L. crispatus + F. vaginae)

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Read the CSV file
data <- read.csv("coculture_vs_monoculture_control_significant_short_labeled.csv", stringsAsFactors = FALSE)

# Label top 21 by combined significance score
data <- data %>%
  mutate(significance_score = -log10(adj_p_value) * abs(Log2FC),
         Label = ifelse(rank(-significance_score) <= 20, label_edit, ""))

# Create the plot
ggplot(data, aes(x = Log2FC, y = -log10(adj_p_value))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_point(aes(size = mix_mean, 
                 fill = ifelse(Log2FC > 0, "Upregulated", "Downregulated")),
             alpha = 0.5,
             stroke = 0.5,  # Border thickness
             shape = 21,     # Shape that supports fill and color separately
             color = "black") +  # Black border for all points) +) +
  geom_text_repel(aes(label = Label),
                  size = 5,
                  max.overlaps = 50,
                  box.padding = 0.5) +
  scale_size_continuous(range = c(1, 20), 
                        name = "Mean normalized\nprotein\nabundance") +
  scale_fill_manual(values = c("Downregulated" = "#b0d5e0", 
                               "Upregulated" = "#4497AF"),
                    name = "Regulation") +
  labs(
    title = expression( ~ italic("Fannyhessea vaginae") ~ "protein changes in the mix LcFv"),
    subtitle = "Top 20 most relevant proteins labeled",
    x = "Log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal()

ggsave("Fv_LcFv.svg", width = 6, height = 6)
