# Log2FC Bar Plot by Condition and Module

# ---- libraries ----
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# ---- load data ----
file_path <- "Gv_regulated_modules_all.csv"
df <- read_csv(file_path, show_col_types = FALSE)

# ---- clean condition into Lc / Li / Fv / Pb ----
df <- df %>%
  mutate(
    cond_short = case_when(
      grepl("Lc", condition) ~ "Lc",
      grepl("Li", condition) ~ "Li",
      grepl("Fv", condition) ~ "Fv",
      grepl("Pb", condition) ~ "Pb",
      TRUE ~ NA_character_
    ),
    cond_key = paste0(cond_short, "_", ifelse(regulation == "up", "up", "dn")),
    regulation = factor(regulation, levels = c("up", "down")),
    cond_key = factor(cond_key,
                      levels = c("Lc_up","Li_up","Fv_up","Pb_up",
                                 "Lc_dn","Li_dn","Fv_dn","Pb_dn"))
  )

# ---- color palette ----
cond_colors <- c(
  Lc_up = "#E3377E", Li_up = "#B02CBA", Fv_up = "#4497AF", Pb_up = "#56BC37",
  Lc_dn = "#f19ec1", Li_dn = "#e5a6ea", Fv_dn = "#b0d5e0", Pb_dn = "#b8e6aa"
)

# ============================================================
# KEY: one row per (protein x condition)
# ============================================================
# Order: regulation -> module -> protein -> condition
# This keeps the conditions of the same protein together.
df <- df %>%
  arrange(regulation, modules, protein, cond_key) %>%
  mutate(
    # unique label per row
    row_label = paste0(protein, "  [", cond_short, "]"),
    row_label = factor(row_label, levels = unique(row_label))
  )

# ============================================================
# MAIN PANEL: diverging bars, ONE bar per row
# ============================================================
p_main <- ggplot(df, aes(x = Log2FC, y = row_label, fill = cond_key)) +
  geom_col(width = 0.8) +                                  # no dodge
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.4) +
  facet_grid(regulation ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = cond_colors, name = "Condition") +
  labs(x = "Log2FC", y = NULL,
       title = "Regulated proteins: Log2FC by condition and module") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 7)
  )

# ============================================================
# LEFT STRIP: module per row
# ============================================================
p_mod <- ggplot(df, aes(x = 1, y = row_label, label = modules)) +
  geom_text(size = 2.6, hjust = 1, color = "grey20") +
  facet_grid(regulation ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(strip.text = element_blank(),
        plot.margin = margin(r = 2))

# ============================================================
# COMBINE
# ============================================================
p_final <- p_mod + p_main +
  plot_layout(widths = c(0.8, 5), guides = "collect")

# ---- save ----
ggsave("plot_bars_modules.pdf", p_final,
       width = 10, height = 10, units = "in", device = "pdf")
