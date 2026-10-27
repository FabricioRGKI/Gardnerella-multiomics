#STRING Rank Score Input Preparation

rm(list = ls())
setwd("your/path")
getwd()

# Read CSV
df <- read.csv("condition-x_locus_tag_fold-change_p-value.csv", stringsAsFactors = FALSE)

# Combined rank score
df$rank_score <- -log10(df$p_value) * sign(df$Log2FC)

# Two-column STRING input
string_input <- df[, c("locus_tag", "rank_score")]

# Save
write.csv(string_input, "string_rank_input.csv", row.names = FALSE, quote = FALSE)
