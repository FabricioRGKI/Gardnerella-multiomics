Pairwise Percent Identity Analysis for Orthologous Proteins

aa_identity contains an R script that calculates pairwise percent identity between orthologous protein sequences from two entities.

The script takes named protein sequence vectors for each entity, finds proteins shared between entity pairs, performs global pairwise alignments, and computes percent identity (identical positions / alignment length × 100). It runs all-vs-all comparisons across entities, then produces summary statistics (mean, median, min, max, SD, number of comparisons) and exports results as CSV files. It also flags highly conserved proteins (>95% identity) and divergent proteins (<80% identity).

Input: named character vectors of protein sequences per entity, combined into a list.

Output: per-comparison CSV files with identity and length data, plus a summary CSV.

Requirements: R with Biostrings, pwalign, and tidyverse.

Note: alignments are global, so identity reflects full-length comparison. The included example uses two bacterial strains and can be extended to any number by adding entries to the strain list.
