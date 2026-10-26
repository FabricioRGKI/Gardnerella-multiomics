Proteomics Analysis Script

This R script performs differential abundance analysis on proteomics CSV data comparing a mixed condition against a control condition for a specified species.

It filters the input by species and peptide count, handles missing abundance values, normalizes each sample by column sum, averages triplicate replicates, computes ratio and Log2 fold change (mix vs control), and runs per-protein t-tests with Benjamini-Hochberg correction. Results are flagged as significant when adjusted p < 0.05 and |Log2FC| > 1.

Input: a proteomics CSV file, mix condition, control condition, and species name.

Output: two CSV files — full results and significant results — plus a printed summary.

Requirements: R with dplyr, tidyr, stringr, and tibble.

Note: expects three replicates per condition (_R1, _R2, _R3). Usage examples for Gardnerella vaginalis and Lactobacillus crispatus are included at the bottom of the script.
