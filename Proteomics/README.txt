Proteomics Analysis Script

This R script performs differential abundance analysis on proteomics CSV data comparing a mixed condition against a control condition for a specified species.

It filters the input by species and peptide count, handles missing abundance values, normalizes each sample by column sum, averages triplicate replicates, computes ratio and Log2 fold change (mix vs control), and runs per-protein t-tests with Benjamini-Hochberg correction. Results are flagged as significant when adjusted p < 0.05 and |Log2FC| > 1.

Input: a proteomics CSV file, mix condition, control condition, and species name.

Output: two CSV files — full results and significant results — plus a printed summary.

Requirements: R with dplyr, tidyr, stringr, and tibble.

Note: expects three replicates per condition (_R1, _R2, _R3). Usage examples for Gardnerella vaginalis and Lactobacillus crispatus are included at the bottom of the script.


______

Proteomic Shift Index Calculation Script

This R script processes all CSV files in a folder, parses protein descriptions into separate fields (locus tag, protein, gene, etc.), assigns labels, filters by p-value ≤ 0.05, splits proteins into upregulated (Log2FC > 0) and downregulated (Log2FC < 0), and calculates the Proteomic Shift Index (PSI) for each group.

PSI is computed as the mean or median Log2FC of each group multiplied by that group's proportion of total significant proteins. Both mean- and median-based PSI values are produced, along with combined totals.

Input: a folder of proteomics CSV files (path set via setwd()).

Output: intermediate _desc_filtered.csv and _labeled.csv files per input, plus a combined proteomic_shift_index_summary.csv.

Requirements: R with dplyr, tidyr, and stringr.

Note: expects p_value and Log2FC columns (with fallback name detection). Change the setwd() path before running.
