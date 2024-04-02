# Preparation and Analysis Scripts for Microbiome Data

This directory contains R scripts for three analyses for the purpose of the data optimisation and standardisation. Each analysis serves a unique purpose and is crucial for ensuring the quality and integrity of the data. Here’s a brief overview of the three analyses:
1.  **Trimming and Optimal Sequence Length Selection**: Determining the longest sequence length that maximizes quality and minimizes error rates.
2.  **Comparative Analysis of Taxonomic Levels of Reference Databases**: Aggregating matched hits by names to consolidate ASVs and analyze taxonomic levels from Kingdom to Species.
3.  **Batch Effect Correction and Normalisation**: Correcting for batch effects to ensure consistent and accurate data across different sequencing runs.

## Contents
-   **Case_00-prep.R**: This script encompasses all preprocessing steps, including trimming, sequence length selection, comparative analysis of taxonomic levels, and batch effect correction.
-   **case00A-length_dist_ref_databases.Rmd**: The R Markdown file for the analysis of sequence length distributions across reference databases.
-   **case00A-length_dist_ref_databases.md**: The markdown documentation detailing the methodology and results of the sequence length distribution analysis.


-   **case00B-Batch_effects.Rmd**: The R Markdown file for the analysis and correction of batch effects in sequencing data.
-   **case00B-Batch_effects.md**: The markdown documentation explaining the batch effect correction process and its results.

-   **case00A-length_dist_ref_databases_files/figure-markdown_strict**: A folder containing figures generated from the sequence length distribution analysis.
-   **case00B-Batch_effects_files/figure-markdown_strict**: A folder containing figures generated from the batch effect analysis.


