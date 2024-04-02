
# 16S Amplicon Sequencing Analysis

In the  **16S Amplicon Sequencing Analysis**  sub-folder, contains information associated with crop plants using 16S rRNA gene sequencing data. 

## Overview

-   **16S Sequencing Process**: The  `16s-sequence-process.md`  file provides a comprehensive guide to the entire workflow for 16S amplicon sequencing. It covers data preparation, quality control, and downstream analysis steps.
    
-   **Scripts and Shell Files**:
    
    -   `16s-step01-startup.sh`: A shell script that sets up the necessary environment by installing required packages.
    -   `16s-step01-figaro.yml`: : A YAML file for setting up the `figaro` conda environment with necessary packages.
    -   `16s-step02-processing.sh`: These scripts handle data preparation during the initial stage of analysis.
    -   `16s-step03-qiime2.sh`: This script performs the actual sequencing data analysis using Qiime2.

## Key Outcomes

We have organised the key outcomes generated from Qiime2 during the 16S sequencing process into separate sub-folders:

1.  **[Qiime2]DADA2_outcomes**:
    
    -   These outcomes were processed by DADA2 and for reference databases mapping (such as GreenGenes and Silva).
    -   Files:
        -   `428_228_220_rep-seqs.qza`: Processed representative sequences.
        -   `428_228_220_stats.qza`: Statistical summary of the data.
        -   `428_228_220_table.qza`: Processed feature table.
2.  **[Qiime2]Silva_138**:
    
    -   These outcomes correspond to Silva version 138, post-DADA2 processing.
    -   Files:
        -   `428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza`: Fungal representative sequences.
        -   `428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza`: Fungal feature table.
        -   `428_228_220_taxonomy_silva138.qza`: Taxonomic assignments for fungal sequences.
