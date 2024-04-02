# ITS Amplicon Sequencing Analysis

In this repository, it contains Internal Transcribed Spacer (ITS) amplicon sequencing data for wheat, which provides insights into the fungal communities inhabiting agricultural ecosystems.

## Overview

-   **ITS Sequencing Process**: The  `ITS-sequencing-process.md`  file explains the entire workflow for 16S amplicon sequencing. It covers data preparation, quality control, and downstream analysis steps.
    
-   **Shell Script and the meta-table**:
    
    -   `ITS-step01-startup.sh`: A shell script that sets up the necessary environment by installing required packages.
    -   `ITS-step02-processing.sh`: These scripts handle data preparation during the initial stage of analysis.
    -   `ITS-step03-qiime2.sh`: This script performs the actual sequencing data analysis using Qiime2.
    -   `ITS-meta-table.txt`: A metadata table that contains information about the samples used in the analysis. It includes details such as sample IDs, sample types, and other relevant characteristics.

## Key Outcomes

We have organised the key outcomes generated from Qiime2 during the ITS1 sequencing process into separate sub-folders:

1.  **[Qiime2]DADA2_outcomes**:
    
    -   These outcomes were processed by DADA2 and are now ready for mapping to reference databases (specifically UNITE Version 9).
    -   Files:
        -   `rep-seqs-its.qza`: Processed representative sequences.
        -   `stats-its.qza`: Statistical summary of the data.
        -   `table-its.qza`: Processed feature table.
2.  **[Qiime2]UNITE_9**:
    
    -   These outcomes correspond to UNITE Version 9.
    -   Files:
        -   `rep-seqs-its-fungi-with-phyla-no-mitochondria-no-chloroplast.qza`: Fungal representative sequences.
        -   `table-its-fungi-with-phyla-no-mitochondria-no-chloroplast.qza`: Fungal feature table.
        -   `taxonomy-fungi.qza`: Taxonomic assignments for fungal sequences.
