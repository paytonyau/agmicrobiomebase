# UK Crop Microbiome Cryobank
## Agmicrobiomebase

The UK Crop Microbiome Cryobank integrates genomic (DNA) data with a cryobank collection of samples for the soil microbiomes of the UK major crop plant systems. For this project, the microbiomes are from the rhizosphere (the soil surrounding the crop plant roots) and from bulk soil (soil outside the rhizosphere). The Cryobank provides a facility for researchers to source data and samples, including cryo-preserved microbial material and genomic and metagenomic sequences from different soil microbiome environments.

The data catalogue for this project can be accessed at https://agmicrobiomebase.org/.
This repository contains the scripts used to process the fastq files from the 16S and ITS amplicon Illumina sequencing and for the metagenomic Illumina sequencing. 

### **Table of contents**
--------------------------------
### *A. Amplicon Sequencing*
#### A1. 16s Amplicon Gene (V3-V4) Sequencing Data Analysis

[01. 16s Sequencing Analysis (markdown)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/16s_sequence_process.md)

[02. 16s - Down Stream Analysis (R markdown)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/16s_Rmarkdown.Rmd)
[*16s - Down Stream Analysis (pdf format)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/16s_Rmarkdown.pdf)

[03. Merged Sequence Length Distribution for Figure 4B (R script)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/R-scripts-figures/Fig_4B_merged_length_distribution.R) 

[04. Reference Databases Comparison for Figure 5 (R script)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/R-scripts-figures/Fig_5_ref_databases_comparison.R) 

------------------------------------
##### 16S Documents/Data/Results
[05. DADA2 analysis outcomes (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BQiime2%5DDADA2_outcomes)

[06. Silva [v138] (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BQiime2%5DSilva_138)

[07. GreenGene1 [v13.8] (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BQiime2%5DGreenGenes_13_8)

[08. GreenGenes2 [v2022.10] (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BQiime2%5DGreenGenes2_2022_10)

[09.  Silva [v138] for Figure 4B (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BQiime2%5DSilva_138_233_226_3_2_Fig_4B) *for the merged sequence length distribution analysis

[10. Figure 4B* sequence length distribution data (tsv file in zip format)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/R-scripts-figures/16s_length_distribution.zip) 

[11. FIGARO (Result)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/16s/%5BResult%5DFIGARO)

[12. Meta-table](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/16s/meta-table.txt)

-----------------------------------------------------------
#### A2. ITS Amplicon Gene (ITS1) Sequencing Data Analysis

[13. ITS Sequencing Analysis (markdown)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/ITS/ITS_sequencing_process.md)

[14. ITS Down Stream Analysis (R markdown)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/ITS/ITS_Rmarkdown.Rmd)
[*ITS Down Stream Analysis (pdf format)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/ITS/ITS_Rmarkdown.pdf)

------------------------------------
##### ITS Documents/Data

[15. DADA2 analysis outcomes (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/ITS/%5BQiime2%5DDADA2_outcomes)

[16. UNITE [v9.0] (Qiime2 data structure)](https://github.com/paytonyau/agmicrobiomebase/tree/main/Amplicon/ITS/%5BQiime2%5DUNITE_9)

[17. Meta-table (txt)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/ITS/ITS-meta-table.txt)

[19. Other R programming scripts for the supplementary figures](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/R-scripts-figures/)

------------------------------------
#### A3. Amplicon Sequencing Fastq files upload

We used Oilseed rape as an example to demonstrate the mapping process for the European Nucleotide Archive (ENA) upload

[20. Fastq checklist mapping (R mark down)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/Fastq_checklist_mapping/fastq_checklist_mapping.Rmd)
[*Fastq checklist mapping (pdf format)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/Fastq_checklist_mapping/fastq_checklist_mapping.pdf)

------------------------------------
##### Documents/Data for the fastq files upload
[21. MD5 CheckSum (txt)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/Fastq_checklist_mapping/md5.txt)

[22. Plant Checklist (tsv)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/Fastq_checklist_mapping/Checklist_GSC-MIxS_16Samplicons_OR_TESTv1.tsv)

[23. Receipt from the Plant Checklist uploaded (txt)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Amplicon/Fastq_checklist_mapping/Webin-accessions-2023-12-07T15_42_52.222Z_OR.txt)

------------------------------------
### *B. Shotgun Metagenomics Sequencing*
####  B1. Shotgun Metagenomics Sequencing Data Analysis

------------------------------------
##### Documents/Data for the fastq files upload
[51. MD5 CheckSum (txt)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Shotgun/Fastq_checklist_mapping/md5.txt)

[52. Meta-table (txt)](https://github.com/paytonyau/agmicrobiomebase/blob/main/Shotgun/Fastq_checklist_mapping/meta-table.txt)