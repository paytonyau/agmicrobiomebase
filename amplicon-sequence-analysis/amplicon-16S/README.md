# UK CROP MICROBIOME CRYOBANK (16S amplicon sequencing analysis)

Below is the scripting documents for each of the process, we included both markdown and `.sh` documents.

16s-sequence-process.md: This markdown file explains the entire process for the 16s amplicon sequencing.
16s-step01-startup.sh: A .sh file that includes the necessary packages to be installed for the environment before analysis.
16s-step02-processing.sh: These scripts handle data preparation during the initial stage.
16s-step03-qiime2.sh: This script is used for sequencing data analysis.

We have organised all the outcomes generated from Qiime2 and Figaro, where is it from the 16s sequencing process, into separate sub-folders.

1. [Qiime2]DADA2_outcomes: These outcomes were processed by DADA2 and are now ready for mapping to reference databases (such as GreenGenes and Silva)
   - 428_228_220_rep-seqs.qza
   - 428_228_220_stats.qza
   - 428_228_220_table.qza
   
2. [Qiime2]GreenGenes_13_8: These outcomes correspond to GreenGenes version 1 (13.8)
   - 428_228_220_rep-seqs_gg-13-8-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_table_gg-13-8-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_taxonomy_gg-13-8.qza
   - 428_228_220_taxonomy_gg-13-8.qzv
   - level_counts_by_group_gg1.csv
   
3. [Qiime2]GreenGenes2_2022_10: These outcomes relate to GreenGenes version 2 (2022.10)
   - 428_228_220_rep-seqs_gg_2022_10-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_table_gg_2022_10-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_taxonomy_gg_2022_10.qza
   - level_counts_by_group_gg2.csv
   
4. [Qiime2]Silva_138: These outcomes pertain to Silva version 138
   - 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_taxonomy_silva138.qza
   - level_counts_by_group_silva138.csv
P.S. The numbers 428, 228, and 220 indicate the longest biologically meaningful read for V3-V4, the nucleotide length of forward reads, and the nucleotide length of reverse reads.

5. [Qiime2]Silva_138_233_226_3_2_Fig_4B: These outcomes specifically address the use of trimming length identification within Silva
   - denoising_stats_233_226_3_2.qza
   - rep-seqs_233_226_3_2.qza
   - table_233_226_3_2.qza
   - rep-seqs_233_226_3_2-with-phyla-no-mitochondria-no-chloroplast.qza
   - table_233_226_3_2-with-phyla-no-mitochondria-no-chloroplast.qza
   
6. [Result]FIGARO: These results were generated from FIGARO
   - figaro.1693245011762037.log
   - forwardExpectedError.png
   - reverseExpectedError.png
   - trimParameters.json