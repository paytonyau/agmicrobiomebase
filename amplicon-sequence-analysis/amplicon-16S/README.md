# UK CROP MICROBIOME CRYOBANK (16S amplicon sequencing analysis)

Below are the scripting documents for each of the process; we included both markdown (`.md`) and `.sh` documents. An .md file is a text file that uses Markdown syntax to format text elements, such as headings, lists, links, and images; it is a lightweight and easy-to-read markup language that can be converted to HTML or other formats. A .sh file is a text file that contains commands to be executed by a Unix shell or a specific interpreter, and it is used to automate tasks in Unix-like operating systems.


16s-sequence-process.md: This markdown file explains the entire process for the 16s amplicon sequencing.
16s-step01-startup.sh: A .sh file that includes the necessary packages to be installed for the environment before analysis.
16s-step02-processing.sh: These scripts handle data preparation during the initial stage.
16s-step03-qiime2.sh: This script is used for sequencing data analysis.


We have organised the key outcomes generated from Qiime2, which is it from the 16s sequencing process, into separate sub-folders to help with using our sequencing analysis process.

1. [Qiime2]DADA2_outcomes: These outcomes were processed by DADA2 and are now ready for mapping to reference databases (such as GreenGenes and Silva)
   - 428_228_220_rep-seqs.qza
   - 428_228_220_stats.qza
   - 428_228_220_table.qza
      
2. [Qiime2]Silva_138: These outcomes pertain to Silva version 138
   - 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
   - 428_228_220_taxonomy_silva138.qza
   - level_counts_by_group_silva138.csv

P.S. The numbers 428, 228, and 220 indicate the longest biologically meaningful read for V3-V4, the nucleotide length of forward reads, and the nucleotide length of reverse reads.

