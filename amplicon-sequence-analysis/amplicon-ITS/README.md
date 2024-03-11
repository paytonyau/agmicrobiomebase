# UK CROP MICROBIOME CRYOBANK (ITS amplicon sequencing analysis)

Below is the scripting documents for each of the process, we included both markdown and `.sh` documents for giving an example for the ITS sequencing analysis. An .md file is a text file that uses Markdown syntax to format text elements, such as headings, lists, links, and images; it is a lightweight and easy-to-read markup language that can be converted to HTML or other formats. A .sh file is a text file that contains commands to be executed by a Unix shell or a specific interpreter, and it is used to automate tasks in Unix-like operating systems.


ITS-sequencing-process.md: This markdown file explains the entire process for the 16s amplicon sequencing.

ITS-step01-startup.sh: A .sh file that includes the necessary packages to be installed for the environment before analysis.
ITS-step02-processing.sh: These scripts handle data preparation during the initial stage.
ITS-step03-qiime2.sh: This script is used for sequencing data analysis


We have organised the key outcomes generated from Qiime2, where is it from the ITS1 sequencing process, into separate sub-folders.
1. [Qiime2]DADA2_outcomes - These outcomes were processed by DADA2 and are now ready for mapping to reference database, here for UNITE Version 9
- rep-seqs-its.qza
- stats-its.qza
- table-its.qza

2. [Qiime2]UNITE_9 - These outcomes correspond to UNITE Version 9
 - rep-seqs-its-fungi-with-phyla-no-mitochondria-no-chloroplast.qza
 - table-its-fungi-with-phyla-no-mitochondria-no-chloroplast.qza
 - taxonomy-fungi.qza