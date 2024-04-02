#!/bin/bash

# 02/04/2024 - Payton Yau

# List of required files:
# 1. manifest.csv: CSV file that follows a specific format required by QIIME2
# 2. silva-138-99-nb-classifier.qza: Classifier for feature-classifier

# List of output files:
# 1. demux.qza: Output from the qiime tools import command
# 2. demux.qzv: Visualization of the demultiplexed data
# 3. primer-trimmed-demux.qza: Output from the qiime cutadapt trim-paired command
# 4. primer-trimmed-demux.qzv: Visualization of the primer-trimmed data
# 5. 428_228_220_rep-seqs.qza, 428_228_220_table.qza, 428_228_220_stats.qza: Outputs from the qiime dada2 denoise-paired command
# 6. 428_228_220_taxonomy_silva138.qza: Output from the qiime feature-classifier classify-sklearn command
# 7. 428_228_220_taxonomy_silva138.qzv: Visualization of the taxonomy
# 8. 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza: Filtered table from the qiime taxa filter-table command
# 9. 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza: Filtered sequences from the qiime taxa filter-seqs command
# 10. 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qzv: Visualization of the representative sequence
# 11. table-exported: Directory containing the exported BIOM file
# 12. sequences: Directory containing the exported FASTA file
# 13. picrust2_out_pipeline: Output directory from the picrust2_pipeline.py command

# Start of the script
# Activate the environment 
# You can activate this conda environment with this command (you may need to swap in source for conda if you get an error):

conda activate qiime2-2023.5

# Import FASTQs as QIIME 2 artifact
# To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an "artifact" (with the extension QZA). The first step is to import the paired-end raw reads as a QZA file. Note that you need to prepare the manifest.csv file following the QIIME2 requirement.

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33

# Summarise FASTQs
# We can also obtain the sequencing data quality information by executing the command below.

qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

# Trim primers using Cutadapt QIIME2 plugin
# Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the 16s V3-V4 region. DADA2's chimera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimeras to be identified.

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-cores 8 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-discard-untrimmed \
  --o-trimmed-sequences primer-trimmed-demux.qza \
  --verbose \
  &> primer_trimming.log

# Summarise FASTQs after the primers trimming
# You can run the demux summarize command after trimming the reads to get a report of the number of reads per sample and quality distribution across the reads. This generates a more basic output compared to FASTQC/MultiQC, but is sufficient for this step.

qiime demux summarize \
  --i-data primer-trimmed-demux.qza \
  --o-visualization primer-trimmed-demux.qzv

# Denoising the reads into amplicon sequence variants using DADA2
# Based on the projected outcomes of the FIGARO programme, the expected parameters - error rates (--p-max-ee-f & --p-max-ee-r) and the optimal truncation length (--p-trunc-len-f & --p-trunc-len-r) - can be utilised in this scenario to achieve the best denoising results.
# From Figaro, we have the trimming position of Forward = 245 and Reverse  = 241, minus the length of the primers (as we previously trimmed). Therefore, 

# Forward: 245 - 17 = 228 nt
# Reverse: 241 - 21 = 220 nt

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed-demux.qza \
  --p-trunc-len-f 228 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 8 \
  --o-representative-sequences 428_228_220_rep-seqs.qza \
  --o-table 428_228_220_table.qza \
  --o-denoising-stats 428_228_220_stats.qza

# Assigning taxonomy to ASVs by running taxonomic classification
# The command for running the taxonomic classification is one of the most time-consuming and memory-intensive commands in the SOP. If you encounter an error due to insufficient memory and cannot increase your memory usage, consider adjusting the --p-reads-per-batch option to a value lower than the default (which is dynamic, depending on sample depth and the number of threads), and try running the command with fewer jobs (for example, set --p-n-jobs to 1). Alternatively, you can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit-learn Python library, along with the SILVA 138 or Greengene2 database. The pre-trained classifiers can be obtained from Qiime2 data-resources.

qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/ref/silva-138-99-nb-classifier.qza \
  --i-reads 428_228_220_rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification 428_228_220_taxonomy_silva138.qza

# Visualising the taxonomy outcomes

qiime metadata tabulate \
  --m-input-file 428_228_220_taxonomy_silva138.qza \
  --o-visualization 428_228_220_taxonomy_silva138.qzv

# Filtering out contaminant and unclassified ASVs
# With taxonomy assigned to our ASVs, we can leverage this information to eliminate ASVs that are likely contaminants or noise, based on their taxonomic labels. Mitochondrial and chloroplast 16S sequences are common contaminants in 16S sequencing data and can be removed by excluding any ASV containing these terms in its taxonomic label. It may also be beneficial to exclude any ASV unclassified at the phylum level, as these sequences are more likely to be noise (e.g., potential chimeric sequences).
# Please note that if your data has not been classified against SILVA, you will need to modify 'P__' to a string that allows phylum-level assignments to be identified, or simply omit that line. If you are studying a poorly characterized environment where there is a good chance of identifying novel phyla, you might also want to omit that line.

qiime taxa filter-table \
  --i-table 428_228_220_table.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza

qiime taxa filter-seqs \
  --i-sequences 428_228_220_rep-seqs.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza

# Visualise the Representative Sequence
qiime metadata tabulate \
  --m-input-file 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza \
  --o-visualization 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qzv

# Functional Prediction - `PICRUSt2`
# Output `BIOMV210DirFmt` (BIOM) file
qiime tools export --input-path 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza  --output-path table-exported

# Output `DNASequencesDirectoryFormat` (FASTA) format file
qiime tools export --input-path 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza --output-path sequences

# Create and activate the environment for `PICRUSt2`
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2

# Activate the new conda environment
conda activate picrust2

# Run the program
picrust2_pipeline.py -s sequences/dna-sequences.fasta -i table-exported/feature-table.biom -o picrust2_out_pipeline -p 1