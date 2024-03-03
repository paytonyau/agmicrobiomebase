#!/bin/bash

# Activate the QIIME2 environment
conda activate qiime2-2023.5

# Import the paired-end raw reads as a QIIME2 artifact
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv \
--output-path demux \
--input-format PairedEndFastqManifestPhred33

# Obtain sequencing data quality information
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

# Trim primers using `Cutadapt` from the QIIME2 plugin
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-cores 4 \
  --p-adapter-f AGATGTGTATAAGAGACAG \
  --p-adapter-r AGATGTGTATAAGAGACAG \
  --p-front-f GAACCWGCGGARGGATCA \
  --p-front-r GCTGCGTTCTTCATCGATGC \
  --p-discard-untrimmed \
  --o-trimmed-sequences primer-trimmed-demux.qza \
  --verbose \
  &> primer_trimming.log

# DADA2 Denoising
qiime dada2 denoise-paired \
--i-demultiplexed-seqs primer-trimmed-demux.qza \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-max-ee-f 2 \
--p-max-ee-r 2 \
--p-n-threads 4 \
--o-representative-sequences rep-seqs-its.qza \
--o-table table-its.qza \
--o-denoising-stats stats-its.qza

# Taxonomic Classification
qiime feature-classifier classify-sklearn \
--i-classifier /mnt/shared/scratch/pyau/qiime2-ref/unite-9-dynamic-s-all-29.11.2022-Q2-2023.5.qza \
--i-reads rep-seqs-its.qza \
--o-classification taxonomy-its.qza

# Unrooted Phylogenetic Tree Construction
qiime phylogeny align-to-tree-mafft-fasttree \  
  --i-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza \  
  --output-dir phylogeny-align-to-tree-mafft-fasttree