#!/bin/bash

# 02/04/2024
# List of required packages/software:
# 1. Anaconda or Miniconda
# 2. fastqc
# 3. multiqc
# 3. Trimmomatic
# 4. 16s-step01-figaro.yml
# 5. FIGARO
# 6. QIIME2

# To perform a sequencing analysis on the 16S amplicon data from the UK Crop Microbiome Cryobank, it is necessary to prepare a Linux environment with specific pre-installed packages. These packages, crucial for processing sequencing data, encompass a range of tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Prior to initiating the analysis, ensure that your system is equipped with all the necessary packages to facilitate a seamless and successful analysis.

# Install Anaconda / Miniconda
# To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda and create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

# Install sequence quality check packages
# For the installation of QIIME2, Anaconda or Miniconda can be utilised as they offer a self-contained environment and package manager. Begin by downloading and installing Miniconda, then create a new environment specifically for QIIME2. This approach simplifies the management of your QIIME2 installation and its dependencies.

# Create a new environment:
conda create --name qc

# Activate the conda environment
conda activate qc

# Install fastqc & multiqc
conda install -c bioconda fastqc multiqc

# Install Trimmomatic
conda install -c conda-forge trimmomatic
# Note that Java may need to be installed before the run.

# Install FIGARO
# FIGARO is a software that assists in estimating the truncation parameters for the QIIME2 DADA2 plugin. A pre-print detailing its usage is readily available. The detailed process for installing the software can be found on the author’s GitHub page. 

wget https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-step01-figaro.yml
conda env create -n figaro -f figaro.yml
git clone https://github.com/Zymo-Research/figaro.git
unzip master.zip
rm master.zip
cd figaro-master/figaro
chmod 755 *.py

# Install QIIME2 
# Run the workflow in a specific conda environment, which makes sure the correct version of the Python required packages are being used for QIIME2.
# Please download and follow the instruction from https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment