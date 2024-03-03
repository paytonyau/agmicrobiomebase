#!/bin/bash

# Install Anaconda / Miniconda
# To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda, then create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

# Install Sequence Quality Check Packages
# For the installation of QIIME2, Anaconda or Miniconda can be utilised as they offer a self-contained environment and package manager. Begin by downloading and installing Miniconda, then create a new environment specifically for QIIME2. This approach simplifies the management of your QIIME2 installation and its dependencies.

conda create --name qc
conda install -c bioconda fastqc multiqc

# Install Trimmomatic
conda install -c conda-forge trimmomatic

# Install FIGARO
# FIGARO is a software that assists in estimating the truncation parameters for the QIIME2 DADA2 plugin. A pre-print detailing its usage is readily available.
git clone https://github.com/Zymo-Research/figaro.git
cd figaro

# Install QIIME2
# Run the workflow in a specific conda environment, which ensures the correct version of the Python required packages are being used for QIIME2.