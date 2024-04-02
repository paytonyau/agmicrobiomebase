#!/bin/bash

# 02/04/2024 - Payton Yau

# List of required files/software:
# 1. FastQC: Quality control tool for high throughput sequence data
# 2. MultiQC: Tool to aggregate results from bioinformatics analyses across many samples into a single report
# 3. Trimmomatic: Flexible read trimming tool for Illumina sequence data
# 4. FIGARO: Software for estimating the truncation parameters for the QIIME2 DADA2 plugin
# 5. QIIME2: Microbiome bioinformatics platform for flexible and reproducible microbial analysis
# 6. Input files: .fastq.gz files for analysis
# 7. Manifest file: Text file listing all the files for analysis
# 8. CSV file: CSV file generated from the manifest file for QIIME2

# STEP 1
# Examine all the fastq files using FastQC and consolidate the results into a single report with the help of MultiQC
# This procedure is designed to visualise the quality of sequences and any remaining adaptors in each individual fastq file.

fastqc *.fastq.gz
multiqc .

# STEP 2
# Sequencing reads trimming for FIGARO

# To maintain consistency as required by the program, sequence reads should be trimmed to a length of 248 nt, and any reads shorter than 248 nt should be discarded. This process can be accomplished using tools such as Cutadapt or Trimmomatic. In this case, Trimmomatic was used. Please note that the trimmed sequences are then ready for further analysis with FIGARO.

java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248

# Identify the optimal denoising setting in DADA2 by using FIGARO

# Following the trimming process, the Figaro program should be utilised. The program involves applying various read lengths to estimate the error rates for DADA2. The aim is to optimize the accuracy of the DADA2 algorithm and improve its performance in accurately identifying and correcting sequencing errors in the dataset.

# We use 428 nt as the longest biologically meaningful read and the setting would be  

python figaro.py -i (path) -o (path) -f 17 -r 21 -a 428 -m 18


# Generate a manifest.txt for QIIME2 artifact (modification required to fit the QIIME2 format)
ls -d "$PWD"/* > manifest.txt

# The prescribed format entails the inclusion of "sample-id, absolute file-path, direction". The following script can be utilised to generate the necessary format. You could copy the script below and paste it into a `.sh` file.

input_file="manifest.txt"
output_file="manifest.csv"

# Create the CSV header
echo "sample-id,absolute-filepath,direction" > "$output_file"

# Read the input file line by line
while IFS= read -r line; do
    # Extract the sample name, direction, and remove "_S123"
    sample=$(echo "$line" | sed 's#.*/##; s/\..*//; s/_L001_.*//; s/_\([0-9]\+\)$//; s/_S[0-9]\+$//')
    direction=$(echo "$line" | grep -o "_R[12]_")

    # Determine the direction and create the CSV lines
    if [ "$direction" = "_R1_" ]; then
        echo "$sample,$line,forward" >> "$output_file"
    elif [ "$direction" = "_R2_" ]; then
        echo "$sample,$line,reverse" >> "$output_file"
    fi
done < "$input_file"

echo "CSV file created: $output_file"
