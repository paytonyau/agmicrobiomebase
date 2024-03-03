#!/bin/bash

# STEP 1
# Examine All the FastQ Files Using FastQC and Consolidate the Results into a Single Report with the Help of MultiQC
# This procedure is designed to visualize the quality of sequences and any remaining adaptors in each individual FastQ file.
conda activate qc
fastqc *.fastq.gz
multiqc .

# STEP 2
# Sequencing Reads Trimming for FIGARO
# To maintain consistency as required by the program, sequence reads should be trimmed to a length of 248 nt, and any reads shorter than 248 nt should be discarded. This process can be accomplished using tools such as Cutadapt or Trimmomatic. In this case, Trimmomatic was used. Please note that the trimmed sequences are then ready for further analysis with FIGARO.
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248

# We use 300 nt long as for the average merged length as ITS1 region has a roughly ranging from 220-520nt.
python figaro.py -i (path) -o (path) -f 20 -r 20 -a 300 -m 12

# Generate a Manifest.txt for QIIME2 Artifact
# Generate a manifest.txt file containing the absolute file paths of the sequencing data, and then convert it to QIIME2-compatible format (manifest.csv).
ls -d "$PWD"/* > manifest.txt

# Input and output file paths
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
