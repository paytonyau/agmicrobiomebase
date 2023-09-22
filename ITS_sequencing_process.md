## QIIME2 Data Processing Pipeline for ITS (first draft)

  

### Install QIIME2

To run the workflow in a specific conda environment, ensuring the correct version of required Python packages for QIIME2, refer to the following documentation: [Install QIIME 2 within a Conda environment](https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment)

### Step 1: Sequence Quality Control

#### FastQC and MultiQC

Visualise sequence quality and adapters using FastQC for each individual fastq file, then integrate the reports using MultiQC.

Activate the QC environment:
`conda activate qc`
Run FastQC and MultiQC:
`fastqc *.fastq.gz`
`multiqc .`

  #### Step 2: Sequencing Reads Trimming for FIGARO

Trim sequence reads to a length of 248nt and discard reads shorter than 248nt using Trimmomatic for FIGARO.
```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248
```

We use 300nt long as for the average merged length as ITS1 region has a roughly ranging from 220-520nt[^1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309391)

`python figaro.py -i (path) -o (path) -f 20 -r 20 -a 300 -m 12`

Please replace `(path)` with actual paths, and ensure that you have the required software and resources available for each step of the pipeline.
### Generate a manifest.txt for QIIME2 artifact


Generate a manifest.txt file containing the absolute file paths of the sequencing data, and then convert it to QIIME2-compatible format (manifest.csv).

##### Generate a manifest.txt for QIIME2 artifact (medification required to fit for the QIIME2 format)

```ls -d "$PWD"/* > manifest.txt```

The prescribed format entails the inclusion of "sample-id, absolute file-path, direction". The following script can be utilised to generate the necessary format, you could copy the script below and paste it into a `.sh` file.

```
#!/bin/bash
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
```
### Step 3: Import FASTQs as QIIME2 Artifact

Activate the QIIME2 environment:

`conda activate qiime2-2023.5`
Import the paired-end raw reads as a QIIME2 artifact:
```
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv \
--output-path demux \
--input-format PairedEndFastqManifestPhred33
```
### Summarize FASTQs
Obtain sequencing data quality information:
```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

##### Trim primers using Cutadapt QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the ITS1 region. DADA2's chimaera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimaeras to be identified.

The primes used in this dataset for ITS1 region
###### TS1-Fl2 Primer-F	
5’ 
TCGTCGGCAGCGTC
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*GAACCWGCGGARGGATCA*  (Primer sequence)
3’

###### ITS2 Primer-R	
5’ 
GTCTCGTGGGCTCGG
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*GCTGCGTTCTTCATCGATGC* (Primer sequence)
  3’ 

```
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

### Step 4: DADA2 Denoising and Taxonomic Classification
Denoise and perform taxonomic classification using DADA2 and the UNITE reference classifier.
##### DADA2 Denoising:
```
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

```
Taxonomic Classification:  

```
qiime feature-classifier classify-sklearn \
--i-classifier /mnt/shared/scratch/pyau/qiime2-ref/unite-9-dynamic-s-all-29.11.2022-Q2-2023.5.qza \
--i-reads rep-seqs-its.qza \
--o-classification taxonomy-its.qza
```

### Step 5: Taxonomy-Based Filtering
Filter tables and sequences based on taxonomy to exclude unwanted taxa (e.g., mitochondria and chloroplast).

Filter Table:
```
qiime taxa filter-table \
--i-table table-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-table table-its-with-phyla-no-mitochondria-no-chloroplast.qza
```


Filter Sequences:


```
qiime taxa filter-seqs \
--i-sequences rep-seqs-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza
```