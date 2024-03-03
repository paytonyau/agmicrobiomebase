
# UK CROP MICROBIOME CRYOBANK (ITS amplicon sequencing analysis)

To perform a sequencing analysis on the ITS amplicon data from the UK Crop Microbiome Cryobank, it is necessary to prepare a Linux environment with specific pre-installed packages. These packages, crucial for processing sequencing data, encompass a range of tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Prior to initiating the analysis, ensure that your system is equipped with all the necessary packages to facilitate a seamless and successful analysis.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [QIIME2 - the cummunity developed suite](https://qiime2.org/)
- [FIGARO](https://github.com/Zymo-Research/figaro)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### 1. Install packages/programs for ITS amplicon sequencing data process

#### Install Anaconda / Miniconda

To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda and create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

#### Install sequence quality check packages

For the installation of QIIME2, Anaconda or Miniconda can be utilised as they offer a self-contained environment and package manager. Begin by downloading and installing Miniconda, then create a new environment specifically for QIIME2. This approach simplifies the management of your QIIME2 installation and its dependencies.

```conda create --name qc```

```conda install -c bioconda fastqc multiqc```

#### Install Trimmomatic

```conda install -c conda-forge trimmomatic```
Noted that Java may need to be installed before the run.

#### Install FIGARO
FIGARO is a software that assists in estimating the truncation parameters for the QIIME2 DADA2 plugin. A pre-print detailing its usage is readily available[^1].

[^1]:Weinstein, M. et al. (2019) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. bioRxiv DOI: 10.1101/610394; Sasada, R. et al. (2020) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. J. Biomol. Tech. 31, S2

```
git clone https://github.com/Zymo-Research/figaro.git
cd figaro
```

#### Install QIIME2

Run the workflow in a specific conda environment, which makes sure the correct version of the Python required packages are being used for QIIME2.
<https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment>

### STEP 1
#### Examine all the fastq files using FastQC and consolidate the results into a single report with the help of MultiQC
This procedure is designed to visualize the quality of sequences and any remaining adaptors in each individual fastq file.
``` conda activate qc ```

``` fastqc *.fastq.gz ```
``` multiqc . ```

### STEP 2
#### Sequencing reads trimming for FIGARO

To maintain consistency as required by the program, sequence reads should be trimmed to a length of 248 nt, and any reads shorter than 248 nt should be discarded[^2]. This process can be accomplished using tools such as Cutadapt or Trimmomatic. In this case, Trimmomatic was used. Please note that the trimmed sequences are then ready for further analysis with FIGARO.
[^2]: https://github.com/Zymo-Research/figaro/issues/37
```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248
```

We use 300 nt long as for the average merged length as ITS1 region has a roughly ranging from 220-520nt[^3].
[^3]:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5309391

`python figaro.py -i (path) -o (path) -f 20 -r 20 -a 300 -m 12`

Please replace `(path)` with actual paths, and ensure that you have the required software and resources available for each step of the pipeline.

### Generate a manifest.txt for QIIME2 artifact

Generate a manifest.txt file containing the absolute file paths of the sequencing data, and then convert it to QIIME2-compatible format (manifest.csv).
#### Generate a manifest.txt for QIIME2 artifact (medification required to fit for the QIIME2 format)

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
#### Summarise FASTQs
Obtain sequencing data quality information:
```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

#### Trim primers using `Cutadapt` from the QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the ITS1 region. DADA2's chimaera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimaeras to be identified.

The primers used in this dataset for ITS1 region
##### TS1-Fl2 Primer-F	
5’ 
TCGTCGGCAGCGTC
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*GAACCWGCGGARGGATCA*  (Primer sequence)
3’

##### ITS2 Primer-R	
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
```
### Step 4: DADA2 Denoising and Taxonomic Classification

Denoise and perform taxonomic classification using DADA2 and the UNITE reference classifier. The pre-trained classifier was obtained from the QIIME forum [^4].
[^4]:[pre-trained UNITE 9.0 classifiers for QIIME 2023.7 (and older!) - Community Contributions / Data resources - QIIME 2 Forum](https://forum.qiime2.org/t/pre-trained-unite-9-0-classifiers-for-qiime-2023-7-and-older/24140)
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

#### Filter Table:
```
qiime taxa filter-table \
--i-table table-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-table table-its-with-phyla-no-mitochondria-no-chloroplast.qza
```
#### Filter Sequences:
```
qiime taxa filter-seqs \
--i-sequences rep-seqs-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza
```

#### Unrooted Phylogenetic Tree Construction
In evolutionary biology, an unrooted phylogenetic tree is used to illustrate the relationships and evolutionary distances between entities such as species or sequences. Unlike their rooted counterparts, unrooted trees do not designate a common ancestor, but rather focus on depicting branching patterns and relative relationships. These trees aid in deducing genetic diversity and shared ancestry by scrutinising molecular data and creating diagrams that underscore evolutionary links. Unrooted trees play a crucial role in elucidating the evolutionary landscape and the interconnections among biological entities.

```
qiime phylogeny align-to-tree-mafft-fasttree \  
  --i-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza \  
  --output-dir phylogeny-align-to-tree-mafft-fasttree
```

The framework structure established with QIIME2 can be utilised for various QIIME2 plugins dedicated to a range of downstream analyses. In this scenario, we plan to use the Phyloseq package within the R programming environment for our subsequent analysis.