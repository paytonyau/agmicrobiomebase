
# UK CROP MICROBIOME CRYOBANK (ITS amplicon sequencing analysis)

To perform a sequencing analysis on the ITS amplicon data from the UK Crop Microbiome Cryobank, it is necessary to prepare a Linux environment with specific pre-installed packages. These packages, crucial for processing sequencing data, encompass a range of tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Prior to initiating the analysis, ensure that your system is equipped with all the necessary packages to facilitate a seamless and successful analysis. There are some key settings that are a bit different from 16s and thus we have a separate document for this reason.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [QIIME2 - the cummunity developed suite](https://qiime2.org/)
- [FIGARO](https://github.com/Zymo-Research/figaro)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)


### 1. Install packages/programs for amplicon sequencing data processing

#### Install Anaconda / Miniconda

To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda and create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

#### Install sequence quality check packages

For the installation of QIIME2, Anaconda or Miniconda can be utilised as they offer a self-contained environment and package manager. Begin by downloading and installing Miniconda, then create a new environment specifically for QIIME2. This approach simplifies the management of your QIIME2 installation and its dependencies.

Create a new environment:
```conda create --name qc```

Activate the conda environment
``` conda activate qc ```

#### Install fastqc & multiqc

```conda install -c bioconda fastqc multiqc```

#### Install Trimmomatic

```conda install -c conda-forge trimmomatic```
Noted that Java may need to be installed before running Trimmomatic.

#### Install FIGARO
FIGARO is a software that assists in estimating the truncation parameters for the QIIME2 DADA2 plugin. A pre-print detailing its usage is readily available[^1]. The detailed process for installing the software can be found on the author’s GitHub page (https://github.com/Zymo-Research/figaro). Further guidance is provided in John Quensen’s tutorial (https://john-quensen.com/tutorials/figaro/).

[^1]:Weinstein, M. et al. (2019) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. bioRxiv DOI: 10.1101/610394; Sasada, R. et al. (2020) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. J. Biomol. Tech. 31, S2


1.  `wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml`

If the provided URL is not available, you may also create a `figaro.yml` file on your own and include the information provided below, or use the `16s-step01-figaro.yml` file in our Github folder. 

```
name: figaro
channels:
  - bioconda
  - defaults
  - conda-forge
dependencies:
  - python >=3.6,<3.7.0
  - numpy
  - scipy==1.2.1
  - matplotlib==3.0.2
``` 

2.  `conda env create -n figaro -f figaro.yml`
    
3.  `git clone https://github.com/Zymo-Research/figaro.git`
    
4.  `unzip master.zip`
    
5.  `rm master.zip`
    
6.  `cd figaro-master/figaro`
    
7.  `chmod 755 *.py`

#### Install QIIME2

Run the workflow from the URL below provided in a specific conda environment, which ensures that the correct version of the Python required packages are being used for QIIME2.
<https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment>

### STEP 1
#### Examine all the fastq files using FastQC and consolidate the results into a single report with the help of MultiQC
This procedure is designed to visualise the quality of sequences and any remaining adaptors in each individual fastq file.

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

Due to the high variability in length of ITS, the trimming process used in 16s will not be adapted, but an estimation of error rates still needs to be used for the denoising process.

The outcome from the FIGARO run will generate four different files as listed below:

-   trimParameters.json
-   reverseExpectedError.png
-   forwardExpectedError.png
-   figaro.1672786258486002.log

#### Generate a manifest.txt for QIIME2 artifact (modification required to fit the QIIME2 format)

To create a QIIME2 artifact for the analysis, a `manifest` file for the corresponding fastq files are required. This requried a format entails the inclusion of “sample-id, absolute file-path, direction” and below is an example,

|sample-id| absolute file-path| direction|
| ----------- | ----------- |----------- |
|CO-CL-BO-1|/raw/Fastq/CO-CL-BO-1_S65_L001_R1_001.fastq.gz|forward|
|CO-CL-BO-1|/raw/Fastq/CO-CL-BO-1_S65_L001_R2_001.fastq.gz|reverse|
|CO-CL-BO-2|/raw/Fastq/CO-CL-BO-2_S76_L001_R1_001.fastq.gz|forward|
|CO-CL-BO-2|/raw/Fastq/CO-CL-BO-2_S76_L001_R2_001.fastq.gz|reverse|
|......|...|...|

**sample-id**: This is the identifier for each sample. In this case,  `CO-CL-BO-1`  and  `CO-CL-BO-2`  are the identifiers for two different samples.
    
**absolute file-path**: This is the full path to the location of the sequencing data files on your system. For example,  `/raw/Fastq/CO-CL-BO-1_S65_L001_R1_001.fastq.gz`  is the path to the forward read of the first sample.
    
**direction**: This indicates whether the sequencing data file is a forward read (`forward`) or a reverse read (`reverse`). In paired-end sequencing, you have both forward and reverse reads for each sample.


First, we will generate a `manifest.csv` file, which is an essential component in the process of creating a QIIME 2 artifact (.QZA) file. The manifest.csv file provides QIIME 2 with the necessary information about your sequencing reads. We also need to create a `manifest.txt` file that lists all directories in the current working directory. This can be done using the following command:

```ls -d "$PWD"/* > manifest.txt```

This command lists all directories (`ls -d`) in the current working directory (`"$PWD"/*`) and redirects the output to a file named  `manifest.txt`  (`> manifest.txt`).

Next, we will use a script to transform the  `manifest.txt`  file into the required  `manifest.csv`  format. You can copy the script below and paste it into a  `.sh`  file:

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
This script reads the `manifest.txt` file line by line, extracts the sample name and direction from each line, and writes this information into the `manifest.csv` file in the required format for QIIME2.

### Step 3: Import FASTQs as QIIME2 Artifact

Once we have identified that we have good quality ITS amplicon sequencing data, and defined the error rates, we can then move to the next step:

#### Activate the environment 
You can now activate this conda environment with this command:

```
conda activate qiime2-2023.5
```

#### Import FASTQs as QIIME 2 artifact
To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an “artifact” (with the extension QZA). The first step is to import the paired-end raw reads as a QZA file. Note that you need the manifest.csv file we mentioned previously and please following the QIIME2 requirement.

```
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv \
--output-path demux \
--input-format PairedEndFastqManifestPhred33
```
> Input file:
manifest.csv
the list of fastq (.fastq.gz) sequencing files

> Output file:
demux.qza

#### Summarise FASTQs
Obtain sequencing data quality information:
```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```
> Input file:
demux.qza

> Output file:
demux.qzv

#### Trim primers using `Cutadapt` from the QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the ITS1 region. DADA2's chimaera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimaeras to be identified. Please note that if you have identified a relatively higher adapter sequence contamination (which can be identified from the multiQC report), those adapter sequences should also be removed.

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

> Input file:
demux.qza

> Output file:
primer-trimmed-demux.qza

Please note that for the step here, we also trimmed adapter sequences, and we used the command `--p-discard-untrimmed` to discard the untrimmed reads that are highly variable in length in the ITS region. This was done to improve the accuracy. 

### Denoising the reads into amplicon sequence variants using DADA2

Denoise and perform taxonomic classification using DADA2 and the UNITE reference classifier. The pre-trained classifier was obtained from the QIIME forum [^4].
[^4]:[pre-trained UNITE 9.0 classifiers for QIIME 2023.7 (and older!) - Community Contributions / Data resources - QIIME 2 Forum](https://forum.qiime2.org/t/pre-trained-unite-9-0-classifiers-for-qiime-2023-7-and-older/24140)

Based on the projected outcomes of the FIGARO programme, the expected parameters - error rates (`--p-max-ee-f` & `--p-max-ee-r`) are 2. However, again, due to the highly  variable nature of ITS, we set 0 (no trimming position applied) for the `denoise-paired` process. 

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
> Input file:
primer-trimmed-demux.qza

> Output files:
rep-seqs-its.qza
table-its.qza
stats-its.qza

#### Assign taxonomy to ASVs by running taxonomic classification

The command for running the taxonomic classification is the most time-consuming and memory-intensive commands in the SOP. If you encounter an error due to insufficient memory and cannot increase your memory usage, consider adjusting the `--p-reads-per-batch` option to a value lower than the default (which is dynamic, depending on sample depth and the number of threads), and try running the command with fewer jobs (for example, set `--p-n-jobs` to `1`).  Instead of using Silva or Greengene 2 for 16s, for ITS we selected UNITE 9 as the reference database for the Taxonomic Classification.

```
qiime feature-classifier classify-sklearn \
--i-classifier /mnt/shared/scratch/pyau/qiime2-ref/unite-9-dynamic-s-all-29.11.2022-Q2-2023.5.qza \
--i-reads rep-seqs-its.qza \
  --p-n-jobs 1 \
--o-classification taxonomy-its.qza
```
> Input files:
unite-9-dynamic-s-all-29.11.2022-Q2-2023.5.qza
rep-seqs-its.qza

> Output file:
taxonomy-its.qza


#### Filtering out contaminant and unclassified ASVs

With taxonomy assigned to our ASVs, we can leverage this information to eliminate ASVs that are likely contaminants or noise, based on their taxonomic labels. Mitochondrial and chloroplast 16S sequences are common contaminants in 16S sequencing data and can be removed by excluding any ASV containing these terms in its taxonomic label. It may also be beneficial to exclude any ASV unclassified at the phylum level, as these sequences are more likely to be noise (e.g., potential chimeric sequences).

#### Filter Table:
```
qiime taxa filter-table \
--i-table table-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-table table-its-with-phyla-no-mitochondria-no-chloroplast.qza
```
> Input files:
table-its.qza
taxonomy-its.qza

> Output files:
table-its-with-phyla-no-mitochondria-no-chloroplast.qza


#### Filter Sequences:
```
qiime taxa filter-seqs \
--i-sequences rep-seqs-its.qza \
--i-taxonomy taxonomy-its.qza \
--p-include p__ \
--p-exclude mitochondria,chloroplast \
--o-filtered-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza
```
> Input files:
rep-seqs-its.qza
taxonomy-its.qza

> Output files:
rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza

#### Unrooted Phylogenetic Tree Construction
In evolutionary biology, an unrooted phylogenetic tree is used to illustrate the relationships and evolutionary distances between entities such as species or sequences. Unlike their rooted counterparts, unrooted trees do not designate a common ancestor, but rather focus on depicting branching patterns and relative relationships. These trees aid in deducing genetic diversity and shared ancestry by scrutinising molecular data and creating diagrams that underscore evolutionary links. Unrooted trees play a crucial role in elucidating the evolutionary landscape and the interconnections among biological entities.

```
qiime phylogeny align-to-tree-mafft-fasttree \  
  --i-sequences rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza \  
  --output-dir phylogeny-align-to-tree-mafft-fasttree
```

> Input file:
rep-seqs-its-with-phyla-no-mitochondria-no-chloroplast.qza

> Output files (under the subfolder phylogeny-align-to-tree-mafft-fasttree):
alignment.qza
masked_alignment.qza
rooted_tree.qza
tree.qza

The framework structure established with QIIME2 can be utilised for various QIIME2 plugins dedicated to a range of downstream analyses. In this scenario, we plan to use the Phyloseq package within the R programming environment for our subsequent analysis.