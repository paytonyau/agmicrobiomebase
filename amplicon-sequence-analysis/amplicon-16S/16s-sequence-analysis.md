# UK CROP MICROBIOME CRYOBANK (16S amplicon sequencing analysis)

To perform a sequencing analysis on the 16S amplicon data from the UK Crop Microbiome Cryobank (https://agmicrobiomebase.org), you need to prepare a Linux environment with specific pre-installed packages. These packages are crucial for processing sequencing data and include a range of tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Before initiating the analysis, ensure that your system is equipped with all the necessary packages to facilitate a seamless and successful analysis.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [QIIME2 - the cummunity developed suite](https://qiime2.org/)
- [FIGARO](https://github.com/Zymo-Research/figaro)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [PICRUSt2](https://github.com/picrust/picrust2)

## 1. Install packages/programs for amplicon sequencing data processing

### Install Anaconda / Miniconda

To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda and create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

### Install sequence quality check packages

For the installation of QIIME2, Anaconda or Miniconda can be utilised as they offer a self-contained environment and package manager. Begin by downloading and installing Miniconda, then create a new environment specifically for QIIME2. This approach simplifies the management of your QIIME2 installation and its dependencies.

Create a new environment:

```conda create --name qc```

Activate the conda environment

``` conda activate qc ```

### Install fastqc & multiqc

```conda install -c bioconda fastqc multiqc```

### Install Trimmomatic

```conda install -c conda-forge trimmomatic```

Noted that Java may need to be installed before running Trimmomatic.

### Install FIGARO
FIGARO is a software that assists in estimating the truncation parameters for the QIIME2 DADA2 plugin. A pre-print detailing its usage is readily available[^1]. The detailed process for installing the software can be found on the author’s GitHub page (https://github.com/Zymo-Research/figaro). Further guidance is provided in John Quensen’s tutorial (https://john-quensen.com/tutorials/figaro/).

[^1]:Weinstein, M. et al. (2019) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. bioRxiv DOI: 10.1101/610394; Sasada, R. et al. (2020) FIGARO: An efficient and objective tool for optimizing microbiome rRNA gene trimming parameters. J. Biomol. Tech. 31, S2


1.  `wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml` or `wget https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-step01-figaro.yml`

2.  `conda env create -n figaro -f 16s-step01-figaro.yml`
    
3.  `git clone https://github.com/Zymo-Research/figaro.git`
    
4.  `unzip master.zip`
    
5.  `rm master.zip`
    
6.  `cd figaro-master/figaro`
    
7.  `chmod 755 *.py`

### Install QIIME2

Run the workflow from the URL below provided in a specific conda environment, which ensures that the correct version of the Python required packages are being used for QIIME2.
<https://docs.qiime2.org/2023.5/install/native/#install-qiime-2-within-a-conda-environment>

## STEP 1
#### Examine all the fastq files using FastQC and consolidate the results into a single report with the help of MultiQC
This procedure is designed to visualise the quality of sequences and any remaining adaptors in each individual fastq file.

``` fastqc *.fastq.gz ```

``` multiqc . ```

## STEP 2
#### Sequencing reads trimming for FIGARO

To maintain consistency as required by the program, sequence reads should be trimmed to a length of 248 nt, and any reads shorter than 248 nt should be discarded[^2]. This process can be accomplished using tools such as Cutadapt or Trimmomatic. In this case, Trimmomatic was used. 
[^2]: https://github.com/Zymo-Research/figaro/issues/37

```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248
```

### Identify the optimal setting denosing in DADA2 by using FIGARO

Following the trimming process, the FIGARO program should be utilised. The program involves applying various read lengths to estimate the error rates for DADA2. The aim is to optimise the accuracy of the DADA2 algorithm and improve its performance in accurately identifying and correcting sequencing errors in the dataset.

As the average length of **16s V3-V4** region is **460 nt** and minus the length of the maximum length of the reads, which is 251 nt, it has a gap of ~209 nt 
**Forward: 460 - 251 = 209 nt**
**Reverse: 460 - 251 = 209 nt**

Therefore, the **overlapped** nt for both sides, which can be used for merging, would be
**460 - (209 + 209) = 42 nt** ......(1)

In DADA2, the **minimum merging length or default setting** is **12 nt** [^3]. Ideally, if we can set the merging overlap nucleotides as 12, there is room (*X*) for 30 nucleotides that can be trimmed on both sides.
**42 nt -** *X* **nt = 12 nt**... (2)
*X* **= 30 nt**
[^3]:https://benjjneb.github.io/dada2/ReleaseNotes_1_8.html

The forward primer, which is used for the 16S V3-V4 region and is 17 nucleotides long, and the reverse primer, which spans 21 nucleotides, both directions need to be removed to minimise the errors induced by the primers. Therefore the length, **after the primer removal**, on average is **460 - (17 + 21) = 422 nt**...(3)

However, to maintain a longer length for more biological identities and considering the low quality of both ends of the sequencing reads, an optimal trimming setting needs to be taken into account.

Several iterations were undertaken to determine the optimal length for the 16S V3-V4 region in our dataset using the length of 38 nt overlapped setting for the initial analysis. We revealed that at a level of 50%, the length of the region is 424 nt, and at a level of 98%, the length is 428 nt. We took 428 nt as the target in order to maintain the most possible biological information and at the same time to trim down the low-quality reads. 

**428 - 422 = 6 nt**...(4)
 
 In other words, we added 6 nt for the merging. 
 From (2) and (4), we have 
 **12 + 6 = 18 nt**

Therefore, in the FIGARO setting we use **428 nt** as the longest biologically meaningful read and the setting would be  
```python figaro.py -i (path) -o (path) -f 17 -r 21 -a 428 -m 18```

Here `-f` indicates the length of the forward primer - **`17` nt**, and `-r` indicates the length of the reverse primer - **`21` nt**. `-a` indicates the length of the merged sequence - **`428` nt**, this required to exclude the length of the primers on both sides with the consideration of the total merging length. `-m` indicates the requirement of the merging base pair - **`18` nt**. 

> The outcome from the FIGARO run will generate four different files as listed below:
> 
> - trimParameters.json
> - reverseExpectedError.png
> - forwardExpectedError.png
> - figaro.1672786258486002.log


The top 5 outcomes from the FIGARO outcomes are as follows:

|Rank | T* position (F, R)| maxEE* (F, R)|Read Retention (%)|Score| 
| ----------- | ----------- |----------- |----------- |----------- |
|**1**|**245, 241**|**2, 2**|**80.49**|**78.4904**|
|**2**|**244, 242**|**2, 2**|**80.49**|**78.4870**|
|**3**|**246, 240**|**2, 2**|**80.47**|**78.4743**|
|**4**|**247, 239**|**2, 2**|**80.46**|**78.4564**|
|**5**|**240, 246**|**2, 2**|**80.44**|**78.4403**|
*T, Trim; **EE, expected error

The settings from the first outcome will be used for the `denoise` process.
-   **Trim Position (F, R):**  245, 241
-   **Max Expected Error (F, R):**  2, 2

Please note that the trimmed sequences from FIGARO will not be used for any other analysis. The original fastq files should be used to create the QIIME2 artifact.

### Generate a manifest.txt for QIIME2 artifact 

To create a QIIME2 artifact for the analysis, a `manifest` file for the corresponding fastq files are required. This requried a format entails the inclusion of “sample-id, absolute file-path, direction” and below is an example,

|sample-id| absolute file-path| direction|
| ----------- | ----------- |----------- |
|CO-CL-BO-1|/raw/Fastq/CO-CL-BO-1_S65_L001_R1_001.fastq.gz|forward|
|CO-CL-BO-1|/raw/Fastq/CO-CL-BO-1_S65_L001_R2_001.fastq.gz|reverse|
|CO-CL-BO-2|/raw/Fastq/CO-CL-BO-2_S76_L001_R1_001.fastq.gz|forward|
|CO-CL-BO-2|/raw/Fastq/CO-CL-BO-2_S76_L001_R2_001.fastq.gz|reverse|
|......|...|...|

-   **sample-id**: This is the identifier for each sample. In this case,  `CO-CL-BO-1`  and  `CO-CL-BO-2`  are the identifiers for two different samples.
    
-   **absolute file-path**: This is the full path to the location of the sequencing data files on your system. For example,  `/raw/Fastq/CO-CL-BO-1_S65_L001_R1_001.fastq.gz`  is the path to the forward read of the first sample.
    
-   **direction**: This indicates whether the sequencing data file is a forward read (`forward`) or a reverse read (`reverse`). In paired-end sequencing, you have both forward and reverse reads for each sample.


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

## STEP 3

Once we have identified that we have good quality 16s amplicon sequencing data, and defined the optimal trimming positions and error rates, we can then move to the next step:

### Activate the environment 
You can now activate this conda environment with this command:

```
conda activate qiime2-2023.5
```

### Import FASTQs as QIIME 2 artifact
To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an “artifact” (with the extension QZA). The first step is to import the paired-end raw reads as a QZA file. Note that you need the manifest.csv file we mentioned previously and please following the QIIME2 requirement.

```
qiime tools import
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33
```
> Input files:
> - manifest.csv
> -  the list of fastq (.fastq.gz) sequencing files
> 
> Output file:
> - demux.qza

### Summarise FASTQs
We can also obtain the sequencing data quality information by executing the command below.
```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

> Input file:
> - demux.qza
>
> Output file:
> - demux.qzv

The `demux.qzv` file can be upload to [QIIME 2 View](https://view.qiime2.org/) for an interactive visualisation

### Trim primers using Cutadapt QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the `cutadapt` QIIME2 plugin. The below primers correspond to the 16s V3-V4 region. DADA2’s chimera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimeras to be identified. Please note that if you have identified a relatively higher adapter sequence contamination (which can be identified from the multiQC report), those adapter sequences should also be removed.

The primes used in this dataset for 16s V3-V4 region
##### V3-V4 Primer-F	
5’ 
TCGTCGGCAGCGTC

*AGATGTGTATAAGAGACAG* (Adaptor sequence)

**CCTACGGGNGGCWGCAG**  (Primer sequence)
3’

##### V3-V4 Primer-R	
5’ 
GTCTCGTGGGCTCGG

*AGATGTGTATAAGAGACAG* (Adaptor sequence)

**GACTACHVGGGTATCTAATCC** (Primer sequence)
3’ 

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-cores 8 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --p-discard-untrimmed \
  --o-trimmed-sequences primer-trimmed-demux.qza \
  --verbose \
  &> primer_trimming.log
  ```

Input file:
- demux.qza

Output file:
- primer-trimmed-demux.qza 

Here is the elements and the explainations that we used for the trimming step.
-   `qiime cutadapt trim-paired`: This is the command to run the  `trim-paired`  method from the  `cutadapt`  plugin in QIIME 2. This method trims adapter sequences from paired-end sequence data.
-   `--i-demultiplexed-sequences demux.qza`: This specifies the input file of demultiplexed sequences (`demux.qza`) that you want to trim.
-   `--p-cores 8`: This sets the number of CPU cores to use in the trimming process to 8.
-   `--p-front-f CCTACGGGNGGCWGCAG`  and  `--p-front-r GACTACHVGGGTATCTAATCC`: These specify the adapter sequences to trim from the forward (`--p-front-f`) and reverse (`--p-front-r`) reads.
-   `--p-discard-untrimmed`: This option discards reads in which the adapter sequence could not be found.
-   `--o-trimmed-sequences primer-trimmed-demux.qza`: This specifies the output file (`primer-trimmed-demux.qza`) to which the trimmed sequences will be written.
-   `--verbose`: This option enables verbose output, which means the program will output more information about its progress.
-   `&> primer_trimming.log`: This redirects both the standard output and standard error to a log file (`primer_trimming.log`), so you can review it later for any potential issues or for record-keeping.

### Summarise FASTQs after the primers trimming

You can run the `demux summarize` command after trimming the reads to get a report of the number of reads per sample and quality distribution across the reads. This generates a more basic output compared to FASTQC/MultiQC, but is sufficient for this step.

```
qiime demux summarize \
  --i-data primer-trimmed-demux.qza \
  --o-visualization primer-trimmed-demux.qzv
  ```

> Input file:
> - primer-trimmed-demux.qza
>
> Output file:
> - primer-trimmed-demux.qzv

### Denoising the reads into amplicon sequence variants using DADA2

Based on the projected outcomes of the FIGARO programme, the expected parameters - error rates (`--p-max-ee-f` & `--p-max-ee-r`) and the optimal truncation length (`--p-trunc-len-f` & `--p-trunc-len-r`) - can be utilised in this scenario to achieve the best denoising results.
From Figaro, we have the trimming position of Forward = **245** and Reverse  = **241**, minus the length of the primers (as we previously trimmed). Therefore, 

Forward: **245 - 17 = 228 nt**
Reverse: **241 - 21 = 220 nt**

```
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
```

> Input file:
> - primer-trimmed-demux.qza
> 
> Output files:
> - [428_228_220_rep-seqs.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_rep-seqs.qza)
> - [428_228_220_table.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_table.qza)
> - [428_228_220_stats.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_stats.qza)

> ###### - 428: Represents the longest combined biologically meaningful read, denoted as 428nt.
> ###### - 228: The trimming position for the forward reads, denoted as 228nt.
> ###### - 220: The trimming position for the reverse reads, denoted as 220nt.

Here is the elements and the explainations that we used for the `denose-paired` step:
-   `qiime dada2 denoise-paired`: This is the command to run the  `denoise-paired`  method from the  `dada2`  plugin in QIIME 2. This method denoises paired-end sequence data.  
-   `--i-demultiplexed-seqs primer-trimmed-demux.qza`: This specifies the input file of demultiplexed sequences (`primer-trimmed-demux.qza`) that you want to denoise.  
-   `--p-trunc-len-f 228`  and  `--p-trunc-len-r 220`: These specify the positions at which to truncate the forward (`--p-trunc-len-f`) and reverse (`--p-trunc-len-r`) reads.   
-   `--p-trim-left-f 0`  and  `--p-trim-left-r 0`: These specify the positions at which to trim the forward (`--p-trim-left-f`) and reverse (`--p-trim-left-r`) reads.
-   `--p-max-ee-f 2`  and  `--p-max-ee-r 2`: These specify the maximum expected error rates for the forward (`--p-max-ee-f`) and reverse (`--p-max-ee-r`) reads.
-   `--p-n-threads 8`: This sets the number of CPU threads to use in the denoising process to 8.
-   `--o-representative-sequences 428_228_220_rep-seqs.qza`,  `--o-table 428_228_220_table.qza`, and  `--o-denoising-stats 428_228_220_stats.qza`: These specify the output files to which the representative sequences, feature table, and denoising stats will be written, respectively.

Please note that taxonomic classifiers exhibit optimal performance when they are tailored to your sample and sequencing processes. This encompasses factors such as the primers employed for amplification and the precise length of sequence reads utilised[^4]. It should be highlighted that employing full-length taxonomic classifiers remains a feasible approach.

[^4]: https://docs.qiime2.org/2023.5/tutorials/feature-classifier/

### Assign taxonomy to ASVs by running taxonomic classification

The command for running the taxonomic classification is the most time-consuming and memory-intensive commands in the SOP. If you encounter an error due to insufficient memory and cannot increase your memory usage, consider adjusting the `--p-reads-per-batch` option to a value lower than the default (which is dynamic, depending on sample depth and the number of threads), and try running the command with fewer jobs (for example, set `--p-n-jobs` to `1`). Additionally, you can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit-learn Python library, along with the SILVA 138 or Greengene2 database.The pre-trained classifiers can be obtained from Qiime2 data-resources [^5].  Here, we used Silva v.138 as the reference database for the process.

[^5]: https://docs.qiime2.org/2023.7/data-resources/
```
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads 428_228_220_rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification 428_228_220_taxonomy_silva138.qza
```

> Input files:
> - [silva-138-99-nb-classifier.qza](https://data.qiime2.org/2024.2/common/silva-138-99-nb-classifier.qza)
> - [428_228_220_rep-seqs.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_rep-seqs.qza)
> 
> Output file:
> - [428_228_220_taxonomy_silva138.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_taxonomy_silva138.qza)

Here is the elements and the explainations that we used for the `qiime feature-classifier classify-sklearn` step:

-   `qiime feature-classifier classify-sklearn`: This is the command to run the  `classify-sklearn`  method from the  `feature-classifier`  plugin in QIIME 2. This method classifies features using a specific classifier.
    
-   `--i-classifier silva-138-99-nb-classifier.qza`: This specifies the input file of the classifier (`silva-138-99-nb-classifier.qza`) that you want to use for classification.
    
-   `--i-reads 428_228_220_rep-seqs.qza`: This specifies the input file of representative sequences (`428_228_220_rep-seqs.qza`) that you want to classify.
    
-   `--p-n-jobs 8`: This sets the number of CPU jobs to use in the classification process to 8.
    
-   `--o-classification 428_228_220_taxonomy_silva138.qza`: This specifies the output file (`428_228_220_taxonomy_silva138.qza`) to which the classification results will be written.

### Visualising the taxonomy outcomes
```
qiime metadata tabulate \
  --m-input-file 428_228_220_taxonomy_silva138.qza \
  --o-visualization 428_228_220_taxonomy_silva138.qzv
```

> Input file:
> - [428_228_220_taxonomy_silva138.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_taxonomy_silva138.qza)

> Output file:
> - 428_228_220_taxonomy_silva138.qzv

### Filtering out contaminant and unclassified ASVs

With taxonomy assigned to our ASVs, we can leverage this information to eliminate ASVs that are likely contaminants or noise, based on their taxonomic labels. Mitochondrial and chloroplast 16S sequences are common contaminants in 16S sequencing data and can be removed by excluding any ASV containing these terms in its taxonomic label. It may also be beneficial to exclude any ASV unclassified at the phylum level, as these sequences are more likely to be noise (e.g., potential chimeric sequences).

Please note that if your data has not been classified against SILVA, you will need to modify ‘P__’ to a string that allows phylum-level assignments to be identified, or simply omit that line. If you are studying a poorly characterised environment where there is a good chance of identifying novel phyla, you might also want to omit that line.

The step below is to filter the table (`428_228_220_table.qza`) 

```
qiime taxa filter-table \
  --i-table 428_228_220_table.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
```
> Input files:
> - [428_228_220_table.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_table.qza)
> - [428_228_220_taxonomy_silva138.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_taxonomy_silva138.qza)
>
> Output files:
> - [428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza)

```
qiime taxa filter-seqs \
  --i-sequences 428_228_220_rep-seqs.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
  ```

> Input files:
> - [428_228_220_rep-seqs.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DDADA2_outputfiles/428_228_220_rep-seqs.qza)
> - [428_228_220_taxonomy_silva138.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_taxonomy_silva138.qza)
> 
> Output files:
> - [428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza)

### Visualise the Representative Sequence
```
qiime metadata tabulate \
  --m-input-file 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza \
  --o-visualization 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qzv
```
> Input files:
> - [428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/%5BQiime2%5DSilva_138_outputfiles/428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza)
> 
> Output files:
> - 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qzv


## Functional Prediction - `PICRUSt2`
`PICRUSt2` is a software tool for predicting functional abundances based solely on marker gene sequences. Please refer to the `PICRUSt2` wiki for comprehensive documentation and tutorials. Remember, “function” typically refers to gene families such as KEGG orthologs and Enzyme Classification numbers, but predictions can be made for a variety of other entities. There is also a tutorial from YouTube (https://www.youtube.com/watch?v=xZ7yc-GKcSk) that could be useful to follow.

**A.** Output `BIOMV210DirFmt` (BIOM) file

`qiime tools export --input-path 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza  --output-path table-exported`


**B.** Output `DNASequencesDirectoryFormat` (FASTA) format file

`qiime tools export --input-path 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza --output-path sequences`

Establish a fresh environment with PICRUSt2 implemented and initiate this environment. Should you encounter an error indicating the absence of default reference files when attempting to execute the software, you may need to consider installing from the source as an alternative.
https://github.com/picrust/picrust2/wiki/Installation

**C.** Create and activate the environment for `PICRUSt2`
Use the following command to create a new conda environment named `picrust2`:
`conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2`

**D.** Activate the new conda environment
`conda activate picrust2`

**E.** Run the program
`picrust2_pipeline.py -s sequences/dna-sequences.fasta -i table-exported/feature-table.biom -o picrust2_out_pipeline -p 1`