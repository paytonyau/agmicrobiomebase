

# UK CROP MICROBIOME CRYOBANK (16S amplicon sequencing analysis)

To perform a sequencing analysis on the 16S amplicon data from the UK Crop Microbiome Cryobank (https://agmicrobiomebase.org), it is necessary to prepare a Linux environment with specific pre-installed packages. These packages, crucial for processing sequencing data, encompass a range of tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Prior to initiating the analysis, ensure that your system is equipped with all the necessary packages to facilitate a seamless and successful analysis.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [QIIME2 - the cummunity developed suite](https://qiime2.org/)
- [FIGARO](https://github.com/Zymo-Research/figaro)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### 1. Install packages/programs for 16s amplicon sequencing data process

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

#### Identify the optimal setting denosing in DADA2 by using FIGARO

Following the trimming process, the Figaro program should be utilised. The program involved applying various read lengths to estimate the error rates for DADA2. The aim was to optimize the accuracy of the DADA2 algorithm and improve its performance in accurately identifying and correcting sequencing errors in the dataset.

As the average length of **16s V3-V4** region is **460 nt** and minus the length of the maximum length of the reads, which is 251 nt, it has a gap of ~209 nt 
**Forward: 460 - 251 = 209 nt**
**Reverse: 460 - 251 = 209 nt**

Therefore, the **overlapped** nt for both sides, which can be used for merging, would be
**460 - (209 + 209) = 42 nt** ......(1)

In DADA2, the **minimum merging length or default setting** is **12 nt** [^3] . Ideally, if we set the merging overlap nucleotides as 12, there is room (*X*) for 30 nucleotides that can be trimmed on both sides.
**42 nt -** *X* **nt = 12 nt**... (2)
*X* **= 30 nt**
[^3]:https://benjjneb.github.io/dada2/ReleaseNotes_1_8.html

The forward primer, which is used for the 16S V3-V4 region and is 17 nucleotides long, and the reverse primer, which spans 21 nucleotides, both directions need to be removed to minimise the errors induced by the primers. Therefore the length, **after the primer removal**, on average is **460 - (17 + 21) = 422 nt**...(3)

However, to maintain a longer length for more biological identities and considering the low quality of both ends of the sequencing reads, an optimal trimming setting needs to be taken into account. 

Several iterations were undertaken to determine the optimal length for the 16S V3-V4 region in our dataset using the length of XXX setting for the initial analysis. The analysis revealed that at a level of 50%, the length of the region is 424 nt, and at a  level of 98%, the length is 428 nt. We took 428 nt as the target in order to maintain the most possible biological information and at the same time to trim down the low-quality reads.
**428 - 422 = 6 nt**...(4)
 
 In another word, we added 6 nt for the merging.
 From (2) and (4), we have 
 **12 + 6 = 18 nt**

Moreover, according to the DADA2 tutorial, it is suggested to set the merging level to **20 nt** for the v3-v4 16S region [^4]. 
[^4]:https://benjjneb.github.io/dada2/tutorial_1_8.html & https://forum.qiime2.org/t/extremely-low-merged-score-after-using-dada2/16085/3

Therefore, the FIGARO setting would be  
```python figaro.py -i (path) -o (path) -f 17 -r 21 -a 428 -m 20```

Here `-f` indicates the length of the forward primer - **`17` nt**, and `-r` indicates the length of the reverse primer - **`21` nt**. `-a` indicates the length of the merged sequence - **`428` nt**, this required to exclude the length of the primers on both sides with the consideration of the total merging length. `-m` indicates the requirement of the merging base pair - **`20` nt**. 


Below is the outcomes from using FIGARO with the setting above
|Rank | T* position (F, R)| maxEE* (F, R)|Read Retention (%)|Score| 
| ----------- | ----------- |----------- |----------- |----------- |
|**1**|**245, 241**|**2, 2**|**80.49**|**78.4904**|
|**2**|**244, 242**|**2, 2**|**80.49**|**78.4870**|
|**3**|**246, 240**|**2, 2**|**80.47**|**78.4743**|
|**4**|**247, 239**|**2, 2**|**80.46**|**78.4564**|
|**5**|**240, 246**|**2, 2**|**80.44**|**78.4403**|
*T, Trim; **EE, expected error


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

#### STEP 3
#### Activate the environment 
You can activate this conda environment with this command (you may need to swap in source for conda if you get an error):

```
conda activate qiime2-2023.5
```

#### Import FASTQs as QIIME 2 artifact
To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an "artefact" (with the extension QZA). The first step is to import the paired-end raw reads as a QZA file. Noted that you need to prepare the manifast.csv file following the QIIME2 requirement.

```
qiime tools import
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33
```

#### Summarise FASTQs
We can also obtain the sequencing data quality information by executing the command below.

```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

#### Trim primers using Cutadapt QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the 16s V3-V4 region. DADA2's chimaera removal step requires primers to have been removed, as otherwise, the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimaeras to be identified.

The primes used in this dataset for 16s V3-V4 region
##### V3-V4 Primer-F	
5’ 
TCGTCGGCAGCGTC
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*CCTACGGGNGGCWGCAG*  (Primer sequence)
3’

##### V3-V4 Primer-R	
5’ 
GTCTCGTGGGCTCGG
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*GACTACHVGGGTATCTAATCC* (Primer sequence)
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

#### Summarise FASTQs after the primers trimming

You can run the demux summarize command after trimming the reads to get a report of the number of reads per sample and quality distribution across the reads. This generates a more basic output compared to FASTQC/MultiQC, but is sufficient for this step.

```
qiime demux summarize \
  --i-data primer-trimmed-demux.qza \
  --o-visualization primer-trimmed-demux.qzv
  ```

#### Denoising the reads into amplicon sequence variants using DADA2

Based on the projected outcomes of the FIGARO programme, the expected parameters - error rates (`--p-max-ee-f` & `--p-max-ee-r`) and the optimal truncation length (`--p-trunc-len-f` & `--p-trunc-len-r`) - can be utilised in this scenario to achieve the best denoising results.
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

Please note that taxonomic classifiers exhibit optimal performance when they are tailored to your sample and sequencing processes. This encompasses factors such as the primers employed for amplification and the precise length of sequence reads utilised[^5]. It should be highlighted that employing full-length taxonomic classifiers remains a feasible approach. 
[^5]: https://docs.qiime2.org/2023.7/tutorials/feature-classifier/

#### Assign taxonomy to ASVs by running taxonomic classification

The command for running the taxonomic classification is one of the most time-consuming and memory-intensive commands in the SOP. If you encounter an error due to insufficient memory and cannot increase your memory usage, consider adjusting the `--p-reads-per-batch` option to a value lower than the default (which is dynamic, depending on sample depth and the number of threads), and try running the command with fewer jobs (for example, set `--p-n-jobs` to `1`). Additionally, you can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit-learn Python library, along with the SILVA 138 or Greengene2 database. The pre-trained classifiers can be obtained from Qiime2 data-resources [^6]. 
[^6]: https://docs.qiime2.org/2023.7/data-resources/
```
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/ref/silva-138-99-nb-classifier.qza \
  --i-reads 428_228_220_rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification 428_228_220_taxonomy_silva138.qza
```
#### visualise the taxonomy outcomes
```
qiime metadata tabulate \
  --m-input-file 428_228_220_taxonomy_silva138.qza \
  --o-visualization 428_228_220_taxonomy_silva138.qzv
```

#### Filter out contaminant and unclassified ASVs

With taxonomy assigned to our ASVs, we can leverage this information to eliminate ASVs that are likely contaminants or noise, based on their taxonomic labels. Mitochondrial and chloroplast 16S sequences are common contaminants in 16S sequencing data and can be removed by excluding any ASV containing these terms in its taxonomic label. It may also be beneficial to exclude any ASV unclassified at the phylum level, as these sequences are more likely to be noise (e.g., potential chimeric sequences).

Please note that if your data has not been classified against SILVA, you will need to modify 'P__' to a string that allows phylum-level assignments to be identified, or simply omit that line. If you are studying a poorly characterized environment where there is a good chance of identifying novel phyla, you might also want to omit that line.
```
qiime taxa filter-table \
  --i-table 428_228_220_table.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
```

```
qiime taxa filter-seqs \
  --i-sequences 428_228_220_rep-seqs.qza \
  --i-taxonomy 428_228_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
  ```

#### Visualise the rep seq
```
qiime metadata tabulate \
  --m-input-file 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza \
  --o-visualization 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qzv
```
### Unrooted phylogenetic tree construction
In evolutionary biology, an unrooted phylogenetic tree is employed to illustrate the relationships and evolutionary distances between entities such as species or sequences. Unlike their rooted counterparts, unrooted trees do not designate a common ancestor, but rather concentrate on depicting branching patterns and relative relationships. These trees aid in deducing genetic diversity and shared ancestry by scrutinising molecular data and creating diagrams that underscore evolutionary links. Unrooted trees play a crucial role in elucidating the evolutionary landscape and the interconnections among biological entities.
```
qiime phylogeny align-to-tree-mafft-fasttree \  
  --i-sequences 428_228_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza \  
  --output-dir phylogeny-align-to-tree-mafft-fasttree
```

The framework structure established with QIIME2 can be utilized for various QIIME2 plugins dedicated to a range of downstream analyses. In this scenario, we plan to use the Phyloseq package within the R programming environment for our subsequent analysis.
