# UK CROP MICROBIOME CRYOBANK (Sequencing Analysis)

If you want to conduct sequencing analysis on the 16S amplicon data obtained from THE UK CROP MICROBIOME CRYOBANK (https://agmicrobiomebase.org), you'll need to have a Linux environment set up and certain packages pre-installed. These packages are essential for the sequencing data process and include various tools for quality control, sequence alignment, taxonomic classification, and diversity analysis. Before starting the analysis, make sure you have all the required packages installed on your system to ensure a smooth and successful analysis.

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [QIIME2 - the cummunity developed suite](https://qiime2.org/)
- [FIGARO](https://github.com/Zymo-Research/figaro)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://multiqc.info/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### 1. Install packages/programs for 16s amplicon sequencing data process

##### Install Anaconda / Miniconda

To install QIIME2, you can use Anaconda or Miniconda, which provide a self-contained environment and package manager. Download and install Miniconda and create a new environment for QIIME2. This allows you to manage your QIIME2 installation and dependencies easily.

##### Install sequence quality check packages

Create a new environment and install quality check packages such as FastQC, MultiQC, and Trimmomatic to ensure high-quality sequencing data. This avoids conflicts with other software and helps to streamline the quality control process.

```conda create --name qc```

```conda install -c bioconda fastqc multiqc```

##### Install Trimmomatic

```conda install -c conda-forge trimmomatic```
Noted that java may need to be pre-installed before the run.

##### Install FIGARO

FIGARO is a program that performs estimations when deciding to use the truncation parameters of the QIIME2 DADA2 plug-in. A pre-print of its usage is available.

```
git clone https://github.com/Zymo-Research/figaro.git
cd figaro
```

##### Install QIIME2

Run the workflow in a specific conda environment, which makes sure the correct version of the Python required packages are being used for QIIME2.
<https://docs.qiime2.org/2022.4/install/native/#install-qiime-2-within-a-conda-environment>

#### STEP 1
##### Screen all the fastq files by using FastQC and integrated into one single report by using multiQC


This process is to visulalie the sequence quality and the adaptors left in eash single individual fastq file. 
``` conda activate qc ```

``` fastqc *.fastq.gz ```
``` multiqc . ```

#### STEP 2
##### Sequencing reads trimming for FIGARO
Trim the sequence reads to a length of 250bp and discard reads shorter than 250bp to ensure consistency on each of the as required by the program.[^1] Cutadapt or Trimmomatic can be applied for the process. Here, we used Trimmomatic and please noted that the trimmed sequencing data are only for FIGARO.
[^1]: https://github.com/Zymo-Research/figaro/issues/37

```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:248 CROP: 248
```

##### Identify the optimal setting denosing in DADA2 by using FIGARO

Following the trimming process, the Figaro program should be utilised. The program involved applying various read lengths to estimate the error rates for DADA2. The aim was to optimize the accuracy of the DADA2 algorithm and improve its performance in accurately identifying and correcting sequencing errors in the dataset.

As the average length of **16s V3-V4** region is **460 nt** and minus the length of the forward and reverse primers nt, it has a 209 nt gap from both sides.
Forward: 460 - 251 = 209 nt
Reverse: 460 - 251 = 209 nt

Therefore, the overlapped nt for both side, which can be used for merging, would be 42 nt
460 - (209 + 209) = 42 nt

In the scope of this study, the 16S V3-V4 region was targeted, employing a forward primer of 17 nucleotides and a reverse primer spanning 21 nucleotides.
By primer removal,
Forward: 251 - 17 = 234 nt
Reverse: 251 - 21 = 230 nt

On DADA2, the minimum merging length is 12, if we set the mering overlap nt as 12,

42 nt - 12 nt = 30 nt

Forward: 234 - 30 = 204 nt
Reverse: 230 - 30 = 200 nt

Therefore, total minimum length
204 + 200 = 404 nt

Several iterations were undertaken to determine the optimal length for the 16S V3-V4 region in our dataset. The analysis revealed that the region's length is 426 nucleotides at 50% and 428 nucleotides at 98% similarity levels. 

The table below is based on the maximum length of 447 nt with the minimum of 12 nucleotides sequence overlap. 

|Percentile (%)|Plate 1*|Plate 2*|Plate 3*|ALL*|
| ----------- | ----------- |----------- |----------- |----------- |
|2| 402|	402|	402|	402|
|9| 402|	402|	402|	402|
|25|	404|	409|	405|	405|
|50|	426|	426|	426|	426|
|75|	427|	427|	427|	427|
|91|	428|	428|	428|	428|
|98|	428|	428|	428|	428|
*Length (nts)

Taking into account these findings along with the favorable sequence quality across the dataset, a maximum length of 428 nucleotides was selected for achieving resolution. This length ensures a 32-nucleotide overlap, contributing to the attainment of the highest resolution in the analysis. 

460 - 428 = 32 nt 

In other word, there is 42 - 32 = 10 nt from both forward and reverse sequence can be discarded.

Furthermore, adhering to the DADA2 recommendation for the v3-v4 16S region, it is advised to consider a minimum length of 20 nucleotides [^2]. Therefore, the FIGARO settings would be  
[^2]:https://benjjneb.github.io/dada2/tutorial_1_8.html

```python figaro.py -i (path) -o (path) -f 17 -r 21 -a 428 -m 20```

Here `-f` indecates the length of the forward primer - 17 bp, and `-r` indecates the length of the reverse primer - 21bp. `-a` indecates the length of the merged sequence, this required to exclude the length of the primers on both sides with the consideration of the total merging length. -m indecates the requirement of the merging base pair. 
Noted that the default DADA2 merging rate is 12bp (not 20bp).[^3] 
[^3]:https://forum.qiime2.org/t/extremely-low-merged-score-after-using-dada2/16085/3 and https://benjjneb.github.io/dada2/ReleaseNotes_1_8.html

Below is the outcomes from using FIGARO with the setting above
|Rank | T* position (F, R)| maxEE* (F, R)|Read Retention (%)|Score| 
| ----------- | ----------- |----------- |----------- |----------- |


*T, Trim; **EE, expected error


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

#### STEP 3
### Activate the envernment
You can activate this conda environment with this command (you may need to swap in source for conda if you get an error):

```
conda activate qiime2-2022.4
```

##### Import FASTQs as QIIME 2 artifact
To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an "artifact" (with the extension QZA). The first step is to import the paried-end raw reads as a QZA file. Noted that you need to prepare the manifast.csv file following the QIIME2 requirement.

```
qiime tools import
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33
```

##### Summarise FASTQs
We can also obtain the sequencing data quality information by executing the command below.

```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

##### Trim primers with using Cutadapt QIIME2 plugin
Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the 16s V3-V4 region. DADA2's chimera removal step requires primers to have been removed, as otherwise the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimeras to be identified.

The primes used in this dataset for 16s V3-V4 region
###### V3-V4 Primer-F	
5’ 
TCGTCGGCAGCGTC
==AGATGTGTATAAGAGACAG== (Adaptor sequence)
*CCTACGGGNGGCWGCAG*  (Primer sequence)
3’

###### V3-V4 Primer-R	
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

##### Summarise FASTQs after the primers trimming

You can run the demux summarize command after trimming the reads to get a report of the number of reads per sample and quality distribution across the reads. This generates a more basic output compared to FASTQC/MultiQC, but is sufficient for this step.

```
qiime demux summarize \
  --i-data primer-trimmed-demux.qza \
  --o-visualization primer-trimmed-demux.qzv
  ```

##### Denoising the reads into amplicon sequence variants using DADA2

According to the projected results of the FAGRIO programme, the anticipated parameters - error rates (--p-max-ee-f & --p-max-ee-r) and the optimal truncation length (--p-trunc-len-f & --p-trunc-len-r) - can be employed in this context to achieve optimal denoising results.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed-demux.qza \
  --p-trunc-len-f 222 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads 8 \
  --o-representative-sequences 422_222_220_rep-seqs.qza \
  --o-table 422_222_220_table.qza \
  --o-denoising-stats 422_222_220_stats.qza
```

Please noted that taxonomic classifiers exhibit optimal performance when they are tailored to your sample and sequencing processes. This encompasses factors such as the primers employed for amplification and the precise length of sequence reads utilised[^4]. It should be highlighted that employing full-length taxonomic classifiers remains a feasible approach. 
[^4]: https://docs.qiime2.org/2023.7/tutorials/feature-classifier/

##### Assign taxonomy to ASVs by running taxonomic classification

You can run the taxonomic classification with this command, which is one of the longest running and most memory-intensive command of the SOP. If you receive an error related to insufficient memory (and if you cannot increase your memory usage) then you can look into the --p-reads-per-batch option and set this to be lower than the default (which is dynamic depending on sample depth and the number of threads) and also try running the command with fewer jobs (e.g. set --p-n-jobs 1). In addition, you can also assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit learn Python library and the SILVA 138 or Greengene2 database.

```
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/shared/scratch/pyau/ref/silva-138-99-nb-classifier.qza \
  --i-reads 422_222_220_rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification 422_222_220_taxonomy_silva138.qza
```
##### visualise the taxonomy outcomes
```
qiime metadata tabulate \
  --m-input-file 422_222_220_taxonomy_silva138.qza \
  --o-visualization 422_222_220_taxonomy_silva138.qzv
```

##### Filter out contaminant and unclassified ASVs

Now that we have assigned taxonomy to our ASVs we can use that information to remove ASVs which are likely contaminants or noise based on the taxonomic labels. Two common contaminants in 16S sequencing data are mitochondrial and chloroplast 16S sequences, which can be removed by excluding any ASV which contains those terms in its taxonomic label. It can also be useful to exclude any ASV that is unclassified at the phylum level since these sequences are more likely to be noise (e.g. possible chimeric sequences). 
Note that if your data has not been classified against SILVA you will need to change P__ to be a string that enables phylum-level assignments to be identified or simply omit that line. You may also want to omit that line if you are studying a poorly characterised environment where you have a good chance of identifying novel phyla.

  ```
qiime taxa filter-table \
  --i-table 422_222_220_table.qza \
  --i-taxonomy 422_222_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table 422_222_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
  ```

```
qiime taxa filter-seqs \
  --i-sequences 422_222_220_rep-seqs.qza \
  --i-taxonomy 422_222_220_taxonomy_silva138.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences 422_222_220_rep-seqs_silva138-with-phyla-no-mitochondria-no-chloroplast.qza
  ```

### Unrooted phylogenetic tree construction
An unrooted phylogenetic tree is used in evolutionary biology to display the relationships and evolutionary distances among entities like species or sequences. Unlike rooted trees, unrooted trees do not specify a common ancestor, focusing instead on showing branching patterns and relative relationships. They help infer genetic diversity and shared ancestry by analyzing molecular data and generating diagrams that highlight evolutionary connections. Unrooted trees are essential for understanding the evolutionary landscape and interconnections among biological entities.
```
qiime phylogeny align-to-tree-mafft-fasttree \  
  --i-sequences rep-seqs-trimmed_231_228_2_2-with-phyla-no-mitochondria-no-chloroplast.qza \  
  --output-dir phylogeny-align-to-tree-mafft-fasttree
```

The qiime2 framework structure established here can applied for various qiime2 plugins dedicated to diverse downstream analyses. In this context, we intend to employ the phyloseq package within the R programming environment for our downstream analysis.
