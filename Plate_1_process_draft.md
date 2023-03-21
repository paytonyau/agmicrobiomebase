
# UK CROP MICROBIOME CRYOBANK

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

<https://docs.qiime2.org/2022.2/install/native/#install-qiime-2-within-a-conda-environment>

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
```
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.8 --file qiime2-2022.2-py38-linux-conda.yml
```

#### STEP 1
##### Screen all the fastq files by using FastQC and integrated into one single report by using multiQC


This process is to visulalie the sequence quality and the adaptors left in eash single individual fastq file. 
``` conda activate qc ```

``` fastqc *.fastq.gz ```
``` multiqc . ```

#### STEP 2
##### Sequencing reads trimming for FIGARO
Trim the sequence reads to a length of 250bp and discard reads shorter than 250bp to ensure consistency on each of the as required by the program. [^1] Cutadapt or Trimmomatic can be applied for the process. Here, we used Trimmomatic and please noted that the trimmed sequencing data are only for FIGARO.
[^1]: https://github.com/Zymo-Research/figaro/issues/37

```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
MINLEN:290 CROP: 290
```

##### Identify the optimal setting denosing in DADA2 by using FIGARO
```python figaro.py -i ~/test_figaro/ -o ~/test_figaro/ -f 17 -r 21 -a 460 -m 12```

Here ```-f``` indecates the length of the forward primer - 17 bp, and ```-r``` indecates the length of the reverse primer - 21bp. ```-a``` indecates the length of the merged sequence.
```-m``` indecates the requirement of the merging base pair. Noted that the default DADA2 merging in QIIME2 is 12bp (not 20bp).[^2]

[^2]:https://forum.qiime2.org/t/extremely-low-merged-score-after-using-dada2/16085/3
##### generate a manifest.txt for QIIME2 artifact (medification to fit for the QIIME2 required format)

```ls -d "$PWD"/* > manifest.txt```

#### STEP 3

##### activate the envernment

You can activate this conda environment with this command (you may need to swap in source for conda if you get an error):

```conda activate qiime2-2022.8```

##### Import FASTQs as QIIME 2 artifact

To standardise QIIME 2 analyses and to keep track of provenance (i.e. a list of what commands were previously run to produce a file) a special format is used for all QIIME2 input and output files called an "artifact" (with the extension QZA). The first step is to import the paried-end raw reads as a QZA file. Noted that you need to prepare the manifast.csv file following the QIIME2 requirement.

```
qiime tools import
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path demux \
  --input-format PairedEndFastqManifestPhred33
```

##### Summarize FASTQs
We can also obtain the sequencing data quality information by executing the command below.
```
qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv
```

##### Trim primers with cutadapt QIIME2 plugin

Screen out reads that do not begin with primer sequence and remove primer sequence from reads using the cutadapt QIIME2 plugin. The below primers correspond to the 16S V3-V4 region.
DADA2's chimera removal step requires primers to have been removed, as otherwise the ambiguous nucleotides in most primer sets cause large numbers of false-positive chimeras to be identified [^3].
[^3]: https://www.researchgate.net/post/Is_it_mandatory_to_remove_primers_and_barcodes_from_reads_before_starting_the_analyzes_in_qiime_2

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-cores 4 \
  --p-front-f CCTACGGGNGGCWGCAG \
  --p-front-r GACTACHVGGGTATCTAATCC \
  --o-trimmed-sequences primer-trimmed-demux.qza \
  --verbose \
  &> primer_trimming.log
  ```

##### Summarize FASTQs

You can run the demux summarize command after trimming the reads to get a report of the number of reads per sample and quality distribution across the reads. This generates a more basic output compared to FASTQC/MultiQC, but is sufficient for this step.

```
qiime demux summarize \
  --i-data primer-trimmed-demux.qza \
  --o-visualization primer-trimmed-demux.qzv
  ```

##### Denoising the reads into amplicon sequence variants using DADA2

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs primer-trimmed-demux.qza \
  --p-n-threads 8 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 247  \
  --p-trunc-len-r 225 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --o-representative-sequences rep-seqs-trimmed_247_225_2_2.qza \
  --o-table table-trimmed_247_225_2_2.qza \
  --o-denoising-stats stats-trimmed_247_225_2_2.qza
```

##### Assign taxonomy to ASVs by running taxonomic classification

You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit learn Python library and the SILVA database.

You can run the taxonomic classification with this command, which is one of the longest running and most memory-intensive command of the SOP. If you receive an error related to insufficient memory (and if you cannot increase your memory usage) then you can look into the --p-reads-per-batch option and set this to be lower than the default (which is dynamic depending on sample depth and the number of threads) and also try running the command with fewer jobs (e.g. set --p-n-jobs 1).

```
qiime feature-classifier classify-sklearn \
  --p-n-jobs 1 \
  --i-classifier /mnt/shared/scratch/pyau/qiime2-ref/silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs-trimmed_247_225_2_2.qza \
  --o-classification taxonomy-trimmed_247_225_2_2.qza
```

##### Filter out contaminant and unclassified ASVs

Now that we have assigned taxonomy to our ASVs we can use that information to remove ASVs which are likely contaminants or noise based on the taxonomic labels. Two common contaminants in 16S sequencing data are mitochondrial and chloroplast 16S sequences, which can be removed by excluding any ASV which contains those terms in its taxonomic label. It can also be useful to exclude any ASV that is unclassified at the phylum level since these sequences are more likely to be noise (e.g. possible chimeric sequences). Note that if your data has not been classified against SILVA you will need to change P__ to be a string that enables phylum-level assignments to be identified or simply omit that line. You may also want to omit that line if you are studying a poorly characterized environment where you have a good chance of identifying novel phyla.

```
qiime taxa filter-seqs \
  --i-sequences rep-seqs-trimmed_247_225_2_2.qza \
  --i-taxonomy taxonomy-trimmed_247_225_2_2.qza \
  --p-include p__ \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences rep-seqs-trimmed_247_225_5_4-with-phyla-no-mitochondria-no-chloroplast.qza
  ```




##### Other resources

There are many other possible QIIME 2 analyses that we recommend you look into. You may also find these other resources useful:

ALDEx2 for testing for differential relative abundances in R

corncob for testing for differential relative abundances in R

Building Random Forest Models in R

PICRUSt2 for metagenome prediction based on amplicon sequences

STAMP for straight-forward visualization (see workflow here)
