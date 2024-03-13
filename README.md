  

# UK Crop Microbiome Cryobank

## Agmicrobiomebase

  

This repository contains the scripts to process the fastq sequence files from the 16S and ITS amplicon Illumina sequencing and for the metagenomic Illumina sequencing for the UK Crop Microbiome Cryobank Project.

  

The UK Crop Microbiome Cryobank integrates genomic (DNA) data with a cryobank collection of samples for the soil microbiomes of the UK major crop plant systems.

For this project, the microbiomes are from the rhizosphere (the soil surrounding the crop plant roots) and from bulk soil (soil outside the rhizosphere).

The Cryobank provides a facility for researchers to source data and samples, including cryo-preserved microbial material and genomic and metagenomic sequences from different soil microbiome environments.

  

The project has sequenced soil microbiomes from 6 different UK crops grown in 9 different soil types (within a pot experiment) from across the United Kingdom.

  

The data catalogue for this project can be accessed at [https://agmicrobiomebase.org/](https://agmicrobiomebase.org/).

The catalogue links the raw sequence data files which have been submitted to the European Nucleotide Archive (ENA) with soil metadata.

  

This repository contains procedural information and scripts for the analysis of Amplicon sequence fastq files to derive amplicon sequence variant taxonomies using qiime2 and the corresponding packages.

1. Analysis of **16S Amplicon** sequence fastq files to derive amplicon sequence variant taxonomies using qiime2

[16S-sequence-analysis.md](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-sequence-analysis.md): A markdown file describing all processing steps from fastq to asv taxa output file

**i**. [16S-script-startup](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-step01-startup.sh): script to install software, set-up environments and downloaded databases

**ii**. [16S-script-preprocessing](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-step02-preprocessing.sh): script to trim reads

**iii**. [16S-script-qiime2](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/16s-step03-qiime2.sh): This script performs the core analysis of 16S Amplicon sequences using qiime2, transforming raw sequence data into meaningful taxonomic information.

[README.md](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-16S/README.md):

2. Analysis of **ITS amplicon** fastq sequence files to derive amplicon sequence variant taxonomies using qiime2

[ITS-sequence-analysis.md](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-ITS/ITS-sequencing-analysis.md): markdown file describing all processing steps from fastq to asv taxa output file

**i**. [ITS-script-startup](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-ITS/ITS-step01-setup.sh) : script to install software and set-up environments and downloaded databases

**ii**. [ITS-script-preprocessing](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-ITS/ITS-step02-preprocessing.sh): script to trim reads

**iii**. [ITS-script-qiime2](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-ITS/ITS-step03-qiime2.sh): This script performs the core analysis of ITS1 Amplicon sequences using qiime2, transforming raw sequence data into meaningful taxonomic information.

[README.md](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/amplicon-ITS/README.md):

3. Preparation and mapping of **amplicon sequence data for submission to ENA using standard templates** for intercative submission via the WebinPortal ([https://www.ebi.ac.uk/ena/submit/webin/login](https://www.ebi.ac.uk/ena/submit/webin/login))

- ENA-upload-procedure.md

- ENA sample template example for one crop

- ENA fastq template example for one crop

- ENA ERS output exmaple file for one crop

- [R markdown](https://github.com/paytonyau/agmicrobiomebase/blob/main/amplicon-sequence-analysis/Fastq_checklist_mapping/fastq_checklist_mapping.Rmd): describing all processing steps to map ERS numbers from sample uplaods to fastq template

- [R-code](): to map ERS numbers from sample uplaods to fastq template

- README.md

4. **Pre-processing of metagenomic sequence data**

  

- procedure.md

- scripts

- README.md

5. Preparation and mapping of **metagenomic sequence data for submission to ENA using standard templates** for intercative submission via the WebinPortal ([https://www.ebi.ac.uk/ena/submit/webin/login](https://www.ebi.ac.uk/ena/submit/webin/login))

- procedure.md

- scripts

- README.md

  

6. Processing and data visualisation of **amplicon sequence data for 3 case studies** in the Crop Microbiome Cryobank Publication 1

**i**. Optimisation/Batch effects corrections

-- procedure.md

-- [scripts](https://github.com/paytonyau/agmicrobiomebase/tree/main/publications-scripts/How-to-paper/case_00-preparation/Case_00-prep.R)

-- README.md

**ii**. Case Study A-Influence of soil type

-- procedure.md

-- [scripts](https://github.com/paytonyau/agmicrobiomebase/tree/main/publications-scripts/How-to-paper/case_01-influence_of_soil_type/Case_01-16s.R)


-- README.md

**iii**. Case Study B The core Microbiome

-- procedure.md

-- [scripts](https://github.com/paytonyau/agmicrobiomebase/tree/main/publications-scripts/How-to-paper/case_02-core_microbiome/Case_02-16s.R)


-- README.md

**iv**. Case Study C ITS taxonomy and Fusarium in wheat

-- procedure.md

-- [scripts](https://github.com/paytonyau/agmicrobiomebase/tree/main/publications-scripts/How-to-paper/case_03-ITS_wheat/Case_03-ITS.R)

-- README.md
