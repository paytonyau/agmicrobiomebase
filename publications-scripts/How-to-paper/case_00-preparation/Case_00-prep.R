##### The scripts in this .R document contain 3 different analysis:
##### Trimming and optimal sequence length selection                  ######
##### Comparative analysis of taxonomic Levels of reference databases ######
##### Batch effect correction and normalisation                       ###### 
## Author : Payton Yau
## Date : 12-03-2024

##### Trimming and optimal sequence length selection. ######
# We used the longest possible length to find the highest accumulated percentage of the merged sequence length frequency. 
# It assists in trimming out sequences of lower quality and minimises the error rate for the sequence analysis steps in DADA2. 
# This process is part of the events in preprocessing and to prevent “over-merging”.

# The data extracted from "rep-seqs_233_226_3_2-with-phyla-no-mitochondria-no-chloroplast.qzv"

# The denoise setting is
# qiime dada2 denoise-paired \
#   -i-demultiplexed-seqs demux.qza \
#   --p-trunc-len-f 233 \
#   --p-trunc-len-r 226 \
#   --p-trim-left-f 0 \
#   --p-trim-left-r 0 \
#   --p-max-ee-f 3 \
#   --p-max-ee-r 2 \
#   --p-n-threads 8 \
#   --o-representative-sequences rep-seqs.qza \
#   --o-table table.qza

# Load the necessary libraries
library("ggplot2")     
library("scales")     

# Unzip your file
unzip("16s_length_distribution.zip") # This zipped file contains the full list of the merged sequences for ASVs.

# Read the data
# This step reads the data from the .tsv file into a data frame
df = read.csv("16s_length_distribution.tsv", header = T, sep = ",")

# Calculate the length of each sequence
# This step adds a new column to the data frame that contains the length of each sequence
df$SequenceLength <- apply(df, 1, function(x) nchar(x[['Sequence']]))

# Sort the data frame by SequenceLength
# This step sorts the data frame in ascending order of sequence length
df <- df[order(df$SequenceLength), ]

# Calculate the frequency of each sequence length
# This step calculates the frequency of each sequence length and stores the result in a new data frame
df2 <- table(df$SequenceLength)

# Create a data frame for plotting
# This step creates a new data frame that contains the sequence length and its corresponding frequency
plot_df <- data.frame(SequenceLength = as.numeric(names(df2)), Frequency = as.numeric(df2))

# Calculate the accumulated percentage
# This step calculates the accumulated percentage of the sequence length frequency
plot_df$AccumulatedPercentage <- cumsum(plot_df$Frequency) / sum(plot_df$Frequency) * 100

# Filter the data frame
# This step filters the data frame to include only rows where the accumulated percentage is between 0.5% and 99.5%
plot_df2 <- plot_df[plot_df$AccumulatedPercentage >= 0.5 & plot_df$AccumulatedPercentage <= 99.5, ]

# Create a line graph
# This step creates a line graph of the accumulated sequence length frequency
p = ggplot(plot_df2, aes(x = SequenceLength, y = AccumulatedPercentage)) +
  labs(x = "Sequence Length", y = "Accumulated Percentage (%)", 
       title = "Accumulated Sequence Length Frequency (%)") +
  theme_classic() + 
  geom_line(aes(color = Frequency), size = 2) +
  theme(
    text = element_text(size = 15, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold"))

# Add vertical lines and text labels at the 2%, 50%, and 98% points
# This step adds vertical lines and text labels at the 2%, 50%, and 98% points of the accumulated sequence length frequency
p + geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))]), 
               linetype="dashed", size = 1.5, color = "wheat4") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))], y = 2, 
                label = paste("2% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 2))], ")", sep = "")), 
            vjust = -6.5, color = "wheat4") +
  geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))]), 
             linetype="dashed", size = 1.5, color = "forestgreen") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))], y = 50, 
                label = paste("50% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 50))], ")", sep = "")), 
            vjust = -2, color = "forestgreen") +
  geom_vline(aes(xintercept = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))]), 
             linetype="dashed", size = 1.5, color = "coral4") +
  geom_text(aes(x = plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))], y = 98, 
                label = paste("98% (", plot_df$SequenceLength[which.min(abs(plot_df$AccumulatedPercentage - 98))], ")", sep = "")), 
            vjust = 1.5, color = "coral4") +
  scale_x_continuous(breaks = c(400, 405, 410, 415, 420, 425, 430))

# Save the plot to a PDF file
pdf(file = "Fig2B_merged seq_freq.pdf", width = 10,height = 6)
print(p)
dev.off()

##### Comparative analysis of taxonomic Levels of reference databases. ######
# The script provided below is designed to perform a comparative analysis of taxonomic levels of reference databases. 
# It operates by counting the number of matched hits, which are then aggregated by names. 
# In this process, multiple ASVs with the same names are consolidated into a single unique name. 
# This analysis spans across various taxonomic levels, ranging from Kingdom to Species.

##### Install required packages
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library("qiime2R")
library("phyloseq")
library("dplyr")
library("tidyverse")


# Data and Path Definitions

# Define the taxonomic levels
genus_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Define the datasets and paths
# We compared GreenGene v.1 (13.8), GreenGene v.2 (2022.10), and Silva v.138 reference databases.

dataset_info <- list(
  list(  # GreenGenes 1
    features = "[Qiime2]GreenGenes_13_8/428_228_220_table_gg-13-8-with-phyla-no-mitochondria-no-chloroplast.qza",
    taxonomy = "[Qiime2]GreenGenes_13_8/428_228_220_taxonomy_gg-13-8.qza",
    output_path = "[Qiime2]GreenGenes_13_8/level_counts_by_group_gg1.csv"
  ),
  list(  # GreenGenes2
    features = "[Qiime2]GreenGenes2_2022_10/428_228_220_table_gg_2022_10-with-phyla-no-mitochondria-no-chloroplast.qza",
    taxonomy = "[Qiime2]GreenGenes2_2022_10/428_228_220_taxonomy_gg_2022_10.qza",
    output_path = "[Qiime2]GreenGenes2_2022_10/level_counts_by_group_gg2.csv"
  ),
  list(  # Silva138
    features = "[Qiime2]Silva_138/428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza",
    taxonomy = "[Qiime2]Silva_138/428_228_220_taxonomy_silva138.qza",
    output_path = "[Qiime2]Silva_138/level_counts_by_group_silva138.csv"
  )
)


## Data Processing and Analysis

# Loop through each dataset
for (dataset in dataset_info) {
  # Convert qiime2 results to phyloseq format
  physeq <- qza_to_phyloseq(
    features = dataset$features,
    taxonomy = dataset$taxonomy,
    metadata = "16s_meta-table.txt"
  )
  
  physeq.sum <- subset_samples(physeq, Analysis == "Include")
  physeq.sum <- merge_samples(physeq.sum, "Type", fun = sum)
  
  # Create an empty list to store genus-level abundance data for each taxonomic level
  gentab_levels <- list()
  
  # Set observation threshold
  observationThreshold <- 1
  
  # loop through all the taxonomic levels
  for (level in genus_levels) {
    
    # create a factor variable for each level
    genfac <- factor(tax_table(physeq.sum)[, level])
    
    # calculate the abundance of each genus within each sample
    gentab <- apply(otu_table(physeq.sum), MARGIN = 1, function(x) {
      tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
    })
    
    # calculate the number of samples in which each genus is observed above the threshold
    level_counts <- apply(gentab > observationThreshold, 2, sum)
    
    # create a data frame of level counts with genus names as row names
    BB <- as.data.frame(level_counts)
    BB$name <- row.names(BB)
    
    # add the data frame to the gentab_levels list
    gentab_levels[[level]] <- BB
  }
  
  # Combine all level counts data frames into one data frame
  B2 <- gentab_levels %>% reduce(full_join, by = "name")
  
  # Set row names and column names
  row.names(B2) <- B2$name
  B2$name <- NULL
  colnames(B2)[1:7] <- genus_levels
  
  # Write the data frame to a file
  write.csv(B2, file = dataset$output_path, row.names = TRUE)
  
  # Clean up by removing unnecessary objects
  rm(gentab_levels, observationThreshold, BB, B2)
}



## Data Visualization
# Load the reshape2 and ggplot2 libraries
library(reshape2)
library(ggplot2)

GreenGenes.v1 = read.csv("[Qiime2]GreenGenes_13_8/level_counts_by_group_gg1.csv")[5,]
GreenGenes.v2 = read.csv("[Qiime2]GreenGenes2_2022_10/level_counts_by_group_gg2.csv")[5,]
Sliva.v138 = read.csv("[Qiime2]Silva_138/level_counts_by_group_silva138.csv")[5,]

combined_df <- rbind(GreenGenes.v1, GreenGenes.v2, Sliva.v138)
combined_df$X <- c("GreenGenes.v1", "GreenGenes.v2", "Sliva.v138")

data_long <- melt(combined_df, id.vars = "X", variable.name = "Dataset", value.name = "Count")

colnames(data_long) = c("Ref.Database","Taxonomic.Level","Count")

# Convert Taxonomic.Level to a factor and specify the desired order of the levels
data_long$Taxonomic.Level <- factor(data_long$Taxonomic.Level,
                                    levels = c("Kingdom", "Phylum", "Class", "Order", 
                                               "Family", "Genus", "Species"))


# Plot the data as a line graph using ggplot
# Open a new PDF graphics device
pdf(file = "line_graph.pdf", width=8,height=5)

ggplot(data_long, aes(x = Taxonomic.Level, y = Count, color = Ref.Database, group = Ref.Database)) +
  geom_line(size = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Sliva.v138" = "cornflowerblue", 
                                "GreenGenes.v2" = "greenyellow", 
                                "GreenGenes.v1" = "forestgreen")) +
  labs(x = "Taxonomic Level", y = "Count", color = "Reference\nDatabase") +
  theme_classic() + 
  theme(
    text = element_text(size = 19, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(size = 13.5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size=unit(0.4,"cm")
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_y_continuous(breaks=seq(0,1500,by=250))


# Close the PDF device and save the plot to a file
dev.off()  


###### Batch effect correction and normalisation ######

# Given that the sequencing was conducted in distinct batches, with our data divided into four separate runs, 
# it was necessary to adjust for batch variation using a correction model.
# ConQuR was selected based on conditional quantile regression. 

# Load necessary libraries
library("qiime2R")
library("phyloseq")
library("RColorBrewer")

# Convert qiime2 results to phyloseq format
# This step converts the qiime2 output files into a phyloseq object
physeq <- qza_to_phyloseq(
  features = "~/GitHub/agmicrobiomebase/amplicon-sequence-analysis/amplicon-16S/[Qiime2]Silva_138/428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
  taxonomy = "~/GitHub/agmicrobiomebase/amplicon-sequence-analysis/amplicon-16S/[Qiime2]Silva_138/428_228_220_taxonomy_silva138.qza",
  metadata = "16s_meta-table.txt"
  # , tree = "rooted-tree.qza"
)

# Print the phyloseq object to confirm its creation
physeq

# Remove unwanted (failed and controls) samples before the normalisation
# This step removes samples that are not included in the analysis
physeq.ori <- subset_samples(physeq, Analysis == "Include")

# Remove the original phyloseq object to save memory
rm(physeq)

# Load the ConQuR library for batch effects correction
library(ConQuR)

# Load the doParallel library for parallel computing
library(doParallel)

# Convert ASV table to a data frame and transpose
# This step prepares the ASV table for the ConQuR function
B <- as.data.frame(physeq.ori@otu_table) # taxa
B <- t(B)
B <- as.data.frame(B)

# Extract batch ID from sample data
# This step prepares the batch ID for the ConQuR function
batchid = physeq.ori@sam_data$Plate # batchid

# Extract covariates
# This step prepares the covariates for the ConQuR function
D = physeq.ori@sam_data[, c('Type', 'Soil', 'Location')] #covar
summary(D)

# Correct for batch effects using ConQuR package
# This step corrects for batch effects in the ASV table
options(warn=-1) # required to call
taxa_correct1 = ConQuR(tax_tab = B,
                       batchid = batchid,
                       covariates = D,
                       batch_ref="1"
) # warning messages may appear & it can be ignored

# Transpose the corrected matrix and convert it to a data frame
taxa_correct2 <- t(taxa_correct1)
taxa_correct2 <- as.data.frame(taxa_correct2)

# Create new ASV table, taxonomy table, and sample data
ASV = otu_table(taxa_correct2, taxa_are_rows = TRUE)
TAXA = tax_table(physeq.ori)
sampledata = sample_data(physeq.ori)

# Repack the objects into a level 4 phyloseq structural data
physeq.norm = phyloseq(ASV, TAXA, sampledata)

# Remove unnecessary objects to save memory
rm(B, D, batchid, taxa_correct1, taxa_correct2, ASV, TAXA, sampledata, to_skip)

### Bata diversity - before and after the normalisation

# Load required packages
# install.packages("ggplot2")
library("ggplot2") # ggplot2 is a plotting system for R, based on the grammar of graphics.
# install.packages("dplyr")
library("dplyr") # dplyr is a grammar of data manipulation, providing a consistent set of verbs that help you solve the most common data manipulation challenges.
# install.packages("ggpubr")
library("ggpubr") # ggpubr provides some easy-to-use functions for creating and customizing 'ggplot2'- based publication ready plots.

# Beta diversity analysis - based on Plate and Type - before the normalisation
# ordinate function performs NMDS ordination on the original phylogenetic sequence data (physeq.ori) using Bray-Curtis dissimilarity.
NMDS1 <- ordinate(
  physeq = physeq.ori,
  method = "NMDS",
  distance = "bray"
)

# Plot ordination
# This block of code creates an NMDS plot of the original phylogenetic sequence data, with points colored by 'Plate' and shaped by 'Type'. The plot is customized with various themes and scales.
plot_ordination(
  physeq = physeq.ori,
  ordination = NMDS1,
  color = "Plate",
  shape = "Type"
) +
  theme_classic() +
  geom_point(aes(color = Plate), alpha = 1, size = 3.5) +
  theme(
    text = element_text(size = 18, colour = "black"),
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 17, 3, 4, 16, 18, 21, 22, 23)) # Set custom shapes

# Beta diversity analysis - based on Plate and Type - after the normalisation
# ordinate function performs NMDS ordination on the normalized phylogenetic sequence data (physeq.norm) using Bray-Curtis dissimilarity.
NMDS2 <- ordinate(
  physeq = physeq.norm,
  method = "NMDS",
  distance = "bray"
)

# Plot ordination
# This block of code creates an NMDS plot of the normalized phylogenetic sequence data, with points colored by 'Plate' and shaped by 'Type'. The plot is customized with various themes and scales.
plot_ordination(
  physeq = physeq.norm,
  ordination = NMDS2,
  color = "Plate",
  shape = "Type"
) +
  theme_classic() +
  geom_point(aes(color = Plate), alpha = 1, size = 3.5) +
  theme(
    text = element_text(size = 18, colour = "black"),
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 17, 3, 4, 16, 18, 21, 22, 23)) # Set custom shapes

# Save the "physeq.norm" object for the other case study analysis
save(physeq.norm, file = "norm.RData")