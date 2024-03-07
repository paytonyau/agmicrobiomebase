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
