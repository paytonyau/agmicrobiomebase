### case 03 - ITS

# Install and load the required packages for the conversion
# Uncomment the lines to install the packages if not already installed
# install.packages("devtools")
# devtools::install_github("jbisanz/qiime2R")
library("qiime2R")

# install.packages("BiocManager")
# BiocManager::install("phyloseq")
library("phyloseq")

# install.packages("tidyverse")
library("tidyverse")

# install.packages("ggpubr")
library("ggpubr")

# The conversion of Qiime2 to Phyloseq datastructure
# Convert qiime2 to phyloseq format
physeq.a <- qza_to_phyloseq(
  features = "~/GitHub/agmicrobiomebase/amplicon-sequence-analysis/amplicon-ITS/[Qiime2]UNITE_9/table-its-fungi-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
  tree = "~/GitHub/agmicrobiomebase/amplicon-sequence-analysis/amplicon-ITS/[Qiime2]UNITE_9/phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza",
  taxonomy = "~/GitHub/agmicrobiomebase/amplicon-sequence-analysis/amplicon-ITS/[Qiime2]UNITE_9/taxonomy-fungi.qza",
  metadata = "ITS-meta-table.txt"
)

# Print the phyloseq object to confirm the conversion
print(physeq.a)

# Merge the replicate samples for each Group
physeq.a.group = merge_samples(physeq.a, "Group") # Sum between replicate samples

# Repair factors in the sample metadata
sample_data(physeq.a.group)$Group <- levels(sample_data(physeq.a)$Group)[get_variable(physeq.a.group, "Group")]
sample_data(physeq.a.group)$Soil.Type <- levels(sample_data(physeq.a)$Soil.Type)[get_variable(physeq.a.group, "Soil.Type")]
sample_data(physeq.a.group)$Soil <- levels(sample_data(physeq.a)$Soil)[get_variable(physeq.a.group, "Soil")]

# Assessing ASV distribution (merged and filtered)
# We evaluated the number of Amplicon Sequence Variants (ASVs) and visualised the raw read counts per sample within each group. 
# This provided insights into the distribution of ASVs and the sequencing depth of each sample, allowing for a comprehensive assessment of data quality and quantity.

# install.packages("ggplot2")
library("ggplot2")

# install.packages("RColorBrewer")
library("RColorBrewer")

# Calculate the sum of ASVs for each sample
ASV <- sample_sums(physeq.a)
ASV <- as.data.frame(ASV)
ASV$Group <- physeq.a@sam_data$Group

# Create a scatter plot of the ASV counts for each group
ggplot(ASV, aes(x = Group, y = ASV,  colour = Group))+
  scale_color_brewer(palette = "Paired") +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() + 
  theme(text = element_text(size=15, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=45, hjust=1, 
                                   colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, 
                                   colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=15,face="bold"),
        legend.position = "none")

# Remove the ASV data frame to free up memory
rm(ASV)


##### Beta diversity #####
# Beta diversity analysis
# This section of the script performs a beta diversity analysis, which compares the microbial community composition between samples

# Perform NMDS ordination
# This step performs non-metric multidimensional scaling (NMDS) ordination on the phyloseq object
# The ordination is based on the Bray-Curtis dissimilarity between samples
NMDS <- ordinate(physeq = physeq.a, method = "NMDS", distance = "bray")

# Plot the ordination
# This step creates a scatter plot of the ordination with points colored and shaped by the "Group" variable
# The plot also includes ellipses around the points of each group
plot_ordination(physeq = physeq.a,
                ordination = NMDS,
                color = "Group",
                shape = "Group"
) +
  theme_classic() + 
  scale_color_brewer(palette = "Paired") + 
  scale_fill_brewer(palette = "Paired") + 
  geom_point(aes(color = Group), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom = "polygon", type="norm", alpha=0.25, aes(fill = Group)) +
  scale_shape_manual(values = c(15, 17, 3, 4, 16, 18, 21, 22, 23)) # Set custom shapes

# Clean up by removing objects that are no longer needed
# This step removes the NMDS object to free up memory
rm(NMDS)


########  alpha diversity (Shannon) ########

# Alpha diversity analysis
# This section of the script performs an alpha diversity analysis, which measures the diversity within each sample

# Calculate alpha diversity (Shannon index) for each sample
# This step calculates the Shannon diversity index for each sample in the phyloseq object
alpha.object <- cbind(
  x = sample_data(physeq.a),
  y = estimate_richness(physeq.a, measures = 'Shannon')
)

# Create a scatter plot of the Shannon diversity index for each group
# This step creates a scatter plot of the Shannon diversity index for each group
# The points are jittered to avoid overplotting, and a boxplot is added to show the distribution of the index within each group
ggplot(data = alpha.object, aes(x = x.Group, y = Shannon, color = x.Group)) + 
  scale_color_brewer(palette = "Paired") +
  theme_classic() + 
  labs(
    x = element_blank(),                     # No x-axis label
    y = "Alpha Diversity (Shannon)"          # y-axis label
  ) + 
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, 
                               hjust = 1, size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, 
                               colour = "black", size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "none"                 # Hide legend
  ) 

# Clean up by removing the alpha.object
# This step removes the alpha.object to free up memory
rm(alpha.object)

##### Number of taxa in groups #####

# Load the necessary libraries
library("dplyr")
library("reshape2")

# Merge the replicate samples for each Group
physeq <- merge_samples(physeq.a, "Group", fun = sum)

# Initialize an empty list to store genus-level abundance data for each taxonomic level
gentab_levels <- list()

# Set observation threshold
observationThreshold <- 1

# Define the taxonomic levels
genus_levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Loop through all the taxonomic levels
for (level in genus_levels) {
  # Create a factor variable for each level
  genfac <- factor(tax_table(physeq)[, level])
  
  # Calculate the abundance of each genus within each sample
  gentab <- apply(otu_table(physeq), MARGIN = 1, function(x) {
    tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
  })
  
  # Calculate the number of samples in which each genus is observed above the threshold
  level_counts <- apply(gentab > observationThreshold, 2, sum)
  
  # Create a data frame of level counts with genus names as row names
  BB <- as.data.frame(level_counts)
  BB$Group <- row.names(BB)
  
  # Add the data frame to the gentab_levels list
  gentab_levels[[level]] <- BB
}

# Combine all level counts data frames into one data frame
B2 <- gentab_levels %>% reduce(full_join, by = "Group")

# Set row names and column names
row.names(B2) <- B2$Group
B2$Group <- NULL
colnames(B2)[1:6] <- genus_levels

# Print the resulting data frame
print(B2)

# Convert the data frame to long format for plotting
B2$Group = row.names(B2)
data_long <- reshape2::melt(B2, id.vars = "Group", 
                            variable.name = "Level", 
                            value.name = "Taxa")

# Plot the data as a line graph using ggplot
ggplot(data_long, aes(x = Level, y = Taxa, 
                      color = Group, group = Group)) +
  geom_line() +
  geom_point(size = 4) +
  labs(x = "Taxonomic Level", y = "Count", color = "Group") +
  theme_classic() +   scale_color_brewer(palette = "Paired") +
  theme(
    text = element_text(size = 19, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(size = 13.5, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size=unit(0.4,"cm")
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_y_continuous(breaks=seq(0,400,by=100))

# Clean up by removing unnecessary objects
rm(physeq, gentab_levels, observationThreshold, B2, data_long, genfac, genus_levels, level, level_counts)


##### Upset plot using UpsetR #####

# Load the necessary libraries
library("UpSetR")
library("reshape2")
library("plyr")
library("dplyr")
library("microbiome")

# Aggregate taxa at the genus level
# This step aggregates the taxa in the phyloseq object at the genus level
B <- aggregate_taxa(physeq.a.group, "Genus", verbose = TRUE)

# Remove undesired genera
# This step removes genera that are not of interest from the phyloseq object
taxa_to_remove <- c("uncultured", "Unknown")
B2 <- subset_taxa(B, !get("Genus") %in% taxa_to_remove)

# Extract relevant data from the phyloseq object
sample_data <- sample_data(B2)
otu_table <- otu_table(B2)
abundance <- as.vector(otu_table)

# Create a tibble with the extracted data
D <- tibble(
  Sample = rep(sample_data$Group, each = nrow(otu_table)),
  ASV = rep(rownames(otu_table), times = ncol(otu_table)),
  Abundance = abundance
) %>%
  group_by(Sample) %>%
  mutate(rank = rank(desc(Abundance))) %>%
  filter(Abundance > 0) %>%
  ungroup() %>%
  select(Sample, Abundance, ASV)

# Remove the Abundance column
D$Abundance <- NULL

# Rename the second column to "ASV"
names(D)[2] <- "ASV"
names(D)[1] <- "Group"

# Convert data from long to wide format
E <- dcast(D, ASV ~ Group)

# Define a binary function
binary_fun <- function(x) {
  x[is.na(x)] <- 0
  ifelse(x > 0, 1, 0)
}

col = brewer.pal(n = 9, name = "Paired")

# Apply the binary function to columns 2 to 10
temp_df <- apply(E[2:10], 2, binary_fun)
temp_df <- as.data.frame(temp_df)

# Create an UpSet plot
upset_plot <- upset(temp_df, 
                    sets = colnames(temp_df), 
                    sets.bar.color = (col),
                    order.by = "freq", 
                    empty.intersections = "on",
                    mainbar.y.label = "Counts by Pattern of Conditions", 
                    sets.x.label = "Counts by Condition",
                    matrix.color="blue", 
                    mb.ratio = c(0.6, 0.4),
                    point.size= 2.75,
                    line.size = 1.25, 
                    text.scale = 1.5
)

# Print the UpSet plot
print(upset_plot)


##### Phylogenetic tree ######

# Load the necessary libraries
library("phyloseq")   # For handling phylogenetic sequencing data
library("ggtree")     # For tree visualisation
library("scales")     # For scaling transformations

# Load the GlobalPatterns dataset and prune taxa
# This step removes taxa with zero sums and subsets data based on Phylum value
GP <- prune_taxa(taxa_sums(physeq.a) > 0, physeq.a)  
GP.chl <- subset_taxa(GP, Genus == "Fusarium")   

# Create a ggtree plot
# This step creates a phylogenetic tree plot with points colored and shaped by the "Group" variable
p <- ggtree(GP.chl, ladderize = TRUE) +
  geom_tiplab(aes(label = Species), as_ylab=TRUE) +
  geom_point(aes(x = x + hjust, color = Group, 
                 shape = Group, size = Abundance), na.rm = TRUE) +
  scale_size_continuous(trans = log_trans(2)) +
  scale_shape_manual(values = c(15, 17, 3, 4, 16, 18, 21, 22, 23)) + # Set custom shapes
  scale_color_brewer(palette = "Paired") +
  theme(text = element_text(size=13, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25) ,
        axis.title.y = element_text(color="black", size=2.5,face="bold"), 
        legend.text = element_text(size = 8),
        legend.key.size=unit(0.4,"cm"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom"
  ) 

# Print the ggtree plot
print(p)

# Subset the taxa to Genus from physeq.wheat 
# This step subsets the data to include only the genus "Fusarium" and its species
physeq.a.genus <- subset_taxa(physeq.a, Genus == "Fusarium")
physeq.a.equiseti <- subset_taxa(physeq.a, Species == "Fusarium_equiseti")
physeq.a.nurragi <- subset_taxa(physeq.a, Species == "Fusarium_nurragi")
physeq.a.waltergamsii <- subset_taxa(physeq.a, Species == "Fusarium_waltergamsii")

# Calculate the total abundance of Fusarium for each sample
# This step calculates the total abundance of each Fusarium species in each sample
meta = physeq.a@sam_data
otudf = as.data.frame(t(as.data.frame(physeq.a.genus@otu_table)))
meta$Fusarium = rowSums(otudf)
otudf = as.data.frame(t(as.data.frame(physeq.a.equiseti@otu_table)))
meta$F.equiseti = rowSums(otudf)
otudf = as.data.frame(t(as.data.frame(physeq.a.nurragi@otu_table)))
meta$F.nurragi = rowSums(otudf)
otudf = as.data.frame(t(as.data.frame(physeq.a.waltergamsii@otu_table)))
meta$F.waltergamsii = rowSums(otudf)

# Plot a graph of the abundance of Fusarium for each sample grouped by Group:
# This step creates a scatter plot of the abundance of Fusarium and its species for each sample, grouped by the "Group" variable
p1 <- ggplot(subset(meta, Group %in% c("CL.BO","CL.YO",
                                       "CY.BU","CY.YO",
                                       "SC.HE","SC.SH",
                                       "SL.AN","SL.BE",
                                       "SL.SH")),
             aes(x = Group, y = Fusarium,  colour = interaction(Group))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   scale_color_brewer(palette = "Paired") +
  labs(x = "", y = "\n Fusarium") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=45, hjust=1, 
                                   colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, 
                                   colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=12,face="bold"), 
        legend.position = "none")

p2 <- ggplot(subset(meta, Group %in% c("CL.BO","CL.YO",
                                       "CY.BU","CY.YO",
                                       "SC.HE","SC.SH",
                                       "SL.AN","SL.BE",
                                       "SL.SH")),
             aes(x = Group, y = F.equiseti,  colour = interaction(Group))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   scale_color_brewer(palette = "Paired") +
  labs(x = "", y = "\n F.nurragi") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=45, hjust=1, 
                                   colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, 
                                   colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=12,face="bold"), 
        legend.position = "none")

p3 <- ggplot(subset(meta, Group %in% c("CL.BO","CL.YO",
                                       "CY.BU","CY.YO",
                                       "SC.HE","SC.SH",
                                       "SL.AN","SL.BE",
                                       "SL.SH")),
             aes(x = Group, y = F.nurragi,  colour = interaction(Group))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   scale_color_brewer(palette = "Paired") +
  labs(x = "", y = "\n F.nurragi") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=45, hjust=1, 
                                   colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, 
                                   colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=12,face="bold"), 
        legend.position = "none")

p4 <- ggplot(subset(meta, Group %in% c("CL.BO","CL.YO",
                                       "CY.BU","CY.YO",
                                       "SC.HE","SC.SH",
                                       "SL.AN","SL.BE",
                                       "SL.SH")),
             aes(x = Group, y = F.waltergamsii,  colour = interaction(Group))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   scale_color_brewer(palette = "Paired") +
  labs(x = "", y = "\n F.waltergamsii") +
  theme(text = element_text(size=15, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=45, hjust=1, 
                                   colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, 
                                   colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=12,face="bold"), 
        legend.position = "none") 

# Combine and Arrange the plots
fig <- ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), size = 8,
                 ncol = 2, nrow = 2)

# Add labels
fig <- annotate_figure(fig)

# Print the figure
print(fig)
