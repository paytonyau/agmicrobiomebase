## The influence of soil type (Sugarbeet vs Bulksoil)
## Author : Payton Yau
## Date : 12-03-2024

# Sugar beet, a member of the brassica family, 
# is cultivated in the UK and other regions for its tuber, which is rich in sucrose. 
# We focused our study on this single crop type to examine the influence of soil type and locale. 
# To investigate the rhizosphere effect, which refers to the selection of specific members from the source soil microbiota, 
# we compared the diversity of the sugar beet rhizosphere to that of a no-plant control of bulk soil.

# Load the tidyverse and reshape2 libraries for data manipulation
# tidyverse includes several packages for data manipulation and visualisation
# reshape2 is for data reshaping
library("tidyverse")
library("reshape2")

load("norm.RData") # load physeq.norm

# Subset the phyloseq object to include only samples of type "Sugarbeet" and "Bulksoil"
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Sugarbeet", "Bulksoil"))

# Subset the phyloseq object to include only samples of type "Sugarbeet"
physeq.Sugarbeet <- subset_samples(physeq.norm, Type=="Sugarbeet")

# Merge the replicate samples for each Group
# This step combines the replicate samples into one sample per group
physeq.Sugarbeet.group = merge_samples(physeq.Sugarbeet, "Group")

# Repair factors in the sample metadata
# This step ensures that the metadata factors match the levels in the original phyloseq object
sample_data(physeq.Sugarbeet.group)$Group <- rownames(sample_data(physeq.Sugarbeet.group))
sample_data(physeq.Sugarbeet.group)$Soil <- levels(sample_data(physeq.norm)$Soil)[get_variable(physeq.Sugarbeet.group, "Soil")]
sample_data(physeq.Sugarbeet.group)$Soil.Location <- levels(sample_data(physeq.norm)$Soil.Location)[get_variable(physeq.Sugarbeet.group, "Soil.Location")]

# Further subset the phyloseq object to include only samples from groups "SU.CL.BO" and "SU.CL.YO"
physeq.Sugarbeet.vs <- physeq.Sugarbeet %>% subset_samples(Group %in% c("SU.CL.BO", "SU.CL.YO"))

# Perform NMDS ordination on the subsetted phyloseq object
# This step creates an ordination object using the Bray-Curtis distance measure
NMDS <- ordinate(physeq = physeq.SU,
                 method = "NMDS",
                 distance = "bray"
)

# Define the groups to be ellipsed in the ordination plot
groups_to_ellipse <- c("SC.HE", "CY.BU") 

# Subset the phyloseq object to include only samples from the defined groups
physeq.SU.subset <- subset_samples(physeq.SU, Soil.Location %in% groups_to_ellipse)

# Convert the sample data of the subsetted phyloseq object to a data frame
df <- sample_data(physeq.SU.subset)

# Define a color palette for the ordination plot
my_palette <- c("darkgoldenrod", "limegreen")

# Create the ordination plot
# This step creates a scatter plot of the ordination with points colored by "Type" and shaped by "Soil.Location"
plot_ordination <- plot_ordination(physeq = physeq.SU,
                                   ordination = NMDS,
                                   color = "Type",
                                   shape = "Soil.Location")

# Extract the ordination scores from the plot
df <- plot_ordination$data

# Subset the data for the groups to be ellipsed
df_subset <- df[df$Soil.Location %in% groups_to_ellipse, ]

# Add the ellipse to the ordination plot
# This step adds an ellipse around the points of each group in the subsetted data
plot_ordination <- plot_ordination +
  stat_ellipse(data = df_subset,
               type="norm",
               alpha=0.25,
               aes(group = Soil.Location),
               linetype = 1,
               size = 0.8,
               colour = "purple4"
  ) +
  theme_classic() +
  geom_point(aes(color = Type), alpha = 1, size = 3.5) +
  theme(
    text = element_text(size = 18, colour = "black"),
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5,
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5,
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold")) +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_shape_manual(values = c(15, 17, 3, 4, 16, 18, 21, 22, 23)) # Set custom shapes

# Plot the ordination
plot_ordination


###### Alpha diversity #######

# Available measurements: "Observed" "Chao1" "ACE" "Shannon" "Simpson" "InvSimpson" "Fisher"

# Calculate alpha diversity (Shannon) and store it in physeq.SU object
# This step calculates the Shannon diversity index for each sample in the phyloseq object
alpha.object <- cbind(
  x = sample_data(physeq.SU),
  y = estimate_richness(physeq.SU, measures = 'Shannon')
)

# Data preparation (formatting)
# This step prepares the data for the subsequent analyses
selected_columns <- alpha.object[, c("x.Soil.Location", "x.Group", "x.Type", "Shannon")]
selected_columns2 <- melt(selected_columns)
names(selected_columns2) <- c("Soil_Location","Group" , "Type", "variable", "Shannon")

# Define the comparisons
# This step defines the pairs of groups to be compared in the subsequent analyses
my_comparisons <- list(
  c("SU.CL.BO", "CO.CL.BO"),
  c("SU.CL.YO", "CO.CL.YO"),
  c("SU.CY.BU", "CO.CY.BU"),
  c("SU.CY.YO", "CO.CY.YO"),
  c("SU.SC.HE", "CO.SC.HE"),
  c("SU.SC.SH", "CO.SC.SH"),
  c("SU.SL.AN", "CO.SL.AN"),
  c("SU.SL.BE", "CO.SL.BE"),
  c("SU.SL.SH", "CO.SL.SH")
)

# Initialise an empty data frame to store the results
# This step prepares for the storage of the results from the subsequent analyses
results <- data.frame()

# Perform t-tests for each pair of groups
# This step performs a t-test for each pair of groups defined above
for (i in seq_along(my_comparisons)) {
  group1_data <- selected_columns2$Shannon[selected_columns2$Group == my_comparisons[[i]][1]]
  group2_data <- selected_columns2$Shannon[selected_columns2$Group == my_comparisons[[i]][2]]
  wilcox_test_result <- wilcox.test(group1_data, group2_data)
  results <- rbind(results, data.frame(
    group1 = my_comparisons[[i]][1],
    group2 = my_comparisons[[i]][2],
    p.value = wilcox_test_result$p.value
  ))
}

# Adjust the p-values for multiple comparisons using the Benjamini-Hochberg procedure
# This step corrects for multiple testing by adjusting the p-values
results$p.adjusted <- p.adjust(results$p.value, method = "BH")

# Add significance levels based on the adjusted p-values
# This step adds a column to the results data frame indicating the significance level of each comparison
results$p.signif <- symnum(results$p.adjusted, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                           symbols = c("***", "**", "*", ".", " "))

# Print the results
print(results)

# Define the color palette for the boxplot
my_palette <- c("darkgoldenrod", "limegreen")

# Create the boxplot
# This step creates a boxplot of the Shannon diversity index for each soil location, colored by type
p <- ggplot(data=selected_columns2, aes(x=Soil_Location, y= Shannon, fill = Type)) +
  geom_boxplot(size = 1.1,
               width = 0.825,
               color = "grey20",
               position = position_dodge(0.9)) +
  scale_fill_manual(values = my_palette) +
  labs(x = element_blank(),
       y = "Alpha Diversity (Shannon)"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 18, colour = "black"),
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 1,
                               size = 13, face = "bold"),
    axis.text.y = element_text(angle = 0, hjust = 0, colour = "black",
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    legend.position = "right") +
  scale_y_continuous(breaks = seq(6, 8, by = 0.5), limits = c(6, 8.2))

# Add the results of the comparisons to the plot
# This step adds the adjusted p-values from the comparisons to the boxplot
for (i in seq_len(nrow(results))) {
  # Set the y position for the label
  y_position <- 8.15
  # Add the label to the plot
  p <- p + annotate("text", x = i,
                    y=y_position ,
                    label=round(results$p.adjusted[i], 4),
                    size= 3.75, face = "bold")
}

# Print the boxplot
print(p)




##### Determine the count of taxa within each level and group #####

# Create an empty list to store genus-level abundance data for each taxonomic level
# This step prepares for the storage of the abundance data
gentab_levels <- list()

# Set observation threshold
# This step sets the minimum abundance for a genus to be considered in the analysis
observationThreshold <- 1

# Define the taxonomic levels
# This step defines the taxonomic levels to be considered in the analysis
genus_levels <- c("Kingdom", "Phylum", "Class", "Order",
                  "Family", "Genus", "Species")

# Loop through all the taxonomic levels
# This step performs the analysis for each taxonomic level
for (level in genus_levels) {
  # Create a factor variable for each level
  # This step prepares the data for the calculation of genus abundance
  genfac <- factor(tax_table(physeq.Sugarbeet.group)[, level])
  
  # Calculate the abundance of each genus within each sample
  # This step calculates the total abundance of each genus in each sample
  gentab <- apply(otu_table(physeq.Sugarbeet.group), MARGIN = 1, function(x) {tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
  })
  
  # Calculate the number of samples in which each genus is observed above the threshold
  # This step calculates the number of samples in which each genus is observed with an abundance above the threshold
  level_counts <- apply(gentab > observationThreshold, 2, sum)
  
  # Create a data frame of level counts with genus names as row names
  # This step prepares the results for storage in the gentab_levels list
  BB <- as.data.frame(level_counts)
  BB$name <- row.names(BB)
  
  # Add the data frame to the gentab_levels list
  # This step stores the results in the gentab_levels list
  gentab_levels[[level]] <- BB
}

# Combine all level counts data frames into one data frame
# This step combines the results from all taxonomic levels into one data frame
B2 <- gentab_levels %>% reduce(full_join, by = "name")

# Set row names and column names
# This step sets the row names and column names of the combined data frame
rownames(B2) <- B2$name
B2$name <- NULL
colnames(B2)[1:7] <- genus_levels

# Print the resulting data frame
# This step prints the final results
print(B2)

# Clean up by removing unnecessary objects
# This step removes the objects that are no longer needed to free up memory
rm(gentab_levels, BB)


##### Pairwise comparison using PERMANOVA #####

# Load the necessary libraries
# pairwiseAdonis is used for pairwise comparisons using PERMANOVA
# GGally is used for creating plots
library("pairwiseAdonis")
library("GGally")

# Convert the sample data and OTU table of the phyloseq object to data frames
metdat = as.data.frame(as.matrix(physeq.SU@sam_data))
dat = as.data.frame(t(as.data.frame(physeq.SU@otu_table)))

# Define the pairs for comparison
# These are the pairs of groups that will be compared in the subsequent analyses
pairs <- list(
  c("SU.CL.BO", "CO.CL.BO"),
  c("SU.CL.YO", "CO.CL.YO"),
  c("SU.CY.BU", "CO.CY.BU"),
  c("SU.CY.YO", "CO.CY.YO"),
  c("SU.SC.HE", "CO.SC.HE"),
  c("SU.SC.SH", "CO.SC.SH"),
  c("SU.SL.AN", "CO.SL.AN"),
  c("SU.SL.BE", "CO.SL.BE"),
  c("SU.SL.SH", "CO.SL.SH")
)

# Initialize an empty list to store the results
# This step prepares for the storage of the results from the subsequent analyses
results <- list()

# Loop over each pair
# This step performs the analysis for each pair of groups
for(i in seq_along(pairs)) {
  pair <- pairs[[i]]
  dat$Group = metdat$Group
  dat_subset <- dat[dat$Group %in% pair, ]
  dat_subset$Group <- NULL
  metdat_subset <- metdat[metdat$Group %in% pair, ]
  
  # Perform the pairwise comparison
  # This step performs a pairwise comparison using PERMANOVA for the current pair of groups
  results[[i]] <- pairwise.adonis(dat_subset,
                                  metdat_subset$Group,
                                  sim.function = "vegdist",
                                  sim.method = "bray",
                                  reduce = NULL, perm = 100000)
}

# Convert the list of results to a data frame
# This step combines the results from all pairwise comparisons into one data frame
results_df <- do.call(rbind, lapply(results, function(x) data.frame(t(unlist(x)))))

# Add the adjusted p-value
# This step corrects for multiple testing by adjusting the p-values using the Benjamini-Hochberg procedure
results_df$p.adjusted <- p.adjust(results_df$p.value, method = "BH")

# Print the resulting data frame
# This step prints the final results
print(results_df)

# Save the "physeq.Sugarbeet.group" object for the other case study analysis
save(physeq.Sugarbeet.group, file = "physeq.Sugarbeet.group.RData")

