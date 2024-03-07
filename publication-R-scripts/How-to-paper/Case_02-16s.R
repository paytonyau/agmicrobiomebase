### case 02

# Load the tidyverse and reshape2 libraries for data manipulation
# tidyverse includes several packages for data manipulation and visualization
# reshape2 is for data reshaping
library("tidyverse")
library("reshape2")

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


##### Plotting the top 10 taxa at family level #####

# Transform normalised ASVs to proportions
# This step converts the counts in the ASV table to proportions
proportions = transform_sample_counts(physeq.Sugarbeet.group, function(x) 100 * x/sum(x))

# Identify the top 10 most abundant ASVs
# This step identifies the ASVs with the highest total abundance across all samples
Top10ASVs = names(sort(taxa_sums(proportions), TRUE)[1:21])

# Add a new column to the taxonomy table for the top 10 ASVs
# This step prepares for the highlighting of the top 10 ASVs in the subsequent plot
Taxtab10 = cbind(tax_table(proportions), Family10 = NA)
Taxtab10[Top10ASVs, "Family10"] <- as(tax_table(proportions)[Top10ASVs, "Family"], "character")
tax_table(proportions) <- tax_table(Taxtab10)

# Subset the phyloseq object to include only the top 10 ASVs
# This step prepares for the creation of the bar plot
Rgsm10 = prune_taxa(Top10ASVs, proportions)

# Define the color palette for the bar plot
my_palette <- brewer.pal(n = 10, name = "Spectral")

# Create the bar plot
# This step creates a bar plot of the proportion of each of the top 10 ASVs in each sample
p <- plot_bar(Rgsm10, "Soil.Location", fill = "Family") + coord_flip() +
  ylab("Percentage of Sequences") + ylim(0, 20) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = my_palette) +
  labs(x = element_blank()) +
  theme_classic() +
  theme(text = element_text(size=15, colour = "black"),
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.line = element_line(colour = 'black', size = 1.25),
        axis.text.x = element_text(angle=0, hjust=0.5, colour = "black", size = 13),
        axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
        axis.title.y = element_text(color="black", size=15,face="bold"),
        legend.position = "right",
        legend.text = element_text(size = 9.5),
        legend.key.height= unit(0.45, 'cm'),
        legend.key.width= unit(0.45, 'cm')
  )

# Print the bar plot
# This step displays the final bar plot
print(p)


##### Upset plot using UpsetR  #####
# Load the necessary libraries
library("microbiome")
library("speedyseq")
library("UpSetR")
library("plyr")
library("reshape2")
library("RColorBrewer")

# Aggregate taxa at the genus level
# This step aggregates the taxa in the phyloseq object at the genus level
B <- aggregate_taxa(physeq.Sugarbeet.group, "Genus", verbose = TRUE)

# Remove undesired genera
# This step removes genera that are not of interest from the phyloseq object
taxa_to_remove <- c("uncultured", "Unknown")
B2 <- subset_taxa(B, !get("Genus") %in% taxa_to_remove)

# Convert to tibble, rename columns, select relevant columns, group by Sample, and keep top 100 abundant ASVs
# This step prepares the data for the subsequent analyses
D <- as_tibble(B2) %>%
  mutate(Sample = Soil.Location, ASV = .otu, Abundance = .abundance) %>%
  select(Sample, Abundance, Genus) %>%
  group_by(Sample) %>%
  filter(rank(desc(Abundance)) <= 100) %>% # Filter <= 100
  ungroup()

# Remove the Abundance column
D$Abundance <- NULL

# Rename the second column to "ASV"
names(D)[2] <- "ASV"
names(D)[1] <- "Soil.Location"

# Convert data from long to wide format
E <- dcast(D, ASV ~ Soil.Location)

# Define a binary function
binary_fun <- function(x) {
  x[is.na(x)] <- 0
  ifelse(x > 0, 1, 0)
}

# Apply the binary function to columns 2 to 10
temp_df <- apply(E[2:10], 2, binary_fun)
temp_df <- as.data.frame(temp_df)
rownames(temp_df) = E$ASV

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

# Print the ggtree plot
print(upset_plot)

# Clean up by removing unnecessary objects
rm(B, B2, D, E, binary_fun, upset_plot, physeq.Sugarbeet.group, upset_plot)

# Load the ggvenn library
library(ggvenn)

# Extract the rows where the value is 1 for each column
CL.YO <- rownames(temp_df)[temp_df$CL.YO == 1]
CY.YO <- rownames(temp_df)[temp_df$CY.YO == 1]
SC.SH <- rownames(temp_df)[temp_df$SC.SH == 1]
SL.SH <- rownames(temp_df)[temp_df$SL.SH == 1]

# Create a list with the extracted data
list_data <- list("CL.YO" = CL.YO, "CY.YO" = CY.YO, "SC.SH" = SC.SH, "SL.SH" = SL.SH)

# Use ggvenn to create the Venn diagram
Venn <- ggvenn(
  list_data,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)

# Print the Venn plot
print(Venn)
