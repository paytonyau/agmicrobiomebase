## The core microbiome (Sugarbeet vs Bulksoil)
## Author : Payton Yau
## Date : 12-03-2024

# To study rhizosphere microbiota, we identified the minimal microbiota interacting with a specific crop, 
# sugar beet, in natural conditions across various soils. 
# This approach helps understand real-world microbial interactions in the rhizosphere.

library("RColorBrewer")

load("physeq.Sugarbeet.group.RData")

##### Plotting the top 10 taxa at family level #####
# The script required from the use of "physeq.Sugarbeet.group" which is from`case_01-16s.R`
 
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

col = brewer.pal(n = 9, name = "Set3")

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
rm(B, B2, D, E, binary_fun, upset_plot, physeq.Sugarbeet.group)

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
