# R version 4.2.0
### Convert qiime2 results to phyloseq format


##### Install required packages
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library("qiime2R")

library("phyloseq")
library("dplyr")
library("tidyverse")
#################################################
# Define the taxonomic levels
genus_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Define the datasets and paths
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

# Loop through each dataset
for (dataset in dataset_info) {
  # Convert qiime2 results to phyloseq format
  physeq <- qza_to_phyloseq(
    features = dataset$features,
    taxonomy = dataset$taxonomy,
    metadata = "meta-table.txt"
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

###########################################################
# Load the reshape2 library
library(reshape2)
library(ggplot2)

GreenGenes.1 = read.csv("[Qiime2]GreenGenes_13_8/level_counts_by_group.csv")[5,]
GreenGenes.2 = read.csv("[Qiime2]GreenGenes2_2022_10/level_counts_by_group_gg2.csv")[5,]
sliva.138 = read.csv("[Qiime2]Silva_138//level_counts_by_group.csv")[5,]

combined_df <- rbind(GreenGenes.1, GreenGenes.2, sliva.138)
combined_df$X <- c("GreenGenes.1", "GreenGenes.2", "sliva.138")

data_long <- melt(combined_df, id.vars = "X", variable.name = "Dataset", value.name = "Count")

colnames(data_long) = c("Database","Taxonomic.Level","Count")

# Convert Taxonomic.Level to a factor and specify the desired order of the levels
data_long$Taxonomic.Level <- factor(data_long$Taxonomic.Level,
                               levels = c("Kingdom", "Phylum", "Class", "Order", 
                                          "Family", "Genus", "Species"))


# Plot the data as a line graph using ggplot
# Open a new PDF graphics device
pdf(file = "line_graph.pdf", width=8,height=5)
ggplot(data_long, aes(x = Taxonomic.Level, y = Count, color = Database, group = Database)) +
  geom_line() +
  geom_point(size = 4) +
  scale_color_manual(values = c("sliva.138" = "red", 
                                "GreenGenes.2" = "blue", 
                                "GreenGenes.1" = "green")) +
  labs(x = "Taxonomic Level", y = "Count") +
  theme_classic() + 
  theme(
    text = element_text(size = 19, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 14, face = "bold")
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_y_continuous(breaks=seq(0,1500,by=250))


# Close the PDF device and save the plot to a file
dev.off()  

############################################## individual samples ######################
# Create a factor corresponding to the Species
genfac = factor(tax_table(physeq)[, "Species"])

# Tabulate the counts for each Species in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1
  
Species = apply(gentab > observationThreshold, 2, sum)
BB = as.data.frame(Species)

#######################
# Create a factor corresponding to the Genus
genfac = factor(tax_table(physeq)[, "Genus"])

# Tabulate the counts for each Genus in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Genus = apply(gentab > observationThreshold, 2, sum)
BB$Genus = Genus

#####################
# Create a factor corresponding to the Family
genfac = factor(tax_table(physeq)[, "Family"])

# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Family = apply(gentab > observationThreshold, 2, sum)
BB$Family = Family

###################
# Create a factor corresponding to the Order
genfac = factor(tax_table(physeq)[, "Order"])

# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Order = apply(gentab > observationThreshold, 2, sum)
BB$Order = Order

###################
# Create a factor corresponding to the Class
genfac = factor(tax_table(physeq)[, "Class"])

# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Class = apply(gentab > observationThreshold, 2, sum)
BB$Class = Class

#################
# Create a factor corresponding to the Phylum
genfac = factor(tax_table(physeq)[, "Phylum"])

# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Phylum = apply(gentab > observationThreshold, 2, sum)
BB$Phylum = Phylum

##################
# Create a factor corresponding to the Kingdom
genfac = factor(tax_table(physeq)[, "Kingdom"])

# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(physeq), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})

observationThreshold = 1

Kingdom = apply(gentab > observationThreshold, 2, sum)
BB$Kingdom = Kingdom

write.csv(BB, "7-levels_samples_stat.csv")




#####################
## https://stackoverflow.com/questions/53032504/combine-otu-and-tax-table-and-replace-actual-sequences-with-otu-ids-phyloseq-da


library(dplyr)
library(UpSetR)

# define 0 and 1 function 
binary_fun <- function(x) {ifelse(x > 0, 1, 0)}

# for (i in c("Barley", "Beans", "Bulksoil", "Oats", "OilseedRape", "Sugarbeet", "Wheat")){

(ps = filter_taxa(physeq2, function(x) sum(x) > 1.75, TRUE)) # total sum greater than 1.75%   # Original 247213 

AyBCode <- subset_samples(physeq2, Type== "Barley" ) ### noted that t() may required

AyBCode = merge_samples(AyBCode, "Group")
group_list <- row.names(sample_data(AyBCode))

genfac = factor(tax_table(AyBCode)[, "Genus"])
genfac = as.data.frame(genfac)

# Extract abundance matrix from the phyloseq object
OTU_matrix = as(otu_table(AyBCode), "matrix")
OTU_matrix = as.data.frame(t(OTU_matrix)) # tramspose the matrix

# Add a new column with the row sums
OTU_matrix$row_sum <- rowSums(OTU_matrix)
OTU_matrix = cbind(genfac, OTU_matrix)

df = OTU_matrix %>% 
  group_by(genfac) %>% 
  summarise(across(everything(), sum))


# Subset the data frame to only keep rows where the row sum is 1 or greater
df2 <- df[df$row_sum >= 1, ]
df_binary <- apply(df2[2:10], 2, binary_fun)
# convert resulting matrix back to dataframe
df_binary <- as.data.frame(df_binary)
df_binary$name <- df2$genfac

write.csv(df_binary, "Barley_Genus_group.csv")





## export the data in pdf format
pdf(paste0("UpSetR_", "Wheat", "_Genus1.pdf"), 7, 4.9)
upset(df_binary, sets = group_list, sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")
dev.off()

# }

################# output OTU ####
physeq
head(otu_table(physeq))

library("dplyr")

Y <- otu_table(physeq)
Y = as.data.frame(Y)

Z <- Y  %>%
      replace(is.na(.), 0) %>%
      summarise_all(sum)

Z <- t(Z)
write.csv(Z, "ASV_sum.csv")

##################################
## https://github.com/joey711/phyloseq/issues/613

# Extract abundance matrix from the phyloseq object

B <- cbind(data.frame(otu_table(physeq)), data.frame(tax_table(physeq)))
write.csv(B, "ASV_sum.csv")
