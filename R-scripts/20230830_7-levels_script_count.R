# https://github.com/joey711/phyloseq/issues/337#issuecomment-42254256

# R version 4.2.0
### Convert qiime2 results to phyloseq format
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(qiime2R)

# Convert qiime2 results to phyloseq format
physeq <- qza_to_phyloseq(
  features = "16s/GreenGenes_13_8/428_228_220_table_gg-13-8-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
  # tree = "inst/artifacts/2020.2_moving-pictures/rooted-tree.qza",
  taxonomy = "16s/GreenGenes_13_8/428_228_220_taxonomy_gg-13-8.qza",
  metadata = "16s/meta-table.txt"
)
physeq ## confirm the object

## Subset the data
sample.remove <- c("SU.SC.SH.4", "SU.CY.YO.5", 
                   "BE.SC.SH.5", "BE.SL.AN.3",
                   "OA.SL.SH.4","OA.SL.BE.2",
                   "OR.SL.SH.1", "OR.SL.SH.2", 
                   "WH.SC.SH.4",
                   "BA.SL.BE.3.RE", "CO.CL.YO.5.RE", "OA.SL.AN.4.RE", 
                   "CO.CY.YO.5.RE2", #
                   "p.ve.1", "p.ve.2","p.ve.3", "p.ve.4", "p.ve.5",
                   "n.ve.ext.1", "n.ve.ext.2", "n.ve.ITS.1",
                   "n.ve.1","n.ve.2","n.ve.3","n.ve.4", "n.ve.5")
physeq.a <- subset_samples(physeq,  !(id %in% sample.remove)) # Durie data

# Visualise data
# Script adapted from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

library("phyloseq")
library("ggplot2")
library("ggpubr")

sample_names(physeq)
rank_names(physeq) # "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species"

# Script adopted from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html


##############################
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

########################## in Crop types #########
AyBCode <- merge_samples(physeq.a, "Group", fun = sum)
physeq = AyBCode
##############################

gentab_levels <- list()
observationThreshold <- 1

for (level in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  genfac <- factor(tax_table(physeq)[, level])
  gentab <- apply(otu_table(physeq), MARGIN = 1, function(x) {
    tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
  })
  
  level_counts <- apply(gentab > observationThreshold, 2, sum)
  BB <- as.data.frame(level_counts)
  BB$name <- row.names(BB)
  # Add the result to the gentab_levels list
  gentab_levels[[level]] <- BB
}

require(tidyverse);
B <- reduce(gentab_levels, full_join, by = "name");

# Write the data frame to a file
write.csv(B, file = "level_counts_by_group.csv", row.names = FALSE)


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
