### Convert qiime2 results to phyloseq format
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(qiime2R)

# Convert qiime2 results to phyloseq format
physeq <- qza_to_phyloseq(
  features = "table_231_228_2_2-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
  # tree = "inst/artifacts/2020.2_moving-pictures/rooted-tree.qza",
  taxonomy = "taxonomy_231_228_2_2.qza",
  metadata = "meta-table_1-4.txt"
)
physeq ## confirm the object

# Visualise data
# Script adapted from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html

library("phyloseq")
library("ggplot2")
library("ggpubr")

sample_names(physeq)
rank_names(physeq) # "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species"

## Subset the data
sample.remove <- c("SU-SC-SH-4", "SU-CY-YO-5", "BE-SC-SH-5", "BE-SL-AN-3", "OA-SL-SH-4","OA-SL-BE-2","WH-SC-SH-4",
                   "BA-SL-BE-3-RE", "CO-CL-YO-5-RE", "OA-SL-AN-4-RE", 
                   "OR-SL-SH-1", "OR-SL-SH-2", # failed samples from plate 4
                   "p-ve-1", "p-ve-2","p-ve-3", "p-ve-4", 
                   "n-ve-1","n-ve-2","n-ve-3","n-ve-4",
                   "n-ve-ext-1", "n-ve-ext-2", "n-ve-ITS-1"
)


physeq.a <- subset_samples(physeq,  !(id %in% sample.remove)) # Durie data

# https://github.com/joey711/phyloseq/issues/960


ps_with_only_genus <- subset_samples(physeq.a, Type=="OilseedRape")
ps_with_only_genus <- subset_taxa(ps_with_only_genus, Class == "Actinobacteria")
taxdf = data.frame(tax_table(ps_with_only_genus))
plot_bar(ps_with_only_genus)
