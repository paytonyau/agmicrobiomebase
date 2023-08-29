# R version 4.2.2
# Convert qiime2 results to phyloseq format
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library("qiime2R")
library("phyloseq")
library("ggplot2")

# Convert qiime2 results to phyloseq format
physeq <- qza_to_phyloseq(
  features = "table_231_228_2_2-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
  # tree = "inst/artifacts/2020.2_moving-pictures/rooted-tree.qza",
  taxonomy = "taxonomy_231_228_2_2.qza",
  metadata = "meta-table_1-4.txt"
)
physeq ## confirm the object

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

# https://wdl2459.github.io/ConQuR/ConQuR.Vignette.html
install.packages("doParallel")
devtools::install_github("wdl2459/ConQuR")
devtools::install_github("wdl2459/ConQuR", build_vignettes = TRUE, force=TRUE)

library(ConQuR)
library(doParallel)

B <- as.data.frame(physeq.a@otu_table) # taxa
B <- t(B)
B <- as.data.frame(B) 

batchid = physeq.a@sam_data$Plate # batchid

D = physeq.a@sam_data[, c('Type', 'Soil', 'Location')] #covar
summary(D)

options(warn=-1) # required to call
taxa_corrected1 = ConQuR(tax_tab = B, batchid = batchid, covariates = D, batch_ref="1")

taxa_corrected2 <- t(taxa_corrected1)
taxa_corrected2 <- as.data.frame(taxa_corrected2)

ASV = otu_table(taxa_corrected2, taxa_are_rows = TRUE)
TAX = tax_table(physeq.a)
sampledata = sample_data(physeq.a)

physeq2 = phyloseq(ASV, TAX, sampledata)

#### Bata diversity
moth_sub_pcoa <- ordinate(physeq = physeq2, method = "NMDS", distance = "bray")

pdf("Beta diversity (Soil_Location) norm_ConQuR_1_1.pdf",6 ,5)
plot_ordination(
  physeq = physeq2,
  ordination = moth_sub_pcoa,
  # title = "NMDS",
  color = "Location",
  shape = "Soil"
) +
  # scale_x_discrete(name ="NMDS1 ()") + 
  # scale_y_discrete(name ="NMDS2 ()") + 
  theme_classic() + 
  geom_point(aes(color = Location), alpha = 1, size = 3) +
  theme(text = element_text(size=18, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.1),
        axis.line = element_line(colour = 'black', size = 1.1),
        axis.text.x = element_text(colour = "black", angle=0, 
                                   hjust=0.5, size = 13, face="bold"),
        axis.text.y = element_text(colour = "black", angle=0, 
                                   hjust=0.5, size = 13, face="bold"),
        axis.title.y = element_text(color="black", size=20,face="bold"), 
        axis.title.x = element_text(color="black", size=20,face="bold")) +
  # stat_ellipse(geom = "polygon", type="norm", alpha=0.25, aes(fill = Group)) + # polygon, path, point
  scale_color_brewer(palette="Dark2") + 
  scale_fill_brewer(palette="Dark2") + 
  scale_shape_manual(values = c(15, 16, 17, 4, 5, 6, 12))
dev.off()

###### Alpha diversity ######
# available measurements [c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")]
# original
tab2 = cbind(x = sample_data(physeq2), 
            y = estimate_richness(physeq2, measures = 'Shannon'))

stat.test <- tab %>%
  # group_by(Neutrophils, GROUP1) %>%
  t_test(Shannon ~ x.Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

pdf("Alpha diversity - Shannon (Location) norm_1_1.pdf", 4.5, 5.5)
ggplot(data = tab2, aes(x = x.Location, y = Shannon, color = x.Location, fill = x.Location)) + 
  theme_classic() + 
  labs(# title = "IBD Patients", 
    x = element_blank(), 
    y = "Alpha Diversity (Shannon)") + 
  geom_point(size = 1.75) + 
  geom_boxplot(alpha = 0.5) + 
  #  stat_pvalue_manual(stat.test, 
  #                   y.position = c(9, 9.5, 10.5, 8.5, 12, 9),
  #                     label = "p.adj.signif",
  #                     face="bold", 
  #                     size = 6, 
  #                     linetype = 1,
  #                     tip.length = 0.02,
  #                     inherit.aes = FALSE) + 
  scale_y_continuous(limits=c(0, 12), breaks = c(0, 2.5, 5, 7.5, 10)) +
  theme(text = element_text(size=18, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.1),
        axis.line = element_line(colour = 'black', size = 1.1),
        axis.text.x = element_text(colour = "black",
                                   angle=270, 
                                   size = 13, face="bold"),
        axis.text.y = element_text(angle=0, hjust=0, colour = "black",
                                   size = 13, face="bold"),
        axis.title.y = element_text(color="black", size=15,face="bold"),
        legend.position = "none") +
  scale_color_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")
dev.off()

## https://github.com/joey711/phyloseq/issues/1487
## https://github.com/joey711/phyloseq/issues/1516
## Calculate the most abundance 

# Agglomerate taxa at genus level
GP.genus <- tax_glom(physeq2, taxrank = "Genus")
physeq_Barley <- subset_samples(GP.genus, Type=="Barley") ### noted that t() may required
physeq_Barley = merge_samples(physeq_Barley, "Group")
# Calculate relative abundance
GP.genus.prop <- transform_sample_counts(physeq_Barley, function(x) x / sum(x))
# Define group/condition type
type = levels(sample_data(GP.genus.prop))

# Retrieve sample name of the defined type
sn = sample_names(sample_data(GP.genus.prop)[sample_data(GP.genus.prop)$Group])

# Retrieve top genera names
top5 = names(head(sort(rowSums(otu_table(GP.genus.prop)), decreasing = TRUE), 9)) 

# Remove unselected samples and taxa
rb = prune_taxa(top5, GP.genus.prop)

# Build data.frame
df = cbind(data.frame(ID = c(taxa_names(rb), "Other"),
                      Genus = c(tax_table(rb)[,"Genus"], "Other")),
           rbind(otu_table(rb), 1 - colSums(otu_table(rb))))

# Order data.frame by genera abundance, "Other" placed at the last
df = df[match(c(top5, "Other"), df$ID),]
library(tidyverse)
df = df %>% mutate_at(vars(Genus), factor) 

# Plot
pivot_longer(df, cols = -c(ID, Genus), names_to = "Sample", values_to = "Abundance") %>%
  ggplot(aes(Sample,  Abundance, fill = Genus)) + geom_col(color = "black") +
  # scale_fill_brewer(palette = "Dark2") + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = type, y = "Average relative abundance")

#######

TOPS = names(sort(taxa_sums(physeq_Barley), TRUE)[1:30])
TOPS
TOPS1 = prune_taxa(TOPS, physeq_Barley)
TOPS1

plot_bar(TOPS1, fill = "Soil", x = "Type", facet_grid=~Alimentacion)
plot_bar(TOPS1, fill = "Genus", x = "Group")
plot_tree(TOPS1)


## ConQuR with the penalized fitting strategy: 
# logistic LASSO regression and quantile LASSO regression
options(warn=-1)
taxa_corrected2 = ConQuR(tax_tab=B, batchid=batchid, covariates=D, batch_ref="1",
                         logistic_lasso=T, quantile_type="lasso", interplt=T)
