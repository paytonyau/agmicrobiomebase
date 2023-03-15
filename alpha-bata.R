# R version 4.1.3
# Convert qiime2 results to phyloseq format
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

View(tax_table(physeq))
# Script adapted from https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
# Visualise data
library("phyloseq")
library("ggplot2")
library("ggpubr")

sample_names(physeq)
rank_names(physeq) # "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species"

## Subset the data
sample.remove <- c("SU-SC-SH-4", "SU-CY-YO-5", "BE-SC-SH-5", "BE-SL-AN-3","OA-SL-SH-4","OA-SL-BE-2","WH-SC-SH-4",
                   "BA-SL-BE-3-RE", "CO-CL-YO-5-RE", "OA-SL-AN-4-RE", 
                   "p-ve-1", "p-ve-2","p-ve-3", "p-ve-4", 
                   "n-ve-ext-1", "n-ve-ext-2", "n-ve-ITS-1",
                   "OR-SL-SH-1", "OR-SL-SH-2", 
                   "n-ve-1","n-ve-2","n-ve-3","n-ve-4")
physeq.a <- subset_samples(physeq,  !(id %in% sample.remove)) # Durie data

# https://www.nicholas-ollberding.com/post/introduction-to-phyloseq/
B <- sample_sums(physeq)
B <- as.data.frame(B)
write.csv(B, "processed_total_reads_physeq.csv")

## Normalised number of reads in each sample using median sequencing depth.
total = median(sample_sums(physeq.a))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq.a.norm = transform_sample_counts(physeq.a, standf)

plot_bar(physeq.a.norm, fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#### Bata diversity
moth_sub_pcoa <- ordinate(physeq = physeq.a, method = "NMDS", distance = "bray")

pdf("Beta diversity (Soil_Location).pdf",6 ,5)
plot_ordination(
  physeq = physeq.a,
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


## Function (Modified) ##
beta_boxplot <- function(physeq, method = "bray", group) {
  
  # physeq: phyloseq-class object
  # method: beta-diversity metric. Default "bray", i.e., Bray-Curtis dissimilarity 
  # group: factorial variable to group
  
  ## Packages
  require("phyloseq") # v.1.30.0
  require("ggplot2") # v.3.3.2
  
  ## Identify the correspondence: group and samples
  group2samp <- list() # list to save the correspondence between group <--> samples
  group_list <- get_variable(sample_data(physeq), group) # list of group elements
  for (groups in levels(group_list)) { # loop over the no. of group levels
    target_group <- which(group_list == groups) # vct pos of the curr group variable 
    group2samp[[ groups ]] <- sample_names(physeq)[target_group] # matching samples: based on vct pos
  }  
  
  ## Calculate beta-diversity
  beta_div_dist <- distance(physeq = physeq, method = method)
  beta_div_dist <- as(beta_div_dist, "matrix")
  
  ## Coerce distance mtx into a tidy data frame by group
  dist_df <- data.frame() # save results in df 
  counter <- 1 
  for (groups in names(group2samp)) { # loop over group fct levels 
    sub_dist <- beta_div_dist[ group2samp[[groups]], group2samp[[groups]] ] # subset dist mtx to curr samples
    #print(sub_dist)
    no_samp_col <- ncol(sub_dist) # n cols: curr sub dist
    no_samp_row <- nrow(sub_dist) # n rows: curr sub dist
    for ( cols in seq(no_samp_col) ) { # loop over cols: curr sub_dist
      if ( cols > 1 ) {
        for ( rows in seq((cols-1)) ) { # loop over rows: curr sub_dist 
          ## Save results
          dist_df[ counter, "sample_pair" ] <- paste0( colnames(sub_dist)[cols], "-",  
                                                       rownames(sub_dist)[rows] ) # sample pair
          dist_df[ counter, "group" ] <- groups # group  
          dist_df[ counter, "beta_div_method" ] <- method # method
          dist_df[ counter, "beta_div_value" ] <- sub_dist[rows, cols] # beta-diversity for the sample pair     
          counter = counter + 1
        }
      }
    }
  }
  return(dist_df)
}

## Run the Function ##
beta_boxplot_result <- beta_boxplot(physeq = physeq.b.norm, method = "bray", group = "Treatment") 

## P-value calculation 
stat.test = as.data.frame(beta_boxplot_result)

library(rstatix)
stat.test <- stat.test %>%
  # group_by(Neutrophils, GROUP1) %>%
  t_test(beta_div_value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

### Plotting
pdf("Beta diversity - Box Plot (Normalised).pdf", 7 ,5.5)
ggplot(data = beta_boxplot_result, aes(x = group, y = beta_div_value, color = group, fill = group)) + 
  geom_violin(trim = FALSE, alpha = 0.3) + 
  geom_boxplot(width = 0.07, alpha = 0.5) +
  geom_jitter() + 
  # theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))+ 
  labs(# title = "IBD Patients", 
    x = element_blank(), 
    y = "Distance (Bray)") + 
  theme_classic() +
  # scale_x_discrete(limits=c("Control", "Plunge", "Rate")) +
  stat_pvalue_manual(stat.test, 
                     y.position = 0.8,
                     label = "p.adj.signif",
                     face="bold", 
                     size = 6, 
                     linetype = 1,
                     tip.length = 0.02,
                     inherit.aes = FALSE) + 
  theme(text = element_text(size=18, colour = "black"), 
        axis.ticks = element_line(colour = "black", size = 1.1),
        axis.line = element_line(colour = 'black', size = 1.1),
        axis.text.x = element_text(colour = "black", angle=0, hjust=0.5, size = 13, face="bold"),
        axis.text.y = element_text(colour = "black", angle=0, hjust=0.5, size = 13, face="bold"),
        axis.title.x = element_text(color="black", size=15,face="bold"),
        axis.title.y = element_text(color="black", size=15,face="bold")) +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
dev.off()


###### Alpha diversity ######
# available measurements [c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")]
# original
tab = cbind(x = sample_data(physeq.a), 
            y = estimate_richness(physeq.a, measures = 'Shannon'))

stat.test <- tab %>%
  # group_by(Neutrophils, GROUP1) %>%
  t_test(Shannon ~ x.Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

pdf("Alpha diversity - Shannon (Type).pdf", 4.5, 5.5)
ggplot(data = tab, aes(x = x.Type, y = Shannon, color = x.Type, fill = x.Type)) + 
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
