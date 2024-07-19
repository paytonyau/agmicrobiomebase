# How-to Guide: Correcting Batch Effects

The UK Crop Microbiome Cryobank integrates genomic (DNA) data with a
cryobank collection of samples for the soil microbiomes of the UK major
crop plant systems. For this project, the microbiomes are from the
rhizosphere (the soil surrounding the crop plant roots) and from bulk
soil (soil outside the rhizosphere). The Cryobank provides a facility
for researchers to source data and samples, including cryo-preserved
microbial material and genomic and metagenomic sequences from different
soil microbiome environments.

### Convert Qiime2 objects to Phyloseq objects

`Qiime2` is a microbial community analysis tool used for sequencing
analysis, while `Phyloseq` is an R package for analyzing high-throughput
sequencing data. The `Qiime2R` package allows conversion of `Qiime2`
data to `Phyloseq` within R. R enhances features through external
packages from sources like CRAN, Bioconductor, and GitHub. After
installation, packages must be loaded into the R session using the
**library()** function.Once data is imported, Phyloseq enables data
manipulation, analysis, and visualization. This conversion leverages R’s
analysis tools like ggplot2 for visualization, dplyr for data
manipulation, and vegan for ecological community analysis.

    # Download phyloseq from Bioconductor
    # if (!require("BiocManager", quietly = TRUE))
    #   install.packages("BiocManager")
    # BiocManager::install("phyloseq")
    library("phyloseq")

    # Download qiime2R from Github
    # if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
    # devtools::install_github("jbisanz/qiime2R")
    library("qiime2R")

    # install.packages("RColorBrewer")
    library("RColorBrewer")

    # install.packages("doParallel")
    library("doParallel")

    # Convert qiime2 results to phyloseq format
    physeq <- qza_to_phyloseq(
      features = "C:/Users/Admin/OneDrive/Desktop/temp/428_228_220_table_silva138-with-phyla-no-mitochondria-no-chloroplast.qza", # table.qza
      taxonomy = "C:/Users/Admin/OneDrive/Desktop/temp/428_228_220_taxonomy_silva138.qza",
      metadata = "C:/Users/Admin/OneDrive/Desktop/temp/16s_meta-table.txt"
      #, tree = "rooted-tree.qza"
    )

    physeq ## confirm the object

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 266974 taxa and 332 samples ]
    ## sample_data() Sample Data:       [ 332 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 266974 taxa by 7 taxonomic ranks ]

### Remove unwanted (failed and controls) samples before the normalisation

Removing unwanted samples before normalisation is a common step in
microbiome data analysis pipelines. In many cases, some samples may fail
during sequencing or quality control, while others may be controls or
blanks that are not of interest. These samples can introduce noise and
bias in downstream analyses if not removed. By removing the unwanted
samples before normalization, the remaining samples can be normalised
based on their true biological variation, allowing for more accurate
comparisons between samples.

    # sample_names(physeq) 
    # rank_names(physeq) # "Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species"

    ## unwanted samples removel (including failed samples)
    physeq.ori <- subset_samples(physeq, Analysis == "Include")

    ## remove object
    rm(physeq)

### Batch effects correction using Constrained Quantile Normalisation (ConQUR)

Constrained Quantile Normalisation (ConQUR,
<https://www.nature.com/articles/s41467-022-33071-9>) is a normalisation
technique used in high-throughput sequencing data, particularly in
microbiome studies. It is a type of **quantile normalisation** approach
that preserves the relative abundances of taxa between samples while
simultaneously removing systematic technical variations that can arise
due to differences in sequencing depth, PCR amplification bias, or other
factors. ConQUR uses kernel density estimation to model the distribution
of taxon abundances across all samples, and then constrains the
normalisation process to maintain the relative position of each taxon
within that distribution.

    # devtools::install_github("wdl2459/ConQuR")
    library(ConQuR)

    # Convert ASV table to a data frame and transpose
    B <- as.data.frame(physeq.ori@otu_table) # taxa
    B <- t(B)
    B <- as.data.frame(B) 

    # Extract batch ID from sample data
    batchid = physeq.ori@sam_data$Plate # batchid

    # Extract covariates
    D = physeq.ori@sam_data[, c('Type', 'Soil', 'Location')] #covar
    summary(D)

    ##           Type    Soil     Location
    ##  Barley     :44   CL: 68   AN:35   
    ##  Beans      :44   CY: 69   BE:33   
    ##  Bulksoil   :45   SC: 70   BO:35   
    ##  Oats       :45   SL:101   BU:35   
    ##  OilseedRape:41            HE:35   
    ##  Sugarbeet  :45            SH:68   
    ##  Wheat      :44            YO:67

    # Correct for batch effects using ConQuR package
    options(warn=-1) # required to call
    taxa_correct1 = ConQuR(tax_tab = B, 
                           batchid = batchid, 
                           covariates = D, 
                           batch_ref="1"
                           ) #  warning messages may appear & it can be ignored

    # Transpose the corrected matrix and convert it to a data frame
    taxa_correct2 <- t(taxa_correct1)
    taxa_correct2 <- as.data.frame(taxa_correct2)

    # Create new ASV table, taxonomy table, and sample data
    ASV = otu_table(taxa_correct2, taxa_are_rows = TRUE)
    TAXA = tax_table(physeq.ori)
    sampledata = sample_data(physeq.ori)

    # repack the objects into a level 4 phyloseq structural data
    physeq.norm = phyloseq(ASV, TAXA, sampledata)

    # Save the "physeq.norm" object for the other case study analysis
    # save(physeq.norm, file = "norm.RData")
    # load("C:/Users/Admin/OneDrive/Desktop/temp/norm.RData")

    # remove
    rm(B, D, batchid, taxa_correct1, taxa_correct2, ASV, TAXA, sampledata, to_skip)

### Bata diversity - before and after the normalisation

Beta diversity quantifies variation in microbial composition among
samples, aiding in identifying patterns in microbial distribution.
Non-Metric Multidimensional Scaling (NMDS) and Principal Coordinates
Analysis (PCoA) are ordination techniques used for beta diversity
analysis.

**NMDS** preserves the rank order of pairwise dissimilarities between
samples in a lower-dimensional space, making it suitable for cases where
distances between samples are not well-preserved. The distances on the
NMDS plot reflect the similarities or dissimilarities between samples
but are not directly interpretable.

**PCoA**, a metric multidimensional scaling technique, attempts to
preserve the actual distances between samples in a lower-dimensional
space. The distances on the PCoA plot reflect the actual dissimilarities
between samples. Unlike NMDS, PCoA may not perform as well with
non-linear or rank-based dissimilarity measures.

Here, we employ NMDS to analyze Beta diversity, allowing us to draw
comparisons between the states before and after normalisation.

    # install.packages("ggplot2")
    library("ggplot2")
    # install.packages("dplyr")
    library("dplyr")
    # install.packages("ggpubr")
    library("ggpubr")

### Bata diversity - before the normalisation

    # method options: NMDS / PCoA
    NMDS1 <- ordinate(physeq = physeq.ori, 
                      method = "NMDS", 
                      distance = "bray"
                      )

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.1860312 
    ## Run 1 stress 0.2021264 
    ## Run 2 stress 0.1900259 
    ## Run 3 stress 0.1953407 
    ## Run 4 stress 0.1907094 
    ## Run 5 stress 0.2128206 
    ## Run 6 stress 0.1945328 
    ## Run 7 stress 0.1907156 
    ## Run 8 stress 0.1932072 
    ## Run 9 stress 0.1915157 
    ## Run 10 stress 0.1877566 
    ## Run 11 stress 0.1883731 
    ## Run 12 stress 0.1856945 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01453486  max resid 0.2225923 
    ## Run 13 stress 0.1870766 
    ## Run 14 stress 0.2008177 
    ## Run 15 stress 0.1861711 
    ## ... Procrustes: rmse 0.009250783  max resid 0.1480255 
    ## Run 16 stress 0.1908439 
    ## Run 17 stress 0.1918221 
    ## Run 18 stress 0.1948028 
    ## Run 19 stress 0.196395 
    ## Run 20 stress 0.1938841 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##      2: no. of iterations >= maxit
    ##      9: stress ratio > sratmax
    ##      9: scale factor of the gradient < sfgrmin

    # Plot ordination
    plot_ordination(physeq = physeq.ori,
                          ordination = NMDS1,
                          color = "Plate",
                          shape = "Type"
                          ) + 
      theme_classic() + 
      geom_point(aes(color = Plate), alpha = 1, size = 3.5) +
      stat_ellipse(aes(color = Plate, group = Plate), geom = "path", alpha = 1.5) + # Add ellipses
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
      scale_color_brewer(palette = "Dark2") + 
      scale_fill_brewer(palette = "Dark2") + 
      scale_shape_manual(values = c(15, 17, 3, 4, 16, 21, 22, 23)) # Set custom shapes

![](case00B-Batch_effects_files/figure-markdown_strict/Beta_1-1.png)

### Bata diversity - after the normalisation

    ##### based on Plate and Type - after the normalisation
    # method options: NMDS / PCoA
    NMDS2 <- ordinate(physeq = physeq.norm, 
                      method = "NMDS", 
                      distance = "bray"
                      )

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.2013146 
    ## Run 1 stress 0.2490011 
    ## Run 2 stress 0.2095548 
    ## Run 3 stress 0.2013147 
    ## ... Procrustes: rmse 9.700457e-06  max resid 0.0001026464 
    ## ... Similar to previous best
    ## Run 4 stress 0.2099039 
    ## Run 5 stress 0.2050828 
    ## Run 6 stress 0.199885 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01324442  max resid 0.1051303 
    ## Run 7 stress 0.198506 
    ## ... New best solution
    ## ... Procrustes: rmse 0.007524845  max resid 0.1039669 
    ## Run 8 stress 0.2226366 
    ## Run 9 stress 0.198506 
    ## ... Procrustes: rmse 1.052236e-05  max resid 9.520903e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.2129426 
    ## Run 11 stress 0.2424379 
    ## Run 12 stress 0.198657 
    ## ... Procrustes: rmse 0.003398883  max resid 0.05779645 
    ## Run 13 stress 0.1987556 
    ## ... Procrustes: rmse 0.004635021  max resid 0.07616554 
    ## Run 14 stress 0.206461 
    ## Run 15 stress 0.2076644 
    ## Run 16 stress 0.2013321 
    ## Run 17 stress 0.1996339 
    ## Run 18 stress 0.2011418 
    ## Run 19 stress 0.1986446 
    ## ... Procrustes: rmse 0.004175996  max resid 0.05802733 
    ## Run 20 stress 0.2128785 
    ## *** Best solution repeated 1 times

    # Plot ordination
    plot_ordination(physeq = physeq.norm,
                    ordination = NMDS2,
                    color = "Plate",
                    shape = "Type"
                    ) +
      theme_classic() + 
      geom_point(aes(color = Plate), alpha = 1, size = 3.5) +
      stat_ellipse(aes(color = Plate, group = Plate), geom = "path", alpha = 1.5) + # Add ellipses
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
      scale_color_brewer(palette = "Dark2") + 
      scale_fill_brewer(palette = "Dark2") + 
      scale_shape_manual(values = c(15, 17, 3, 4, 16, 21, 22, 23)) # Set custom shapes

![](case00B-Batch_effects_files/figure-markdown_strict/Beta_2-1.png)

#### Convert `physeq.norm` object to ASV matrix and Save the “physeq.norm” object for the other case study analysis

### Session Info

    sessionInfo()

    ## R version 4.4.1 (2024-06-14 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 22631)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggpubr_0.6.0       dplyr_1.1.4        ggplot2_3.5.1      ConQuR_2.0        
    ##  [5] doParallel_1.0.17  iterators_1.0.14   foreach_1.5.2      RColorBrewer_1.1-3
    ##  [9] qiime2R_0.99.6     phyloseq_1.48.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] tensorA_0.36.2.1        rstudioapi_0.16.0       jsonlite_1.8.8         
    ##   [4] shape_1.4.6.1           magrittr_2.0.3          NADA_1.6-1.1           
    ##   [7] farver_2.1.2            rmarkdown_2.27          zlibbioc_1.48.2        
    ##  [10] vctrs_0.6.5             multtest_2.58.0         ROCR_1.0-11            
    ##  [13] RCurl_1.98-1.16         base64enc_0.1-3         rstatix_0.7.2          
    ##  [16] htmltools_0.5.8.1       broom_1.0.6             truncnorm_1.0-9        
    ##  [19] Rhdf5lib_1.24.2         Formula_1.2-5           rhdf5_2.46.1           
    ##  [22] KernSmooth_2.23-24      htmlwidgets_1.6.4       plyr_1.8.9             
    ##  [25] igraph_2.0.3            lifecycle_1.0.4         pkgconfig_2.0.3        
    ##  [28] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
    ##  [31] GenomeInfoDbData_1.2.12 clue_0.3-65             digest_0.6.36          
    ##  [34] colorspace_2.1-0        spatial_7.3-17          S4Vectors_0.40.2       
    ##  [37] Hmisc_5.1-3             vegan_2.6-6.1           labeling_0.4.3         
    ##  [40] randomForest_4.7-1.1    fansi_1.0.6             abind_1.4-5            
    ##  [43] httr_1.4.7              mgcv_1.9-1              compiler_4.4.1         
    ##  [46] withr_3.0.0             htmlTable_2.4.2         backports_1.5.0        
    ##  [49] inline_0.3.19           carData_3.0-5           fastDummies_1.7.3      
    ##  [52] highr_0.11              gplots_3.1.3.1          ggsignif_0.6.4         
    ##  [55] MASS_7.3-60.2           bayesm_3.1-6            quantreg_5.98          
    ##  [58] biomformat_1.32.0       caTools_1.18.2          gtools_3.9.5           
    ##  [61] fBasics_4032.96         permute_0.9-7           tools_4.4.1            
    ##  [64] foreign_0.8-86          ape_5.8                 nnet_7.3-19            
    ##  [67] glue_1.7.0              stabledist_0.7-1        nlme_3.1-164           
    ##  [70] rhdf5filters_1.14.1     grid_4.4.1              checkmate_2.3.1        
    ##  [73] cluster_2.1.6           reshape2_1.4.4          ade4_1.7-22            
    ##  [76] generics_0.1.3          gtable_0.3.5            tidyr_1.3.1            
    ##  [79] data.table_1.15.4       car_3.1-2               utf8_1.2.4             
    ##  [82] XVector_0.42.0          rmutil_1.1.10           BiocGenerics_0.50.0    
    ##  [85] ggrepel_0.9.5           pillar_1.9.0            stringr_1.5.1          
    ##  [88] robustbase_0.99-3       splines_4.4.1           lattice_0.22-6         
    ##  [91] survival_3.6-4          GUniFrac_1.8            SparseM_1.84           
    ##  [94] compositions_2.0-8      tidyselect_1.2.1        Biostrings_2.70.3      
    ##  [97] knitr_1.48              gridExtra_2.3           IRanges_2.36.0         
    ## [100] zCompositions_1.5.0-4   stats4_4.4.1            xfun_0.45              
    ## [103] Biobase_2.62.0          statmod_1.5.0           timeDate_4032.109      
    ## [106] matrixStats_1.3.0       DEoptimR_1.1-3          DT_0.33                
    ## [109] stringi_1.8.4           UCSC.utils_1.0.0        yaml_2.3.9             
    ## [112] evaluate_0.24.0         codetools_0.2-20        timeSeries_4032.109    
    ## [115] cqrReg_1.2.1            tibble_3.2.1            cli_3.6.3              
    ## [118] rpart_4.1.23            munsell_0.5.1           Rcpp_1.0.12            
    ## [121] GenomeInfoDb_1.40.1     stable_1.1.6            modeest_2.4.0          
    ## [124] MatrixModels_0.5-3      bitops_1.0-7            glmnet_4.1-8           
    ## [127] scales_1.3.0            statip_0.2.3            purrr_1.0.2            
    ## [130] crayon_1.5.3            rlang_1.1.4
