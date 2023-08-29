# https://github.com/EvaYiwenWang/PLSDAbatch

# CRAN
cran.pkgs <- c('pheatmap', 'vegan', 'ruv', 'UpSetR', 'gplots', 'ggplot2', 'gridExtra', 'performance')

for(c in seq_len(length(cran.pkgs))){
  if (!requireNamespace(cran.pkgs[c], quietly = TRUE))
    install.packages(cran.pkgs[c])
}

# Bioconductor
bioc.pkgs <- c('mixOmics', 'sva', 'limma', 'Biobase', 'metagenomeSeq')

for(b in seq_len(length(bioc.pkgs))){
  if (!requireNamespace(bioc.pkgs[b], quietly = TRUE))
    BiocManager::install(bioc.pkgs[b])
}

# install PLSDAbatch
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch", dependencies = T, build_vignettes = T)

library("PLSDAbatch")

