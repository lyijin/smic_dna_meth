#!/bin/env Rscript

"> plot_venn.26v29.26v32.R <

Plots a proportional Venn to display shared differentially expressed genes
between 26v29 and 26v32.
" -> doc

library(eulerr)

# read data from sleuth results
twonine <- read.delim('sleuth_results.26Cv29C.tsv', stringsAsFactors=FALSE)
twonine <- twonine[complete.cases(twonine), ]
twonine <- twonine[twonine$qval < 0.05, ]$target_id

threetwo <- read.delim('sleuth_results.26Cv32C.tsv', stringsAsFactors=FALSE)
threetwo <- threetwo[complete.cases(threetwo), ]
threetwo <- threetwo[threetwo$qval < 0.05, ]$target_id

# plot the venn diagram into a pdf file
venn_data <- list(`26 v 29`=twonine,
                  `26 v 32`=threetwo)

pdf('venn.26v29.26v32.pdf', width=4, height=3)
plot(euler(venn_data),
     fills=c('#bdbdbd', '#636363'),
     edges=FALSE,
     quantities=TRUE)
dev.off()
