#!/bin/env Rscript

"> plot_venn.GO_bp.R <

Plots a proportional Venn to display shared GO terms between four different
gene groups:
  1. upregulated/methylated
  2. upregulated/unmethylated
  3. downregulated/methylated
  4. downregulated/unmethylated
" -> doc

library(eulerr)

read_summarised_topGO_BP <- function (filename) {
  conn <- file(filename, open='r')
  temp <- readLines(conn)
  temp2 <- c()
  
  # skip the first line (which would be "-- bp_... --"), and stop the loop
  # when the first blank line is encountered (i.e., include BP terms, and 
  # exclude CC and MF)
  for (i in 2:length(temp)) {
    if (temp[i] == '') break
    
    temp2 <- c(temp2, temp[i])
  }
  
  close(conn)
  
  # convert the data into a dataframe
  data <- read.delim(text=temp2)
  
  # then subselect important columns
  data <- data[, c('GO.ID', 'Term', 'P_value')]
  
  data
}

# read summarised data
up_meth <- read_summarised_topGO_BP('topGO_summarised_results/summary_deg32.up_meth.txt')
up_unmeth <- read_summarised_topGO_BP('topGO_summarised_results/summary_deg32.up_unmeth.txt')
down_meth <- read_summarised_topGO_BP('topGO_summarised_results/summary_deg32.down_meth.txt')
down_unmeth <- read_summarised_topGO_BP('topGO_summarised_results/summary_deg32.down_unmeth.txt')

# plot the venn diagram into a pdf file
venn_data <- list(UM=up_meth$GO.ID,
                  UU=up_unmeth$GO.ID, 
                  DM=down_meth$GO.ID,
                  DU=down_unmeth$GO.ID)

#pdf('venn.GO_bp.pdf', width=4, height=3)
set.seed(1)
plot(euler(venn_data),
     fills=c('#fddbc7', '#ef8a62', '#d1e5f0', '#67a9cf'),
     edges=FALSE,
     quantities=TRUE)
#dev.off()

# check what terms overlapped in the venn diagram
merge(up_meth, up_unmeth, by=c('GO.ID', 'Term'), suffixes=c('.UM', '.UU'))

merge(up_meth, down_meth, by=c('GO.ID', 'Term'), suffixes=c('.UM', '.DM'))

merge(down_meth, down_unmeth, by=c('GO.ID', 'Term'), suffixes=c('.DM', '.DU'))
