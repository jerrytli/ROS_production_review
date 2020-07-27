#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 4/15/2020
# Script Purpose: Plotting ROS data from divergent cirtus responding to different PAMPS
# Inputs Necessary: 
# Outputs: 
#-----------------------------------------------------------------------------------------------


system("mafft --retree 2 --thread 12 --maxiterate 1000 --localpair NADPH_oxidase_homologs.fasta > 'NADPH_oxidase_homologs_alignment2'")



######################################################################
#library packages need to load
######################################################################

#make sure to set path to the same place where the figure 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Biostrings)
devtools::install_version("latticeExtra", version = "0.6-28")
library(latticeExtra)
library(Hmisc)
library(TraMineR)
library(seqinr)


#########################################################
# funciton - convert DNAstringset attribute to dataframe
#########################################################

aa2df <- function(dss){
  return(data.frame(names = rownames(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}

dss2df <- function(dss){
  return(data.frame(width = width(dss), names = names(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}

######################################################################
#losd file
######################################################################

hold_file <- Biostrings::readAAMultipleAlignment(filepath = file.choose(), format = 'fasta')
holdSeq <-  aa2df(hold_file)
referenceSeq_AtRBOHD <- holdSeq[grepl("AtRBOHD",holdSeq$names),2]



## alternate method
hold_file <- Biostrings::readAAStringSet(filepath = file.choose(), format = "fasta")
holdSeq2 <- dss2df(hold_file)
reference_AtRBOHD <- holdSeq2[grepl("AtRBOHD",holdSeq2$names),3]

hold_data <- data.frame(holdSeq2$names, stringsAsFactors = FALSE)
pb = txtProgressBar(min = 0, max = (max(holdSeq2$width) - window_size), initial = 0, style = 3) 
window_size <- 15


for (j in 1:max(holdSeq2$width)){
  if((j+window_size) <= max(holdSeq2$width)){
    hold_scores <- list()
    for (i in 1:nrow(hold_data)){
      hold_sequence <- holdSeq2[i,3]
      alignment_pair <- pairwiseAlignment(substr(reference_AtRBOHD,j,(j+window_size)),
                                          substr(hold_sequence, j, (j+window_size)))
      alignment_score <- pid(alignment_pair, type = "PID4")
      hold_scores[[i]] <- alignment_score
    }
    hold_data <- cbind(hold_data, unlist(hold_scores))
  }
  setTxtProgressBar(pb,j)
}


for (j in 1:ncol(hold_data)){
  if(j < ncol(hold_data)){
    colnames(hold_data)[j+1] <- j
  }
}

average_data <- as.data.frame(colMeans(hold_data[,100:1000]))


ggplot(average_data, aes(x = as.numeric(rownames(average_data)), y = colMeans(hold_data[, 100:1000]))) +
  geom_point()
