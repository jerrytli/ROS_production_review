#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 7/27/2020
# Script Purpose: Processing and Plotting similarity scores for NADPH oxidase homologs
# Inputs Necessary: Theme_ggplot.R, NADPH_oxidase_homologs.fasta, C-term-hits.fasta
# Outputs: 
#-----------------------------------------------------------------------------------------------


######################################################################
# library packages need to load
######################################################################

# make sure to set path to the same place where the figure 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Biostrings)
library(seqinr)
library(ggmsa)
library(reshape2)
library(cowplot)
library(ggpubr)
library(ggstatsplot)
library(postscriptFonts)
library(ggseqlogo)
library(extrafont)
loadfonts(device = "pdf", quiet = FALSE)

##############################################
# Load GGplot style
##############################################


source("./Theme_ggplot.R")


#########################################################
# funciton - convert AAstringset attribute to dataframe
#########################################################

# turning AAmultiplesequence alignment into dataframe
aa2df <- function(dss){
  return(data.frame(names = rownames(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}

# turning AAsequences (fasta) into dataframe
dss2df <- function(dss){
  return(data.frame(width = width(dss), names = names(dss), seq = as.character(dss), stringsAsFactors = FALSE))
}


#########################################################
# funciton - calculatee similarity scores of full length 
# compared to similarity score of highlighted region
#########################################################


calculate_sim_score <- function(full_sequence_in, blast_hit_in){
  seq_global <- Biostrings::pairwiseAlignment(reference_AtRBOHD, full_sequence_in, type = 'global', substitutionMatrix = "BLOSUM62")
  seq_local <- Biostrings::pairwiseAlignment(AtRBODH_Cterm_highlight_region, blast_hit_in, type = 'global', substitutionMatrix = "BLOSUM62")
  my_list <- data.frame("Full Length RBOHD" = pid(seq_global, "PID1"), "RBOHD C-terminus" = pid(seq_local, "PID1"))
  return(my_list)
}




######################################################################
#losd file
######################################################################

## alternate method
hold_file <- Biostrings::readAAStringSet(filepath = file.choose(), format = "fasta") #NADPH_oxidase_homologs.fasta
holdSeq2 <- dss2df(hold_file)
reference_AtRBOHD <- holdSeq2[grepl("AtRBOHD",holdSeq2$names),3]
AtRBODH_Cterm_highlight_region <- substr(reference_AtRBOHD,680,1000)


holdBlast_hits <- Biostrings::readAAStringSet(filepath = file.choose(), format = 'fasta') #c-term-hits.fasta
holdBlast_hits <- dss2df(holdBlast_hits)
rownames(holdBlast_hits) <- NULL

###### need to FIX BUG!!!!! - fixed :)

hold_local_global_Scores <- data.frame("Full Length" = as.numeric(), "C-terminus" = as.numeric())
for (i in 1:nrow(holdSeq2)){
  find_blast_hit_seq <- holdBlast_hits[grepl(holdSeq2[i,2], holdBlast_hits$names, fixed = TRUE),]
  if (nrow(find_blast_hit_seq) == 0){
    #two in original list didn't meet high enough standards after blast search, so we skip their name
    next    
  }
  if (nrow(find_blast_hit_seq) == 1){
    hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_hit_seq[,3]))
  }
  if (nrow(find_blast_hit_seq) > 1){
    for (j in 1:nrow(find_blast_hit_seq)){
      hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_hit_seq[j,3]))
    }
  }
}


# alter data to plot
colnames(hold_local_global_Scores) <- c("Full Length","C-terminus")
melt_hold_local_global_Scores <- melt(hold_local_global_Scores)


# plot data
ggplot(melt_hold_local_global_Scores, aes(x = variable, y = value)) +
  geom_boxplot(fill = "grey60", alpha = 0.6, outlier.alpha = 0) +
  geom_jitter(fill = "black", alpha = 0.4, width = 0.2) +
  my_ggplot_theme +
  ylab("Percent AA Similarity") +
  ylim(0,99) +
  xlab("")




######################################################################
# funciton - calculatee similarity scores of different kmer style
# fragement of each c-term homologs relative to c-term AtRBOHD
######################################################################


calculate_sim_scan_score <- function(sequence_in, window_size){
  hold_scores <- list()
  for (j in 1:nchar(sequence_in)){
    if((j + window_size) <= nchar(sequence_in)){
      test_seq_window <- substr(sequence_in, j, (j + window_size))
      alignment_match <- pairwiseAlignment(holdBlast_hits[1,3], test_seq_window, type = 'local-global')
      hold_scores[[j]] <- pid(alignment_match)
    }
  }
  hold_scores <- unlist(hold_scores)
  return(hold_scores)
}


######################################################################
# Run through all protein hits
######################################################################

# Define dataframe and calulcate similarity scores for each NADPH oxidase homolog 
hold_data_v2 <- data.frame("Position" = seq(1:300))
pb = txtProgressBar(min = 0, max = nrow(holdBlast_hits), initial = 0, style = 3) 
for (i in 2:nrow(holdBlast_hits)){
  save <- data.frame(c(calculate_sim_scan_score_v2(holdBlast_hits[i,3],20), 
                       c(rep(NA, 
                             (nrow(hold_data_v2)-length(calculate_sim_scan_score_v2(holdBlast_hits[i,3],20)))))))
  colnames(save) <- holdBlast_hits[i,2]
  hold_data_v2 <- cbind(hold_data_v2, save)
  setTxtProgressBar(pb, i)
}



######################################################################
# Run through all protein hits
######################################################################


hold_data_melt2 <- melt(hold_data_v2, id = c("Position"))

cutoff <- data.frame( x = c(-Inf, Inf), y = 59, cutoff = factor(59) )
cutoff2 <- data.frame( x = c(-Inf, Inf), y = 52, cutoff2 = factor(52) )
ggplot(hold_data_melt2, aes(x = Position, y = value)) +
  annotate(geom = "rect", xmin = 20, xmax = 28, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.7) +
  annotate(geom = "rect", xmin = 179, xmax = 187, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.7) +
  annotate(geom = "rect", xmin = 229, xmax = 237, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.7) +
  my_ggplot_theme +
  ylim(0,100) +
  ylab("Percent AA Similarity") +
  xlab("Position Along C-terminus") +
  stat_summary(aes(y= value), fun = mean, colour="red", geom="line", size = 1) +
  scale_x_continuous(breaks = c(1, 31, 61, 91, 121, 151, 181, 211, 241),
                     labels = c(680, 710,740,770, 800, 830, 860, 890, 920),
                     limits = c(0,242), expand = c(0,0)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),
        axis.title.x = element_text(vjust = -1),
        axis.text.x = element_text(vjust = -0.2))


######################################################################
# align hits of RBOHD C-term region - for creation of weblogos
######################################################################

  
system("mafft --thread 12 --maxiterate 1000 --localpair C-term-hits.fasta > 'C-term-hits_alignment'")


  
######################################################################
# weblogos
######################################################################  
  
# after much discussion, we decided to use bardo's weblogos instead. however, I'm choosing
# to leave this code base here for my own refernce
  
x <- readAAStringSet(filepath = "./NADPH-oxidase-cTerm-reference.fasta", format = "fasta")

  
S703_file <- readLines(file.choose())
S862_file <- readLines(file.choose())
T912_file <- readLines(file.choose())

  
# S862 
S703_string <- ggmsa::ggmsa(x, seq_name = FALSE, char_width = 0.6, color = "Chemistry_AA", 20, 28, posHighligthed = c(24)) + 
  my_ggplot_theme + 
  theme(panel.border = element_rect(color = "black", size = 0),
        legend.position = "none", axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin=unit(c(-1,1,1,1), "cm")) + 
  scale_x_continuous(breaks = c(20,21,22,23,24,25,26,27,28),
                     labels = c(699,700,701,702,703,704,705,706,708))


S703_logo <- ggseqlogo(S703_file,  seq_type='aa', method = 'prob') + 
  my_ggplot_theme +
  theme(panel.border = element_rect(color = "black", size = 0.5),
        legend.position = "none", axis.text.x = element_blank(),
        plot.margin=unit(c(1,1,-1,1), "cm"))

plot_grid(S703_logo, S703_string,  ncol = 1, align = 'v')


# S862
S862_string <- ggmsa::ggmsa(x, seq_name = FALSE, char_width = 0.6, color = "Chemistry_AA", 179, 187, posHighligthed = c(183)) + 
  my_ggplot_theme + 
  theme(panel.border = element_rect(color = "black", size = 0),
        legend.position = "none", axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin=unit(c(-1,1,1,1), "cm")) + 
  scale_x_continuous(breaks = c(179,180,181,182,183,184,185,186,187),
                     labels = c(858,859,860,861,862,863,864,865,866))



S862_logo <- ggseqlogo(S862_file,  seq_type='aa',method = 'prob') + 
  my_ggplot_theme +
  theme(panel.border = element_rect(color = "black", size = 0.5),
        legend.position = "none", axis.text.x = element_blank(),
        plot.margin=unit(c(1,1,-1,1), "cm"))

plot_grid(S862_logo, S862_string,  ncol = 1, align = 'v')



# T912
T912_string <- ggmsa::ggmsa(x, seq_name = FALSE, char_width = 0.6, color = "Chemistry_AA", 229, 237, posHighligthed = c(233)) + 
  my_ggplot_theme + 
  theme(panel.border = element_rect(color = "black", size = 0),
        legend.position = "none", axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.margin=unit(c(-1,1,1,1), "cm")) + 
  scale_x_continuous(breaks = c(229,230,231,232,233,234,235,236,237),
                     labels = c(908,909,910,911,912,913,914,915,916))



T912_logo <- ggseqlogo(T912_file,  seq_type='aa') + 
  my_ggplot_theme +
  theme(panel.border = element_rect(color = "black", size = 0.5),
        legend.position = "none", axis.text.x = element_blank(),
        plot.margin=unit(c(1,1,-1,1), "cm"))

plot_grid(T912_logo, T912_string,  ncol = 1, align = 'v')

