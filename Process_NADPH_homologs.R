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

#BiocManager::install('Biostrings')
library(Biostrings)
library(seqinr)
library(ggmsa)
library(reshape2)
library(cowplot)
library(ggpubr)
#library(postscriptFonts)
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


######################################################################
# funciton - calculatee similarity scores of full length 
# compared to similarity score of highlighted region
######################################################################


calculate_sim_score <- function(full_sequence_in, blast_N_hit_in, blast_C_hit_in){
  number_of_hits <- sum(1, nrow(blast_N_hit_in), nrow(blast_C_hit_in))
  
  # not similar hits 
  if(number_of_hits == 1){
    seq_global <- Biostrings::pairwiseAlignment(reference_AtRBOHD, full_sequence_in, type = 'global', substitutionMatrix = "BLOSUM62")
    seq_global <- Biostrings::pid(seq_global, "PID1")
    seq_local_N_term <- NA
    seq_local_C_term <- NA
  }
  
  # if an N-temrinal OR C-terminal hit was present
  if(number_of_hits == 2){
    seq_global <- Biostrings::pairwiseAlignment(reference_AtRBOHD, full_sequence_in, type = 'global', substitutionMatrix = "BLOSUM62")
    seq_global <- Biostrings::pid(seq_global, "PID1")
    if(nrow(blast_N_hit_in) == 1 && nrow(blast_C_hit_in) == 0){
      seq_local_N_term <- Biostrings::pairwiseAlignment(AtRBODH_Nterm_highlight_region, blast_N_hit_in$seq, type = 'global', substitutionMatrix = "BLOSUM62")
      seq_local_N_term <- Biostrings::pid(seq_local_N_term, "PID1")
      seq_local_C_term <- NA
    }
    if(nrow(blast_N_hit_in) == 0 && nrow(blast_C_hit_in) == 1){
      seq_local_N_term <- NA
      seq_local_C_term <- Biostrings::pairwiseAlignment(AtRBODH_Cterm_highlight_region, blast_C_hit_in$seq, type = 'global', substitutionMatrix = "BLOSUM62")
      seq_local_C_term <- Biostrings::pid(seq_local_C_term, "PID1")
    }
  }
  
  # if both the N-terminus and the C-terminus has a hit
  if(number_of_hits > 2){
    seq_global <- Biostrings::pairwiseAlignment(reference_AtRBOHD, full_sequence_in, type = 'global', substitutionMatrix = "BLOSUM62")
    seq_local_N_term <- Biostrings::pairwiseAlignment(AtRBODH_Nterm_highlight_region, blast_N_hit_in$seq, type = 'global', substitutionMatrix = "BLOSUM62")
    seq_local_C_term <- Biostrings::pairwiseAlignment(AtRBODH_Cterm_highlight_region, blast_C_hit_in$seq, type = 'global', substitutionMatrix = "BLOSUM62")
    
    seq_global <- Biostrings::pid(seq_global, "PID1")
    seq_local_N_term <- Biostrings::pid(seq_local_N_term, "PID1")
    seq_local_C_term <- Biostrings::pid(seq_local_C_term, "PID1")
  }
  #IGNORE- for debugging orginally
  #print(paste(number_of_hits, nrow(blast_N_hit_in), nrow(blast_C_hit_in)))
  #print(paste(seq_local_N_term,seq_local_C_term))
  
  my_list <- data.frame("Full Length RBOHD" = seq_global, "N-terminus" = seq_local_N_term, "C-terminus" = seq_local_C_term)
  
  return(my_list)
}




######################################################################
#load fasta and hits files file
######################################################################

## alternate method
hold_file <- Biostrings::readAAStringSet(filepath = "./Data_files/NADPH_oxidase_homologs.fasta", format = "fasta") #NADPH_oxidase_homologs.fasta
holdSeq2 <- dss2df(hold_file)
reference_AtRBOHD <- holdSeq2[grepl("AtRBOHD",holdSeq2$names),3]
AtRBODH_Cterm_highlight_region <- substr(reference_AtRBOHD,680,nchar(reference_AtRBOHD))
AtRBODH_Nterm_highlight_region <- substr(reference_AtRBOHD,1,376)


holdBlast_N_hits <- Biostrings::readAAStringSet(filepath = "./Data_files/N-term-hits.fasta", format = 'fasta') #N-term-hits.fasta
holdBlast_N_hits <- dss2df(holdBlast_N_hits)
rownames(holdBlast_N_hits) <- NULL



holdBlast_C_hits <- Biostrings::readAAStringSet(filepath = "./Data_files/C-term-hits.fasta", format = 'fasta') #C-term-hits.fasta
holdBlast_C_hits <- dss2df(holdBlast_C_hits)
rownames(holdBlast_C_hits) <- NULL



######################################################################
# Loop through all proteins and C-terminal qurty hits to calculate 
# similarity score in relationship to its equivalent reference
######################################################################


hold_local_global_Scores <- data.frame("Full Length" = as.numeric(), "N-terminus" = as.numeric(), "C-terminus" = as.numeric())

for (i in 1:nrow(holdSeq2)){
  #find hits from c-termi
  find_blast_N_hit_seq <- holdBlast_N_hits[grepl(holdSeq2[i,2], holdBlast_N_hits$names, fixed = TRUE),]
  find_blast_C_hit_seq <- holdBlast_C_hits[grepl(holdSeq2[i,2], holdBlast_C_hits$names, fixed = TRUE),]
  #print(paste(i,nrow(find_blast_N_hit_seq),nrow(find_blast_C_hit_seq)))
  
  if (nrow(find_blast_N_hit_seq) == 0 && nrow(find_blast_C_hit_seq) == 0){
    #two in original list didn't meet high enough standards after blast search, so we skip their name
    next    
    print("TRUE")
  }
  
  if (nrow(find_blast_N_hit_seq) == 0 && nrow(find_blast_C_hit_seq) == 1){
    hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_N_hit_seq, find_blast_C_hit_seq))
  }
  
  if (nrow(find_blast_N_hit_seq) == 1 && nrow(find_blast_C_hit_seq) == 0){
    hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_N_hit_seq, find_blast_C_hit_seq))
  }
  
  if (nrow(find_blast_N_hit_seq) == 1 && nrow(find_blast_C_hit_seq) == 1 ){
    hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_N_hit_seq, find_blast_C_hit_seq))
  }  

  
  if (nrow(find_blast_C_hit_seq) > 1 || nrow(find_blast_N_hit_seq) > 1){
    for (j in 1:max(nrow(find_blast_N_hit_seq), nrow(find_blast_C_hit_seq))){
      hold_local_global_Scores <- rbind(hold_local_global_Scores, calculate_sim_score(holdSeq2[i,3], find_blast_N_hit_seq[j,], find_blast_C_hit_seq[j,]))
    }
  }
}


######################################################################
# melt data for plotting global comparison to each other
######################################################################

# alter data to plot
colnames(hold_local_global_Scores) <- c("Full Length","N-terminus","C-terminus")
melt_hold_local_global_Scores <- reshape2::melt(hold_local_global_Scores)


# plot data - For Figure 3A
ggplot(melt_hold_local_global_Scores, aes(x = variable, y = value)) +
  geom_boxplot(fill = "grey60", alpha = 0.6, outlier.alpha = 0) +
  geom_jitter(fill = "black", alpha = 0.4, width = 0.2, height = 0) +
  my_ggplot_theme +
  ylab("Percent AA Similarity") +
  # I used 99 instead of 100 to remove AtRBOHD comparison against itself
  ylim(0,99.9) +
  xlab("") +
  theme(axis.ticks.length.x = unit(0, "mm"), 
        axis.text.x.bottom = element_text(vjust = -1, size = 11), 
        axis.text.y.left = element_text(hjust = 0.8))

#export at device demensions of 6.5 x 2.5 inches



#This warning message is expected as ggplot will not plot NA values and AtRBOHD values comparing against itself
#Warning messages:
 # 1: Removed 12 rows containing non-finite values (stat_boxplot). 
 # 2: Removed 12 rows containing missing values (geom_point). 



######################################################################
# funciton - calculatee similarity scores of different kmer style
# fragement of each c-term homologs relative to c-term AtRBOHD
######################################################################



# update which terminus before running: either holdBlast_N_hits or holdBlast_C_hits
calculate_sim_scan_score <- function(sequence_in, window_size, hit_terminus){
  hold_scores <- list()
  for (j in 1:nchar(sequence_in)){
    if((j + window_size) <= nchar(sequence_in)){
      test_seq_window <- substr(sequence_in, j, (j + window_size))
      alignment_match <- pairwiseAlignment(hit_terminus[1,3], test_seq_window, type = 'local-global')
      #print(alignment_match)
      hold_scores[[j]] <- pid(alignment_match)
    }
  }
  hold_scores <- unlist(hold_scores)
  return(hold_scores)
}


######################################################################
# Run through all protein hits
######################################################################



# Define dataframe and calulcate similarity scores for each NADPH oxidase homolog along the N-terminus
hold_data_NTerm <- data.frame("Position" = seq(1:400)) 
  
# update which terminus before running: either holdBlast_N_hits or holdBlast_C_hits
pb = txtProgressBar(min = 0, max = nrow(holdBlast_N_hits), initial = 0, style = 3) 
for (i in 2:nrow(holdBlast_N_hits)){
  save <- data.frame(c(calculate_sim_scan_score(holdBlast_N_hits[i,3], 10, holdBlast_N_hits), c(rep(NA, (nrow(hold_data_NTerm)-length(calculate_sim_scan_score(holdBlast_N_hits[i,3], 10, holdBlast_N_hits)))))))
  colnames(save) <- holdBlast_N_hits[i,2]
  hold_data_NTerm <- cbind(hold_data_NTerm, save)
  setTxtProgressBar(pb, i)
}
  


# Define dataframe and calulcate similarity scores for each NADPH oxidase homolog along the C-terminus 
hold_data_CTerm <- data.frame("Position" = seq(1:300))

# update which terminus before running: either holdBlast_N_hits or holdBlast_C_hits
pb = txtProgressBar(min = 0, max = nrow(holdBlast_C_hits), initial = 0, style = 3) 
for (i in 2:nrow(holdBlast_C_hits)){
  save <- data.frame(c(calculate_sim_scan_score(holdBlast_C_hits[i,3],7, holdBlast_C_hits), c(rep(NA, (nrow(hold_data_v2)-length(calculate_sim_scan_score(holdBlast_C_hits[i,3],7,holdBlast_C_hits)))))))
  colnames(save) <- holdBlast_C_hits[i,2]
  hold_data_CTerm <- cbind(hold_data_CTerm, save)
  setTxtProgressBar(pb, i)
}




######################################################################
# Akter data and plot scanning window analysis
######################################################################

# N-terminus plotting along amino acid position - For Supplemental Figure
hold_data_NTermM <- melt(hold_data_NTerm, id = c("Position"))

ggplot(hold_data_NTermM, aes(x = Position, y = value)) +
  annotate(geom = "rect", xmin = 19, xmax = 28, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  annotate(geom = "rect", xmin = 36, xmax = 42, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  annotate(geom = "rect", xmin = 130, xmax = 136, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  annotate(geom = "rect", xmin = 145, xmax = 151, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  annotate(geom = "rect", xmin = 160, xmax = 166, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  annotate(geom = "rect", xmin = 344, xmax = 340, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
  my_ggplot_theme +
  ylim(0,100) +
  ylab("Percent AA Similarity") +
  xlab("Position Along N-terminus") +
  stat_summary(geom = 'ribbon', fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95), alpha = 0.3, fill = 'red') +
  stat_summary(aes(y = value), fun = mean, colour = "red", geom="line", size = 1) +
  scale_x_continuous(breaks = c(1, 41, 81, 121, 161, 201, 241, 281, 321, 361),
                     labels = c(1, 41, 81, 121, 161, 201, 241, 281, 321, 361),
                     limits = c(0,376), expand = c(0,0)) +
  theme(plot.margin=unit(c(8,10,8,8),"pt"),
        axis.title.x = element_text(vjust = -1),
        axis.text.x = element_text(vjust = -0.2))


#export at device demensions of 6.5 x 2.5 inches



# C-terminus plotting along amino acid position - For Figure 3B
hold_data_CTermM <- melt(hold_data_CTerm, id = c("Position"))
  
ggplot(hold_data_CTermM, aes(x = Position, y = value)) +
    annotate(geom = "rect", xmin = 20, xmax = 25, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
    annotate(geom = "rect", xmin = 177, xmax = 186, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
    annotate(geom = "rect", xmin = 230, xmax = 241, ymin = -Inf, ymax = Inf, fill = "light grey", alpha = 0.55) +
    my_ggplot_theme +
    ylim(0,100) +
    ylab("Percent AA Similarity") +
    xlab("Position Along C-terminus") +
  stat_summary(geom = 'ribbon', fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95), alpha = 0.3, fill = 'red') +
  stat_summary(aes(y = value), fun = mean, colour = "red", geom="line", size = 1) +
    scale_x_continuous(breaks = c(1, 31, 61, 91, 121, 151, 181, 211, 241),
                       labels = c(680, 710,740,770, 800, 830, 860, 890, 920),
                       limits = c(0,242), expand = c(0,0)) +
    theme(plot.margin=unit(c(8,10,8,8),"pt"),
          axis.title.x = element_text(vjust = -1),
          axis.text.x = element_text(vjust = -0.2)) 


#export at device demensions of 5.7 x 2.5 inches


#old assthetics that were later removed
#cutoff <- data.frame( x = c(-Inf, Inf), y = 59, cutoff = factor(52))
#geom_line(aes( x, y, linetype = cutoff ), cutoff)



######################################################################
# align hits of RBOHD C-term region - for creation of weblogos
######################################################################


system("mafft --thread 12 --maxiterate 1000 --localpair N-term-hits.fasta > 'N-term-hits_alignment'")

  
system("mafft --thread 12 --maxiterate 1000 --localpair C-term-hits.fasta > 'C-term-hits_alignment'")


  

######################################################################
# determing the number of cystein residues in N- and C-terminus
######################################################################


N_term_Cysteines <- list()
for (i in 1:nrow(holdBlast_N_hits)){
  N_term_Cysteines[[i]] <- stringr::str_count(holdBlast_N_hits$seq[i], "C")
}
N_term_Cysteines <- cbind(as.data.frame(holdBlast_N_hits$names), as.data.frame(unlist(N_term_Cysteines)))
colnames(N_term_Cysteines) <- c("Homolog Name", "N-terminus")
N_term_Cysteines <- melt(N_term_Cysteines)




C_term_Cysteines <- list()
for (i in 1:nrow(holdBlast_C_hits)){
  C_term_Cysteines[[i]] <- stringr::str_count(holdBlast_C_hits$seq[i], "C")
}
C_term_Cysteines <- cbind(as.data.frame(holdBlast_C_hits$names), as.data.frame(unlist(C_term_Cysteines)))
colnames(C_term_Cysteines) <- c("Homolog Name", " C-terminus")
C_term_Cysteines <- melt(C_term_Cysteines)


test <- rbind(N_term_Cysteines, C_term_Cysteines)
ggplot(test, aes(x = variable, y = value)) +
  geom_boxplot(fill = "grey60", alpha = 0.6, outlier.alpha = 0) +
  geom_jitter(fill = "black", alpha = 0.4, width = 0.25, height = 0.25) +
  #xlim(0,8) +
  ggtitle("Number of Cysteine Residues") +
  my_ggplot_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))



