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
library(postscriptFonts)
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
  ylab("Percent Similarity\n") +
  ylim(0,98) +
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

#### Def
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

ggplot(hold_data_melt2, aes(x = Position, y = value)) +
  geom_rect(aes(xmin = 20, xmax = 28, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
  geom_rect(aes(xmin = 179, xmax = 187, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
  geom_rect(aes(xmin = 229, xmax = 237, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
  my_ggplot_theme +
  ylim(0,100) +
  ylab("Percent Similarity\n") +
  xlab("Position") +
  stat_summary(aes(y= value), fun = mean, colour="red", geom="line", size = 1) +
  scale_x_continuous(breaks = c(0, 30, 60, 90, 120, 150, 180, 210, 240),
                     labels = c(680, 710,740,770, 800, 830, 860, 890, 920),
                     limits = c(0,240), expand = c(0,0)) 



######################################################################
# align hits of RBOHD C-term region
######################################################################

system("mafft --thread 12 --maxiterate 1000 --localpair C-term-hits-fixed.fasta > 'C-term-hits_alignment'")


gmsa(Cm_ncppb382_MSA_cps22, seq_name = FALSE, char_width = 0.6, color = "Chemistry_AA", posHighligthed = c(2,10,18,19,21,22)) + 
  my_ggplot_theme + 
  theme(panel.border = element_rect(color = "black", size = 0.2),
        legend.position = "none") +
  scale_y_discrete(limits = rev(rownames(Cm_ncppb382_MSA_cps22)))


############################################################################################################

gmsa(substr(), seq_name = TRUE, char_width = 0.6, color = "Chemistry_AA", posHighligthed = c(2,10,18,19,21,22)) + 
  my_ggplot_theme + 
  theme(panel.border = element_rect(color = "black", size = 0.2),
        legend.position = "none") +
  scale_y_discrete(limits = rev(rownames(Cm_ncppb382_MSA_cps22)))




######################################################################
# test reverse scanning windpw alignment
######################################################################

calculate_sim_scan_score_v2_rev <- function(sequence_in, window_size){
  hold_scores <- list()
  for (j in nchar(sequence_in):1){
    if((j - window_size) >= 1){
      test_seq_window <- substr(sequence_in, (j - window_size), j)
      alignment_match <- pairwiseAlignment(holdBlast_hits[1,3], test_seq_window, type = 'local-global')
      print(paste(refernce_seq_window, test_seq_window, sep = "_"))
            print(alignment_match)
      print(pid(alignment_match))
      hold_scores[[j]] <- pid(alignment_match)
    }
  }
  hold_scores <- unlist(hold_scores)
  hold_scores <- rev(hold_scores)
  return(hold_scores)
}



####
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
#losd file
######################################################################

#alignment_of_RBOHD_cterm <- Biostrings::readAAMultipleAlignment(filepath = file.choose(), format = "fasta")
#alignment_of_RBOHD_cterm <- aa2df(alignment_of_RBOHD_cterm)
#rownames(alignment_of_RBOHD_cterm) <- NULL
#hold_data <- data.frame(holdSeq2$names, stringsAsFactors = FALSE)


#hold_data <- data.frame(calculate_sim_scan_score(holdBlast_hits[2,3],20))
#colnames(hold_data) <- c(holdBlast_hits[2,2])

#hold_data <- cbind(as.numeric(rownames(hold_data)), hold_data)
#colnames(hold_data) <- c("position",holdBlast_hits[2,2])

hold_data_v3 <- calculate_sim_scan_score_v2_rev(holdBlast_hits[2,3], 20)
pb = txtProgressBar(min = 0, max = nrow(holdBlast_hits), initial = 0, style = 3) 
for (i in 3:nrow(holdBlast_hits)){
  save <- calculate_sim_scan_score_v2_rev(holdBlast_hits[i,3], 20)
  hold_data_v3 <- rbind(hold_data_v3, save)
  #save <- data.frame(calculate_sim_scan_score(holdBlast_hits[i,3],20))
  #colnames(save) <- holdBlast_hits[i,2]
  #hold_data <- cbind(hold_data, save)
  setTxtProgressBar(pb, i)
}





g



#colnames(hold_data_v3)[1] <- c("Position")
#hold_data_melt_v3 <- melt(hold_data_v3, id = c("Position"))
#ggplot(hold_data_v3, aes(x = start.point, y = score)) +
# geom_rect(aes(xmin = 21, xmax = 27, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
#  geom_rect(aes(xmin = 180, xmax = 186, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
#  geom_rect(aes(xmin = 230, xmax = 236, ymin = -Inf, ymax = Inf), fill = "light grey", alpha = 0.01) +
#  my_ggplot_theme +
# ylim(0,100) +
#xlim(0,250) +
#stat_summary(aes(y = score), fun = mean, colour="blue", geom="line", size = 1)


