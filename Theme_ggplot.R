#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 4/20/2020
# Script Purpose: Hold information of style of ggplots 
# Inputs Necessary: n/a
# Outputs: n/a
#-----------------------------------------------------------------------------------------------


######################################################################
#library packages need to load
######################################################################

library(ggplot2)
library(extrafont)

######################################################################
#library packages need to load
######################################################################

my_ggplot_theme <- theme_bw() +
  theme(axis.title.x = element_text(size = 12, color = "black", family = "Arial"),
        axis.title.y = element_text(size = 12, color = "black", family = "Arial"),
        axis.text.x = element_text(size = 11,  color = "black", family = "Arial"),
        axis.text.y = element_text(size = 11, color = "black", family = "Arial"),
        panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", size = 0.8))
