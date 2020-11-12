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

my_ggplot_theme <- theme_linedraw() +
  theme(text=element_text(size = 16,family="Helvetica"),
          axis.title.x = element_text(size = 14, color = "black", family = "Helvetica"),
          axis.title.y = element_text(size = 14, color = "black", family = "Helvetica"),
        axis.text.x = element_text(size = 12,  color = "black", family = "Helvetica"),
        axis.text.y = element_text(size = 14, color = "black", family = "Helvetica"),
        axis.ticks = element_line(size = 0.8),
        axis.ticks.length.y = unit(1.5, "mm"),
        plot.background = element_rect(fill = 'white'),
        panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", size = 0.8))
