#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 1/24/2020
# Script Purpose: 
# Inputs Necessary:
# Outputs: 
#-----------------------------------------------------------------------------------------------

##############################################
# To run python on rserver
# load the reticulate package and 
# activate using minicoda
##############################################

reticulate::repl_python()
#conda activate r-reticulate


##############################################
# if you do not have these libraries, 
# use the following compands on bash command line to import
##############################################

# NOTE: Use pip2 if you're working with python2 (note you'll need to make syntatical adjustments)

sudo apt-get install python-pip
sudo pip3 install numpy
sudo pip3 install pandas
sudo pip3 install sys
sudo pip3 install glob
sudo pip3 install biopython

##############################################
# libraries to load
##############################################

import sys
import os
import glob
import numpy
import pandas as pd
import os.path
from Bio.Seq import Seq


##############################################
# check path before importing files
##############################################


#to check relative path:
os.path.realpath() 

#if this doesn't work, try:
os.getcwd()
os.chdir()
#alter path in glob.glob argument as needed


##############################################
# Step 1: Import files to parse - may need to edit path
##############################################

files_to_import = []
files_to_import = glob.glob("/home/danimstevens/Desktop/ROS_production_review/*.txt")

##############################################
# Step 2: Build up DataFrame:
##############################################

df = pd.DataFrame()
for each_file in files_to_import:
    frame = pd.read_csv(each_file, sep= "\t", header=None)
    #frame['filename'] = each_file # add a column with the filename
    df = df.append(frame)

#print (df)



##############################################
# Step 4: filter data by each gene to write out to file
##############################################

#relabel column of dataframe
df.columns = ['Query','Query_Hit','Percent_Identity','E_value','start','stop','Sequence']

#obtain list of genes to filter by
gene_list = df.Query_Hit



##############################################
# Run this version if you want only top hit for each gene - Use this version for protease Blast results!!
############################################## 

#loop through each gene and seperate into dataframes
for each_gene in gene_list:
  each_gene_array = []
  each_gene_df = pd.DataFrame()
  each_gene_df = df[['Query_Hit','Sequence']]
  each_gene_df['Query_Hit'] = ('>' + each_gene_df['Query_Hit'])
  each_gene_df['Sequence'] = (each_gene_df['Sequence'].map(lambda x: x.replace('-','')))
  each_gene_array = each_gene_df.to_numpy(copy=True)
  numpy.savetxt("C-term-hits.fasta", each_gene_array, delimiter = "\n", newline = "\n", fmt = '%s')


# for debugging
#hits_for_gene.loc[:,['Percent_Identity','Query_length','Strain_Name']]
 
##############################################
# Notes
############################################## 
  
# to increase the number of rows to visualize
#pd.set_option('display.max_rows', 500)
  

