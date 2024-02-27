#This is the main script to run the analysis of the data
#In this script we will explore and analyse the genomic environment of a given protein
#Have fun exploring!


# Load required libraries ####
library(tidyverse)
library(ape)
library(dplyr)
library(gggenes)
library(ggplot2)

# Sourcing the scripts ####
source('01_open.R')
source('02_clean.R')
source('03_functions.R')

# Loading the Data ####
main <- function(){
  # Read protein and assembly data
  BASEPAIRS <- 300
  MAX_NEIGHBORS <- 15
  
  # Generate protein alias data from files in representative folder
  PROTEIN_ALIAS <- read_represatatives() 

  # Get neighboring proteins
  all.neighbours <- collec.all.neigbour(protein.assembly, BASEPAIRS, MAX_NEIGHBORS)
  
  #check if results are saved
  if(file.exists(paste('output/all_neighbours_bp',BASEPAIRS, '_n', MAX_NEIGHBORS,'.csv', sep = ""))){
    all.neighbours <- read_csv(paste('output/all_neighbours_bp',BASEPAIRS, '_n', MAX_NEIGHBORS,'.csv', sep = ""))
  }
  
  # Plot the neighbors
  if(protein_of_interest==""|is.null(protein_of_interest)|protein_of_interest %in% PROTEIN_ALIAS$alias){
    print('No valid protein of interest giving, skipping plotting')
  }else{
    plot.neighbours(all.neighbours, protein_of_interest)
  }
  
  # Geting Clade information from text files
  clades <- read_clades()
  
  # Getting Cluster Domain information from text files
  cluster_domains <- read_cluster_domain()
  
  # Get the types of neighbours and their counts
  amount_of_neighbours(cluster_domains)
  
  # Re-import the types of neighbours after manual annotation
  annotated_neighbours <- read_csv("output/types_of_neighbours_annotated.csv", show_col_types = FALSE)
  
  # Combine and plot
  df <- combine_and_plot(all.neighbours, cluster_domains, annotated_neighbours, clades)
}