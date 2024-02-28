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
  
  #check if results are already present and if not generate them, it will take a long while if the results need to be generated 
  if(file.exists(paste('output/all_neighbours_bp',BASEPAIRS, '_n', MAX_NEIGHBORS,'.csv', sep = ""))){
    all.neighbours <- read_csv(paste('output/all_neighbours_bp',BASEPAIRS, '_n', MAX_NEIGHBORS,'.csv', sep = ""))
  }else{
    # Get neighboring proteins
    all.neighbours <- collec.all.neigbour(protein.assembly, BASEPAIRS, MAX_NEIGHBORS)
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
  combined_df <- combine_and_plot(all.neighbours, cluster_domains, annotated_neighbours, clades)

  # Plot data
  plot_neighbours_per_clade(combined_df)
}

main()