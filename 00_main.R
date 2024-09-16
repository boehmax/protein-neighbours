#This is the main script to run the analysis of the data
#In this script we will explore and analyse the genomic environment of a given protein
#Have fun exploring!


# Load required libraries ####
library(tidyverse)
library(ape)
library(dplyr)
library(gggenes)
library(ggplot2)
library(RColorBrewer)

# Sourcing the scripts ####
source('R/01_open.R')
source('R/02_clean.R')
source('R/03_function.R')

# Loading the Data ####
main <- function(BASEPAIRS = 300, MAX_NEIGHBORS = 15, PATH = 'data', date = NULL) {
  current_date <- format(Sys.Date(), "%Y-%m-%d") 
  # Ask for Running, reloading or exiting while output file path is created
  output_file <- create_or_validate_output_file_path(basepairs=BASEPAIRS, max_neighbors=MAX_NEIGHBORS, date=date)
  
  # Read protein and assembly data
  protein_assembly_data <- read_protein_assembly_data(PATH = PATH)

  # Access the variables
  protein_of_interest <- protein_assembly_data$protein_of_interest
  protein <- protein_assembly_data$protein
  assembly <- protein_assembly_data$assembly
  protein.assembly <- protein_assembly_data$protein_assembly

  # Generate protein alias data from files in representative folder
  PROTEIN_ALIAS <- read_representatives(PATH = PATH) 
  
  
  # Check if the output file exists and if its already present read it
  if (file.exists(output_file)) {
    all.neighbours <- read_csv(output_file)
  } else {
   # Get neighboring proteins
    all.neighbours <- collec_all_neigbour(protein.assembly, BASEPAIRS, MAX_NEIGHBORS, PATH)
   # Save the results
    write_csv(all.neighbours, output_file)
  }

  # Plot the neighbors if a valid protein of interest is provided
  if (is.null(protein_of_interest) || protein_of_interest == "" || !(protein_of_interest %in% PROTEIN_ALIAS$alias)) {
    print('No valid protein of interest given, skipping plotting.')
  } else {
    plot_neighbours(all.neighbours, protein_of_interest)
  }
  
  # Geting Clade information from text files
  clades <- read_clades(PATH = PATH)
  
  # Generate COG/CDD information from all neighbours proteins based on sequence 
  cog_neighbours <- analyze_proteins_cog_classifier(df = all.neighbours, column= 'ID' )
  
  # Get the types of neighbours and their counts
  # This doesnt seem to work so far... the code that is printed in the terminal needs to be copied and run "by hand"
  # Also not very elegant, since i inserted a sleeper function to make the download more stable it now takes the program
  # 0.5 sec per Accession, which translates to 30 min per 3600 accession, which is easily reached 
  amount_of_neighbours(cog_neighbours)
  
  # Re-import the types of neighbours after manual annotation if present
  annotated_neighbours <- read_annotations(current_date)
  
  # Combine and plot
  combined_df <- combine_and_plot(neighbours_data = all.neighbours, cog_data = cog_neighbours, clade_assign = clades, neighbour_annotations = annotated_neighbours)

  # Plot data
  plot_neighbours_per_clade2(combined_df) #more pretty plot without unkown clades and unkown neighbours
  plot_neighbours_per_clade(combined_df) #plot with unkown clades and unkown neighbours
}

#main()
