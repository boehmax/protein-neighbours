# Genomic Environment Analysis

This script is designed to explore and analyse the genomic environment of a given protein. It loads necessary libraries, sources other scripts, reads in data, performs analysis, and generates plots.

## Libraries Used
- tidyverse
- ape
- dplyr
- gggenes
- ggplot2

## Scripts Sourced
- 01_open.R
- 02_clean.R
- 03_functions.R

## Main Function
The main function of the script reads protein and assembly data, generates protein alias data, collects neighboring proteins, checks if results are saved, plots the neighbors, reads clade information, reads cluster domain information, gets the types of neighbours and their counts, re-imports the types of neighbours after manual annotation, and finally combines and plots the data.

## Parameters
- BASEPAIRS: The amount of base pairs between two proteins to still consider them as a neighbour (default is 300)
- MAX_NEIGHBORS: The maximum number of neighbors to consider (default is 15)
- protein_of_interest: The protein to focus on for closer analysis (optional)

## Outputs
The script generates a CSV file with the collected neighboring proteins and a plot of the neighbors. It also reads in annotated neighbours from a CSV file and combines this with the collected data for final output.

## Usage
To run the script, simply source it in your R environment and call the main function.

## Note
Please ensure that all the necessary data files are in the correct directories as specified in the script. Also, make sure that the necessary libraries are installed in your R environment.
