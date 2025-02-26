#' Main script for Protein Neighborhood Analysis
#'
#' This script coordinates the analysis of the genomic environment of proteins.
#' It loads configuration, sets up logging, and runs the analysis pipeline.
#'
#' @author Maximilian Böhm
#' @version 0.2.0

# Load required libraries
library(tidyverse)
library(ape)
library(dplyr)
library(gggenes)
library(ggplot2)
library(RColorBrewer)
library(yaml) # For configuration

# Source R scripts
source('R/utils.R')        # Utility functions including configuration
source('R/01_io.R')        # Input/output functions
source('R/02_neighbors.R') # Neighbor identification functions
source('R/03_annotation.R')# Annotation functions
source('R/04_analysis.R')  # Analysis functions
source('R/05_plotting.R')  # Visualization functions

#' Main function to run the protein neighborhood analysis
#'
#' @param config_file Path to the configuration YAML file (default is "config/config.yaml")
#' @param override_params Named list of parameters that override the config file values
#' @return A list containing the analysis results
#' @export
main <- function(config_file = "config/config.yaml", override_params = NULL) {
  # Load configuration
  config <- load_config(config_file, override_params)
  
  # Set up logging
  pn_setup_logging(config)
  
  # Log start of analysis and input file information
  pn_info("Starting protein neighborhood analysis")
  pn_info(paste("Using configuration file:", config_file))
  
  # Log input file information for reproducibility
  input_file_info <- pn_input_files(config)
  
  # Create output directory for the current date
  current_date <- config$analysis$date
  output_dir <- file.path(config$paths$output_dir, current_date)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Verify output file path for neighbors
  output_file_neighbours <- create_or_validate_output_file_path(
    basepairs = config$analysis$basepairs,
    max_neighbors = config$analysis$max_neighbors,
    date = current_date,
    config = config
  )
  
  # Output file for protein information
  output_file_proteins <- file.path(output_dir, 'all_protein_info.csv')
  
  # Read protein and assembly data
  pn_info("Reading protein and assembly data")
  protein_assembly_data <- read_protein_assembly_data(
    protein_file = config$files$proteins,
    assembly_file = config$files$assemblies,
    protein_assembly_file = config$files$protein_assembly,
    PATH = config$paths$base_dir
  )
  
  # Access the variables
  protein_of_interest <- protein_assembly_data$protein_of_interest
  protein <- protein_assembly_data$protein
  assembly <- protein_assembly_data$assembly
  protein_assembly <- protein_assembly_data$protein_assembly
  
  # Generate protein alias data from files in representative folder
  pn_info("Reading protein representatives")
  PROTEIN_ALIAS <- read_representatives(
    PATH = config$paths$base_dir,
    path = dirname(config$files$representative_files$ipg),
    ipg_file = basename(config$files$representative_files$ipg),
    pdb_file = basename(config$files$representative_files$pdb),
    cluster_file = basename(config$files$representative_files$cluster)
  )
  
  # Check if the output files exist and if they're already present, read them
  if (file.exists(output_file_neighbours)) {
    pn_info("Reading existing neighbor data from file")
    all_neighbours <- read_csv(output_file_neighbours)
    
    if(file.exists(output_file_proteins)) {
      pn_info("Reading existing protein data from file")
      all_protein <- read_csv(output_file_proteins)
    } else {
      pn_info("Collecting protein information")
      all_protein <- collect_all_protein_info(protein_assembly, PATH = config$paths$base_dir)
    }
  } else {
    # Get neighboring proteins
    pn_info("Collecting neighbor information")
    all_neighbours <- collect_all_neighbours(
      protein_assembly, 
      basepairs = config$analysis$basepairs, 
      max_neighbors = config$analysis$max_neighbors, 
      PATH = config$paths$base_dir
    )
    
    pn_info("Collecting protein information")
    all_protein <- collect_all_protein_info(protein_assembly, PATH = config$paths$base_dir)
    
    # Save the results
    pn_info("Saving neighbor and protein data")
    write_csv(all_neighbours, output_file_neighbours)
    write_csv(all_protein, output_file_proteins)
  }
  
  # Combine all data
  all_data <- rbind(all_neighbours, all_protein)
  
  # Plot the neighbors if a valid protein of interest is provided
  if (is.null(protein_of_interest) || protein_of_interest == "" || !(protein_of_interest %in% PROTEIN_ALIAS$alias)) {
    pn_warn('No valid protein of interest given, skipping plotting.')
  } else {
    pn_info(paste("Plotting neighbors for protein:", protein_of_interest))
    plot_neighbours(all_neighbours, protein_of_interest)
  }
  
  # Get clade information from text files
  pn_info("Reading clade information")
  clades <- read_clades(
    PATH = config$paths$base_dir,
    clade_dir = config$files$clade_dir,
    pattern = config$files$clade_pattern
  )
  
  # Annotate proteins using eggNOG or COG
  pn_info("Annotating proteins with", config$annotation$tool)
  annotation_results <- analyze_proteins(
    df = all_neighbours, 
    column = 'ID',
    config = config
  )
  
  # Get the types of neighbours and their counts
  pn_info("Analyzing neighbor types")
  amount_of_neighbours(annotation_results)
  
  # Re-import the types of neighbours after manual annotation if present
  pn_info("Reading any manual annotations")
  annotated_neighbours <- read_annotations(current_date)
  
  # Combine and plot
  pn_info("Combining data and generating plots")
  combined_df <- combine_and_plot(
    neighbours_data = all_data, 
    cog_data = annotation_results, 
    clade_assign = clades, 
    neighbour_annotations = annotated_neighbours
  )
  
  # Plot data with various options
  pn_info("Generating visualization plots")
  
  # Standard plots
  plot_neighbours_per_clade(
    combined_df, 
    exclude_unknown_clade = config$visualization$exclude_unknown_clade, 
    exclude_unknown_cog = config$visualization$exclude_unknown_cog
  )
  
  plot_neighbours_per_clade(combined_df)
  
  # Plots with CODH count
  plot_neighbours_per_clade(
    combined_df, 
    exclude_unknown_clade = config$visualization$exclude_unknown_clade, 
    exclude_unknown_cog = config$visualization$exclude_unknown_cog, 
    plot_count_codh = TRUE
  )
  
  plot_neighbours_per_clade(combined_df, plot_count_codh = TRUE)
  
  # Correlation plot
  pn_info("Generating correlation matrix")
  make_correlation_matrix(
    combined_df %>%
      select(PIGI, assembly, clade) %>%
      unique() %>%
      select(assembly, clade), 
    unique(combined_df$clade)
  )
  
  # Clade histograms
  pn_info("Generating clade histograms")
  create_clade_histograms2(
    combined_df %>%
      select(PIGI, assembly, clade) %>%
      unique(),
    clade_colors = config$visualization$clade_colors
  )
  
  # Neighbour plot with annotated neighbours if available
  if (!is.null(annotated_neighbours) && nrow(annotated_neighbours) > 0) {
    pn_info("Generating plots with manual annotations")
    plot_neighbours_per_clade(
      combined_df %>% mutate(COG_LETTER = ANNOTATION), 
      exclude_unknown_clade = TRUE, 
      exclude_unknown_cog = TRUE, 
      output_path = "annotated_neighbours"
    )
    
    plot_neighbours_per_clade(
      combined_df %>% mutate(COG_LETTER = ANNOTATION), 
      output_path = "annotated_neighbours"
    )
    
    plot_neighbours_per_clade(
      combined_df %>% mutate(COG_LETTER = ANNOTATION), 
      exclude_unknown_clade = TRUE, 
      exclude_unknown_cog = TRUE, 
      output_path = "annotated_neighbours", 
      plot_count_codh = TRUE
    )
    
    plot_neighbours_per_clade(
      combined_df %>% mutate(COG_LETTER = ANNOTATION), 
      output_path = "annotated_neighbours", 
      plot_count_codh = TRUE
    )
  }
  
  # Generate HTML report of the analysis
  pn_info("Generating analysis report")
  generate_analysis_report(combined_df, config, output_dir)
  
  pn_info("Analysis complete. Results available in", output_dir)
  
  # Return the combined data frame and other important objects
  return(list(
    combined_df = combined_df,
    all_neighbours = all_neighbours,
    all_protein = all_protein,
    clades = clades,
    annotation_results = annotation_results,
    config = config
  ))
}

#' Generate an HTML report of the analysis
#'
#' @param combined_df The combined data frame.
#' @param config The configuration.
#' @param output_dir The output directory.
#' @return The path to the generated report.
#' @importFrom rmarkdown render
#' @export
generate_analysis_report <- function(combined_df, config, output_dir) {
  # Check if rmarkdown is installed
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    pn_warn("Package 'rmarkdown' is needed for report generation. Skipping report.")
    return(NULL)
  }
  
  # Create a temporary Rmd file
  report_template <- system.file("rmd", "analysis_report.Rmd", package = "proteinNeighbours")
  
  if (!file.exists(report_template)) {
    # If not found in package, use local template
    report_template <- "inst/rmd/analysis_report.Rmd"
    if (!file.exists(report_template)) {
      pn_warn("Report template not found. Skipping report generation.")
      return(NULL)
    }
  }
  
  # Generate the report
  report_file <- file.path(output_dir, "analysis_report.html")
  
  tryCatch({
    rmarkdown::render(
      input = report_template,
      output_file = report_file,
      params = list(
        combined_df = combined_df,
        config = config,
        output_dir = output_dir
      ),
      quiet = TRUE
    )
    pn_info("Analysis report generated:", report_file)
  }, error = function(e) {
    pn_error("Failed to generate analysis report:", e$message)
  })
  
  return(report_file)
}

# Run the main function if not sourced
if (sys.nframe() == 0) {
  main()
}
