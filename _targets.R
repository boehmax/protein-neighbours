#' Targets Pipeline for Protein Neighborhood Analysis
#'
#' This file defines the targets pipeline using the {targets} package.
#' The pipeline orchestrates the analysis of the genomic environment of proteins
#' with automatic dependency tracking, caching, and parallel execution support.
#'
#' To run the pipeline:
#'   targets::tar_make()
#'
#' To visualize the pipeline:
#'   targets::tar_visnetwork()
#'
#' To load results:
#'   targets::tar_read(target_name)
#'
#' @author Maximilian Böhm
#' @version 0.2.0

# Load required packages
library(targets)
library(tarchetypes)

# Source all R functions
source("R/utils.R")
source("R/01_io.R")
source("R/02_neighbors.R")
source("R/03_annotation.R")
source("R/04_analysis.R")
source("R/05_plotting.R")

# Source targets-specific functions
source("R/targets_functions.R")

# Set target-specific options
tar_option_set(
  packages = c(
    "tidyverse", "dplyr", "ggplot2", "svglite", "ape",
    "gggenes", "RColorBrewer", "yaml", "rmarkdown", "logger"
  ),
  format = "rds",  # Default format for storing targets
  error = "continue"  # Continue with other targets if one fails
)

# Define the pipeline
list(
  # Target 1: Load configuration
  tar_target(
    config,
    load_config("config/config.yaml"),
    deployment = "main"
  ),
  
  # Target 2: Setup logging
  tar_target(
    logging_setup,
    {
      pn_setup_logging(config)
      pn_info("Starting protein neighborhood analysis with targets pipeline")
      pn_info(paste("Using configuration file: config/config.yaml"))
      TRUE
    },
    deployment = "main"
  ),
  
  # Target 3: Create output directory
  tar_target(
    output_dir,
    {
      current_date <- config$analysis$date
      out_dir <- file.path(config$paths$output_dir, current_date)
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }
      pn_info(paste("Output directory:", out_dir))
      out_dir
    },
    deployment = "main"
  ),
  
  # Target 4: Define output file paths
  tar_target(
    output_files,
    {
      output_file_neighbours <- create_or_validate_output_file_path(
        basepairs = config$analysis$basepairs,
        max_neighbors = config$analysis$max_neighbors,
        date = config$analysis$date,
        config = config
      )
      output_file_proteins <- file.path(output_dir, 'all_protein_info.csv')
      list(
        neighbours = output_file_neighbours,
        proteins = output_file_proteins
      )
    },
    deployment = "main"
  ),
  
  # Target 5: Read protein and assembly data
  tar_target(
    protein_assembly_data,
    {
      logging_setup  # Ensure logging is set up
      pn_info("Reading protein and assembly data")
      read_protein_assembly_data(
        protein_file = config$files$proteins,
        assembly_file = config$files$assemblies,
        protein_assembly_file = config$files$protein_assembly,
        PATH = config$paths$base_dir
      )
    }
  ),
  
  # Target 6: Read protein representatives
  tar_target(
    protein_alias,
    {
      pn_info("Reading protein representatives")
      read_representatives(
        PATH = config$paths$base_dir,
        path = dirname(config$files$representative_files$ipg),
        ipg_file = basename(config$files$representative_files$ipg),
        pdb_file = basename(config$files$representative_files$pdb),
        cluster_file = basename(config$files$representative_files$cluster)
      )
    }
  ),
  
  # Target 7: Collect neighbor information (with file caching)
  tar_target(
    all_neighbours,
    {
      if (file.exists(output_files$neighbours)) {
        pn_info("Reading existing neighbor data from file")
        readr::read_csv(output_files$neighbours, show_col_types = FALSE)
      } else {
        pn_info("Collecting neighbor information")
        neighbours <- collect_all_neighbours(
          protein_assembly_data$protein_assembly,
          basepairs = config$analysis$basepairs,
          max_neighbors = config$analysis$max_neighbors,
          PATH = config$paths$base_dir,
          overlap = config$analysis$overlap
        )
        pn_info("Saving neighbor data")
        readr::write_csv(neighbours, output_files$neighbours)
        neighbours
      }
    }
  ),
  
  # Target 8: Collect protein information (with file caching)
  tar_target(
    all_protein,
    {
      if (file.exists(output_files$proteins)) {
        pn_info("Reading existing protein data from file")
        readr::read_csv(output_files$proteins, show_col_types = FALSE)
      } else {
        pn_info("Collecting protein information")
        protein_info <- collect_all_protein_info(
          protein_assembly_data$protein_assembly,
          PATH = config$paths$base_dir
        )
        pn_info("Saving protein data")
        readr::write_csv(protein_info, output_files$proteins)
        protein_info
      }
    }
  ),
  
  # Target 9: Combine all data
  tar_target(
    all_data,
    {
      rbind(all_neighbours, all_protein)
    }
  ),
  
  # Target 10: Plot neighbors for protein of interest (optional)
  tar_target(
    neighbor_plot,
    {
      protein_of_interest <- protein_assembly_data$protein_of_interest
      if (is.null(protein_of_interest) || protein_of_interest == "" || 
          !(protein_of_interest %in% protein_alias$alias)) {
        pn_warn('No valid protein of interest given, skipping plotting.')
        NULL
      } else {
        pn_info(paste("Plotting neighbors for protein:", protein_of_interest))
        plot_neighbours(all_neighbours, protein_of_interest)
        protein_of_interest
      }
    }
  ),
  
  # Target 11: Read clade information
  tar_target(
    clades,
    {
      pn_info("Reading clade information")
      read_clades(
        PATH = config$paths$base_dir,
        clade_dir = config$files$clade_dir,
        pattern = config$files$clade_pattern
      )
    }
  ),
  
  # Target 12: Annotate proteins
  tar_target(
    annotation_results,
    {
      pn_info("Annotating proteins with", config$annotation$tool)
      analyze_proteins(
        df = all_neighbours,
        column = 'ID',
        config = config
      )
    }
  ),
  
  # Target 13: Analyze neighbor types
  tar_target(
    neighbor_types,
    {
      pn_info("Analyzing neighbor types")
      amount_of_neighbours(annotation_results)
    }
  ),
  
  # Target 14: Read manual annotations (if available)
  tar_target(
    annotated_neighbours,
    {
      pn_info("Reading any manual annotations")
      read_annotations(config$analysis$date)
    }
  ),
  
  # Target 15: Combine and plot data
  tar_target(
    combined_df,
    {
      pn_info("Combining data and generating plots")
      combine_and_plot(
        neighbours_data = all_data,
        cog_data = annotation_results,
        clade_assign = clades,
        neighbour_annotations = annotated_neighbours
      )
    }
  ),
  
  # Target 16: Generate standard plots
  tar_target(
    standard_plots,
    {
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
      
      TRUE
    }
  ),
  
  # Target 17: Generate correlation matrix
  tar_target(
    correlation_matrix,
    {
      pn_info("Generating correlation matrix")
      make_correlation_matrix(
        combined_df %>%
          select(PIGI, assembly, clade) %>%
          unique() %>%
          select(assembly, clade),
        unique(combined_df$clade)
      )
    }
  ),
  
  # Target 18: Create clade histograms
  tar_target(
    clade_histograms,
    {
      pn_info("Generating clade histograms")
      create_clade_histograms2(
        combined_df %>%
          select(PIGI, assembly, clade) %>%
          unique(),
        clade_colors = config$visualization$clade_colors
      )
    }
  ),
  
  # Target 19: Generate annotated plots (if manual annotations exist)
  tar_target(
    annotated_plots,
    {
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
        
        TRUE
      } else {
        pn_info("No manual annotations available, skipping annotated plots")
        NULL
      }
    }
  ),
  
  # Target 20: Generate HTML report
  tar_target(
    analysis_report,
    {
      pn_info("Generating analysis report")
      generate_analysis_report(combined_df, config, output_dir)
    }
  ),
  
  # Target 21: Final summary
  tar_target(
    analysis_summary,
    {
      pn_info("Analysis complete. Results available in", output_dir)
      list(
        combined_df = combined_df,
        all_neighbours = all_neighbours,
        all_protein = all_protein,
        clades = clades,
        annotation_results = annotation_results,
        config = config,
        output_dir = output_dir,
        report_file = analysis_report
      )
    }
  )
)
