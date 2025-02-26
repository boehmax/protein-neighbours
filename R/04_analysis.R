#' Analysis Functions
#'
#' This file contains functions for analyzing protein neighbors,
#' combining data, and calculating statistics.
#'
#' @author Maximilian Böhm
#' @modified Improved version with enhanced documentation and error handling

#' Get the Types of Neighbours and Their Counts
#'
#' This function calculates the types of neighbours and their counts, and generates a plot and CSV file.
#' If no cluster domain assignments are found, the function uses protein product information for NCBI annotation.
#'
#' @param cog_data Data frame with COG information.
#' @param output_dir Directory to save output files (default is based on current date).
#' @return A data frame with types of neighbors and their counts.
#' @export
amount_of_neighbours <- function(cog_data, output_dir = NULL) {
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_dir)) {
    output_dir <- file.path('output', current_date)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pn_info("Analyzing neighbor types and counts")
  
  # Check if cog_data is valid
  if (is.null(cog_data) || nrow(cog_data) == 0) {
    pn_warn("No COG data provided. Cannot analyze neighbor types.")
    return(NULL)
  }
  
  # Check if COG_NAME and COG_LETTER columns exist
  if (!all(c("COG_NAME", "COG_LETTER") %in% colnames(cog_data))) {
    pn_warn("Required columns (COG_NAME, COG_LETTER) not found in COG data.")
    pn_warn("Available columns:", paste(colnames(cog_data), collapse = ", "))
    return(NULL)
  }
  
  # Calculate types of neighbors and their counts
  tryCatch({
    types_of_neighbours <- cog_data %>%
      dplyr::select(COG_NAME, COG_LETTER) %>%
      dplyr::add_count(COG_NAME) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      unique()
    
    pn_info("Identified", nrow(types_of_neighbours), "unique neighbor types")
    
    # Plot the types of neighbours in a range (top 100)
    plot_data <- types_of_neighbours
    if (nrow(plot_data) > 100) {
      plot_data <- plot_data[1:100, ]
      pn_info("Plotting top 100 neighbor types")
    }
    
    tryCatch({
      p <- ggplot2::ggplot(plot_data, 
                        ggplot2::aes(x = reorder(COG_LETTER, -n), y = n)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = "Types of Neighbors by COG Category",
                  x = "COG Category",
                  y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      # Save the plot
      output_file <- file.path(output_dir, 'types_of_neighbours.png')
      ggplot2::ggsave(output_file, plot = p, width = 10, height = 8)
      pn_info("Saved neighbor types plot to:", output_file)
    }, error = function(e) {
      pn_error("Failed to create or save neighbor types plot:", e$message)
    })
    
    # Save the types of neighbors as CSV
    output_file <- file.path(output_dir, 'types_of_neighbours.csv')
    readr::write_csv(types_of_neighbours, output_file)
    pn_info("Saved neighbor types data to:", output_file)
    
    # Print instructions for manual annotation
    message('Please open the output folder to see the types of neighbours and annotate them as per the instructions in the README.md file (optional)')
    
    return(types_of_neighbours)
    
  }, error = function(e) {
    pn_error("Failed to analyze neighbor types:", e$message)
    return(NULL)
  })
}

#' Combine and plot neighbors
#'
#' This function combines data from different sources and prepares it for visualization.
#'
#' @param neighbours_data A data frame with neighbors data.
#' @param cog_data A data frame with COG information from analyze_proteins function.
#' @param clade_assign A data frame with clade assignments.
#' @param neighbour_annotations A data frame with manual neighbor annotations.
#' @param output_dir Directory to save output files (default is based on current date).
#' @return A combined data frame.
#' @export
combine_and_plot <- function(neighbours_data, cog_data, clade_assign, neighbour_annotations, output_dir = NULL) {
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_dir)) {
    output_dir <- file.path('output', current_date)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pn_info("Combining data from multiple sources")
  
  # Check if data frames exist and have necessary columns
  if (is.null(neighbours_data) || nrow(neighbours_data) == 0) {
    pn_error("No neighbor data provided")
    return(data.frame())
  }
  
  # Ensure ID column exists for joining
  if (!"ID" %in% colnames(neighbours_data)) {
    pn_error("ID column missing from neighbors_data")
    pn_info("Available columns:", paste(colnames(neighbours_data), collapse=", "))
    return(neighbours_data)
  }
  
  # Merge the data with the clade assignment
  if (is.null(clade_assign) || nrow(clade_assign) == 0) {
    pn_warn("No clade assignment found. Using 'unknown' for all clades.")
    neighbours_with_clades <- neighbours_data %>%
      dplyr::mutate(clade = 'unknown')
  } else {
    # Check if PIGI column exists in both dataframes
    if (!"PIGI" %in% colnames(neighbours_data)) {
      pn_warn("PIGI column missing from neighbors_data. Cannot join with clade data.")
      neighbours_with_clades <- neighbours_data %>%
        dplyr::mutate(clade = 'unknown')
    } else if (!"PIGI" %in% colnames(clade_assign)) {
      pn_warn("PIGI column missing from clade_assign. Cannot join with neighbor data.")
      neighbours_with_clades <- neighbours_data %>%
        dplyr::mutate(clade = 'unknown')
    } else {
      # Merge the clade assignment
      pn_info("Merging clade assignments")
      neighbours_with_clades <- neighbours_data %>%
        dplyr::left_join(clade_assign, by = 'PIGI') %>%
        tidyr::replace_na(list(clade = 'unknown'))
      
      pn_info("Merged", nrow(clade_assign), "clade assignments with", 
              nrow(neighbours_data), "neighbor records")
    }
  }
  
  # Combine all information into one dataframe if possible
  if (is.null(cog_data) || nrow(cog_data) == 0) {
    pn_warn("No additional COG data found. Using default values.")
    combined_data <- neighbours_with_clades %>%
      dplyr::mutate(COG_ID = NA,
                    CDD_ID = NA,
                    EVALUE = NA,
                    GENE_NAME = NA,
                    COG_NAME = 'unknown',
                    COG_LETTER = 'unknown',
                    COG_DESCRIPTION = 'unknown'
      )
  } else {
    # Check column names in cog_data
    pn_info("COG data columns:", paste(colnames(cog_data), collapse=", "))
    
    # Determine the correct join column in cog_data
    join_col <- NULL
    if ("QUERY_ID" %in% colnames(cog_data)) {
      join_col <- "QUERY_ID"
    } else if ("query" %in% colnames(cog_data)) {
      join_col <- "query"
    } else {
      pn_warn("Neither QUERY_ID nor query column found in cog_data. Cannot join.")
      combined_data <- neighbours_with_clades %>%
        dplyr::mutate(COG_ID = NA,
                      CDD_ID = NA,
                      EVALUE = NA,
                      GENE_NAME = NA,
                      COG_NAME = 'unknown',
                      COG_LETTER = 'unknown',
                      COG_DESCRIPTION = 'unknown'
        )
    }
    
    if (!is.null(join_col)) {
      # Merge the COG data
      pn_info("Merging COG annotation data")
      pn_info("Joining on", "ID from neighbors_data and", join_col, "from cog_data")
      
      # Use distinct to avoid duplicate rows (many-to-many relationship warning)
      distinct_cog_data <- cog_data %>% dplyr::distinct(!!as.name(join_col), .keep_all = TRUE)
      
      combined_data <- neighbours_with_clades %>%
        dplyr::left_join(distinct_cog_data, by = c("ID" = join_col)) %>%
        tidyr::replace_na(list(COG_ID = NA,
                               CDD_ID = NA,
                               EVALUE = NA,
                               GENE_NAME = NA,
                               COG_NAME = 'unknown',
                               COG_LETTER = 'unknown',
                               COG_DESCRIPTION = 'unknown'
        ))
      
      pn_info("Merged", nrow(distinct_cog_data), "COG annotations")
    }
  } 
  
  # Add manual annotations if available
  if (!is.null(neighbour_annotations) && nrow(neighbour_annotations) > 0) {
    pn_info("Adding manual neighbor annotations")
    
    # Check if COG_NAME exists in both dataframes
    if (!"COG_NAME" %in% colnames(combined_data)) {
      pn_warn("COG_NAME column missing from combined_data. Cannot add manual annotations.")
    } else if (!"COG_NAME" %in% colnames(neighbour_annotations)) {
      pn_warn("COG_NAME column missing from neighbour_annotations. Cannot add manual annotations.")
    } else {
      # Add manual annotations
      combined_data <- combined_data %>%
        dplyr::left_join(neighbour_annotations %>% 
                           dplyr::select(COG_NAME, ANNOTATION), by = "COG_NAME") %>%
        tidyr::replace_na(list(ANNOTATION = 'unknown'))
      
      pn_info("Added", nrow(neighbour_annotations), "manual annotations")
    }
  }
  
  # Write the combined data to a CSV file
  output_file <- file.path(output_dir, "combined_df_all_neighbours_assigned.csv")
  readr::write_csv(combined_data, output_file)
  pn_info("Saved combined data to:", output_file)
  
  return(combined_data)
}

#' Calculate the Number of Clades per Organism
#'
#' This function calculates the number of clades per organism from the provided data.
#'
#' @param fasta_data A data frame containing data with columns for organism and clade.
#' @return A data frame with the number of clades per organism.
#' @export
how_many_clades_per_assembly <- function(fasta_data) {
  pn_info("Calculating number of clades per assembly")
  
  # Check if required columns exist
  if (!all(c("assembly", "clade") %in% colnames(fasta_data))) {
    pn_error("Required columns (assembly, clade) not found in input data")
    return(NULL)
  }
  
  # Create a data frame with assembly and clade
  assembly.and.clade <- data.frame(
    assembly = fasta_data$assembly, 
    clade = fasta_data$clade,
    stringsAsFactors = FALSE
  )
  
  # Create a matrix from the data frame
  tryCatch({
    assembly.and.clade.matrix <- as.data.frame.matrix(table(assembly.and.clade))
    
    # Calculate the total for each clade
    assembly.and.clade.matrix$total <- rowSums(assembly.and.clade.matrix)
    
    pn_info("Calculated clade distribution for", nrow(assembly.and.clade.matrix), "assemblies")
    
    return(assembly.and.clade.matrix)
    
  }, error = function(e) {
    pn_error("Failed to calculate clades per assembly:", e$message)
    return(NULL)
  })
}

#' Calculate Correlation Matrix
#'
#' This function calculates a correlation matrix based on the co-occurrence of elements in a data frame.
#'
#' @param df A data frame containing the data to calculate the correlation matrix.
#' @param vector A vector of column names to be used in the correlation matrix.
#' @return A correlation matrix.
#' @export
calculate_correlation <- function(df, vector) {
  pn_info("Calculating correlation matrix")
  
  tryCatch({
    # Create a contingency table
    matrix <- as.data.frame.matrix(table(na.omit(df)))
    
    # Initialize the correlation matrix
    correlation.matrix <- matrix(nrow = length(vector), ncol = length(vector))
    colnames(correlation.matrix) <- vector
    rownames(correlation.matrix) <- vector
    
    # Calculate correlation values
    for (i in seq_along(vector)) {
      for (l in seq_along(vector)) {
        # Calculate probability of having clade l given clade i
        var1 <- nrow(matrix %>% dplyr::filter(matrix[, i] >= 1 & matrix[, l] >= 1)) / 
              nrow(matrix %>% dplyr::filter(matrix[, i] >= 1))
        
        # Handle edge cases (0/0)
        if (is.nan(var1)) {
          var1 <- 0
        }
        
        correlation.matrix[i, l] <- var1
      }
    }
    
    pn_info("Successfully calculated correlation matrix")
    return(correlation.matrix)
    
  }, error = function(e) {
    pn_error("Failed to calculate correlation matrix:", e$message)
    return(NULL)
  })
}

#' Generate Summary Statistics
#'
#' This function calculates summary statistics from the combined data.
#'
#' @param combined_data The combined data frame.
#' @param output_dir Directory to save output files (default is based on current date).
#' @return A list of summary statistics.
#' @export
generate_summary_statistics <- function(combined_data, output_dir = NULL) {
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_dir)) {
    output_dir <- file.path('output', current_date)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pn_info("Generating summary statistics")
  
  # Initialize results list
  results <- list()
  
  tryCatch({
    # Count proteins and neighbors
    results$total_proteins <- nrow(combined_data %>% 
                                  dplyr::filter(!is.neighbour) %>% 
                                  dplyr::distinct(PIGI))
    
    results$total_neighbors <- nrow(combined_data %>% 
                                   dplyr::filter(is.neighbour))
    
    results$unique_neighbor_types <- length(unique(combined_data$COG_LETTER))
    
    results$clades <- unique(combined_data$clade)
    results$clade_counts <- combined_data %>%
      dplyr::filter(!is.neighbour) %>%
      dplyr::count(clade, sort = TRUE)
    
    results$assemblies <- unique(combined_data$assembly)
    results$assembly_count <- length(results$assemblies)
    
    # Proteins per assembly
    results$proteins_per_assembly <- combined_data %>%
      dplyr::filter(!is.neighbour) %>%
      dplyr::count(assembly, sort = TRUE)
    
    # Neighbors per protein
    results$neighbors_per_protein <- combined_data %>%
      dplyr::filter(is.neighbour) %>%
      dplyr::count(PIGI, sort = TRUE)
    
    # Neighbor type distribution
    results$neighbor_type_dist <- combined_data %>%
      dplyr::filter(is.neighbour) %>%
      dplyr::count(COG_LETTER, sort = TRUE)
    
    # Create summary data frame
    summary_df <- data.frame(
      Metric = c(
        "Total Proteins", 
        "Total Neighbors", 
        "Unique Neighbor Types", 
        "Number of Clades",
        "Number of Assemblies"
      ),
      Value = c(
        results$total_proteins, 
        results$total_neighbors, 
        results$unique_neighbor_types, 
        length(results$clades),
        results$assembly_count
      ),
      stringsAsFactors = FALSE
    )
    
    # Save summary statistics
    output_file <- file.path(output_dir, "summary_statistics.csv")
    readr::write_csv(summary_df, output_file)
    pn_info("Saved summary statistics to:", output_file)
    
    # Save detailed statistics
    output_file <- file.path(output_dir, "clade_distribution.csv")
    readr::write_csv(results$clade_counts, output_file)
    
    output_file <- file.path(output_dir, "neighbors_per_protein.csv")
    readr::write_csv(results$neighbors_per_protein, output_file)
    
    output_file <- file.path(output_dir, "neighbor_type_distribution.csv")
    readr::write_csv(results$neighbor_type_dist, output_file)
    
    pn_info("Successfully generated summary statistics")
    return(results)
    
  }, error = function(e) {
    pn_error("Failed to generate summary statistics:", e$message)
    return(NULL)
  })
}
