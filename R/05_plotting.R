#' Visualization Functions
#'
#' This file contains functions for creating visualizations of protein neighbor data,
#' including bar plots, correlation matrices, and histograms.
#'
#' @author Maximilian Böhm
#' @modified Improved version with enhanced documentation and error handling

#' Plot neighbors per clade
#'
#' This function creates a bar plot showing the distribution of neighbor types per clade.
#'
#' @param combined_data A combined data frame with neighbors and clade information.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades. Default is FALSE.
#' @param exclude_unknown_cog Logical indicating whether to exclude "unknown" COG_LETTER. Default is FALSE.
#' @param output_path A string representing the output path. Default is NULL.
#' @param plot_count_codh Logical indicating whether to plot the count of CODH per clade or count of
#'        Neighbours per clade. Default is FALSE.
#' @param width Plot width in cm (default is 30).
#' @param height Plot height in cm (default is 15).
#' @return Invisible NULL, called for side effects.
#' @export
plot_neighbours_per_clade <- function(combined_data, exclude_unknown_clade = FALSE, 
                                     exclude_unknown_cog = FALSE, output_path = NULL, 
                                     plot_count_codh = FALSE, width = 30, height = 15) {
  pn_info("Plotting neighbors per clade")
  pn_info(paste("Parameters: exclude_unknown_clade =", exclude_unknown_clade, 
                ", exclude_unknown_cog =", exclude_unknown_cog,
                ", plot_count_codh =", plot_count_codh))
  
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_path)) {
    output_path <- file.path("output", current_date)
  } else {
    output_path <- file.path("output", current_date, output_path)
  }
  
  # Add subdirectory for excluded unknowns if needed
  if (exclude_unknown_clade || exclude_unknown_cog) {
    output_path <- file.path(output_path, "exclude_unknowns")
  }
  
  # Add subdirectory for CODH count if needed
  if (plot_count_codh) {
    output_path <- file.path(output_path, "count_codh")
  }
  
  # Create output directory if it doesn't exist
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data for plotting
  tryCatch({
    # Get counts of neighbors by clade and COG letter
    neighbour_count <- prepare_neighbour_count(combined_data, exclude_unknown_clade, 
                                             exclude_unknown_cog, plot_count_codh)
    
    # Get representative counts by clade
    representatives <- prepare_representatives(combined_data, exclude_unknown_clade)
    
    # Get total neighbors by clade
    total_neighbours <- prepare_total_neighbours(combined_data, exclude_unknown_clade)
    
    # Get neighbor counts per protein
    neighbours_per_protein <- prepare_neighbour_count_per_protein(combined_data, 
                                                               exclude_unknown_clade, 
                                                               exclude_unknown_cog)
    
    # Save prepared data
    readr::write_csv(neighbour_count, file.path(output_path, "neighbour_count_per_clade.csv"))
    readr::write_csv(representatives, file.path(output_path, "representatives_per_clade.csv"))
    readr::write_csv(total_neighbours, file.path(output_path, "total_neighbours_per_clade.csv"))
    readr::write_csv(neighbours_per_protein, file.path(output_path, "how_many_neighbours_per_protein.csv"))
    
    pn_info("Saved prepared data files to:", output_path)
    
    # Create clade labels
    clade_labels <- create_clade_labels(representatives)
    
    # Calculate plot height
    plot_height <- calculate_plot_height(neighbour_count)
    
    # Create the plot
    p <- create_plot(neighbour_count, clade_labels, plot_height)
    
    # Generate color vector
    col_vector <- generate_color_vector()
    
    # Save the plot
    plot_file <- file.path(output_path, "neighbour_distribution_by_clade.png")
    save_plot_with_colors(p, plot_file, col_vector, width, height)
    
    # Save SVG version
    svg_file <- file.path(output_path, "neighbour_distribution_by_clade.svg")
    save_plot_with_colors(p, svg_file, col_vector, width, height)
    
    pn_info("Saved plots to:", output_path)
    return(invisible(NULL))
    
  }, error = function(e) {
    pn_error("Failed to create neighbors per clade plot:", e$message)
    return(invisible(NULL))
  })
}

#' Prepare neighbor count data
#'
#' This function prepares the neighbor count data for plotting.
#'
#' @param data A data frame with combined neighbor data.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades.
#' @param exclude_unknown_cog Logical indicating whether to exclude "unknown" COG_LETTER.
#' @param plot_count_codh Logical indicating whether to plot CODH counts.
#' @return A data frame with neighbor count data.
#' @keywords internal
prepare_neighbour_count <- function(data, exclude_unknown_clade, exclude_unknown_cog, plot_count_codh = FALSE) {
  # Filter data based on options
  if (exclude_unknown_clade) {
    data <- data %>% dplyr::filter(clade != "unknown")
  }
  
  if (exclude_unknown_cog) {
    data <- data %>% dplyr::filter(COG_LETTER != "unknown")
  }
  
  if (plot_count_codh) {
    # Use PIGI as ID to count CODH proteins with specific neighbors
    data <- data %>% dplyr::mutate(ID = PIGI)
  }
  
  # Count proteins with distinct COG_LETTER when is.neighbour is TRUE
  neighbours_count <- data %>%
    dplyr::filter(is.neighbour == TRUE) %>%
    dplyr::count(clade, COG_LETTER, .drop = FALSE) %>%
    tidyr::replace_na(list(clade = "unknown", COG_LETTER = "unknown"))
  
  # Count proteins when is.neighbour is FALSE
  protein_count <- data %>%
    dplyr::filter(is.neighbour == FALSE) %>%
    dplyr::count(clade, .drop = FALSE) %>%
    tidyr::replace_na(list(clade = "unknown"))
  
  # Calculate proteins with no neighbours
  no_neighbours_count <- protein_count %>%
    dplyr::left_join(neighbours_count %>% 
                  dplyr::group_by(clade) %>% 
                  dplyr::summarise(total_neighbours = sum(n)), by = "clade") %>%
    dplyr::mutate(n = n - dplyr::coalesce(total_neighbours, 0)) %>%
    dplyr::select(clade, n) %>%
    dplyr::mutate(COG_LETTER = "no neighbours")
  
  # Combine the two dataframes
  combined_count <- dplyr::bind_rows(neighbours_count, no_neighbours_count)
  
  # Select and arrange the columns as required
  final_output <- combined_count %>%
    dplyr::select(clade, COG_LETTER, n) %>%
    dplyr::arrange(clade, COG_LETTER) %>%
    dplyr::mutate(COG_LETTER = as.factor(COG_LETTER))
  
  return(final_output)
}

#' Prepare representatives data
#'
#' This function prepares the representatives data for plotting.
#'
#' @param data A data frame with combined neighbor data.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades.
#' @return A data frame with representative counts by clade.
#' @keywords internal
prepare_representatives <- function(data, exclude_unknown_clade) {
  if (exclude_unknown_clade) {
    data <- data %>% dplyr::filter(clade != "unknown")
  }
  
  data <- data %>%
    dplyr::select(ID, clade, COG_LETTER, PIGI) %>%
    tidyr::replace_na(list(clade = "unknown")) %>%
    dplyr::distinct(PIGI, clade)
  
  data %>%
    dplyr::count(clade, .drop = FALSE) %>%
    dplyr::mutate(clade = as.factor(clade))
}

#' Prepare total neighbors data
#'
#' This function prepares the total neighbors data for plotting.
#'
#' @param data A data frame with combined neighbor data.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades.
#' @return A data frame with total neighbor counts by clade.
#' @keywords internal
prepare_total_neighbours <- function(data, exclude_unknown_clade) {
  if (exclude_unknown_clade) {
    data <- data %>% dplyr::filter(clade != "unknown")
  }
  
  data <- data %>%
    dplyr::filter(is.neighbour == TRUE) %>%
    dplyr::select(ID, clade, COG_LETTER, PIGI) %>%
    dplyr::distinct(COG_LETTER, clade)
  
  data %>%
    dplyr::count(clade, .drop = FALSE) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character)) %>%
    tidyr::replace_na(list(clade = "unknown"))
}

#' Prepare neighbor count per protein data
#'
#' This function prepares the neighbor count per protein data.
#'
#' @param data A data frame with combined neighbor data.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades.
#' @param exclude_unknown_cog Logical indicating whether to exclude "unknown" COG_LETTER.
#' @return A data frame with neighbor counts per protein.
#' @keywords internal
prepare_neighbour_count_per_protein <- function(data, exclude_unknown_clade, exclude_unknown_cog) {
  if (exclude_unknown_clade) {
    data <- data %>% dplyr::filter(clade != "unknown")
  }
  
  if (exclude_unknown_cog) {
    data <- data %>% dplyr::filter(COG_LETTER != "unknown")
  }
  
  data %>%
    dplyr::count(PIGI, clade, .drop = FALSE) %>%
    tidyr::replace_na(list(clade = "unknown")) %>%
    dplyr::mutate(n = n - 1)  # Subtract 1 to account for the protein itself
}

#' Create clade labels
#'
#' This function creates labeled strings for clades.
#'
#' @param representatives A data frame with representative counts by clade.
#' @return A named vector with clade labels.
#' @keywords internal
create_clade_labels <- function(representatives) {
  # Create labels with clade and count
  clade_labels <- paste(c(LETTERS[1:6], 'hcp', 'unknown clade'), 
                      ", n = ", representatives$n, sep = "")
  names(clade_labels) <- c(LETTERS[1:7], 'unknown')
  return(clade_labels)
}

#' Calculate plot height
#'
#' This function calculates an appropriate plot height based on the maximum count.
#'
#' @param neighbour_count A data frame with neighbor counts.
#' @return A numeric value for plot height, or NULL if invalid.
#' @keywords internal
calculate_plot_height <- function(neighbour_count) {
  # Round up to nearest 100
  plot_height <- round_any(max(neighbour_count$n), 100, f = ceiling)
  
  # Check if plot_height is valid
  if (!is.finite(plot_height) || plot_height <= 0) {
    pn_info("Calculated plot height is not a finite positive number. Using ggplot2 default.")
    plot_height <- NULL
  }
  
  return(plot_height)
}

#' Create bar plot of neighbor counts
#'
#' This function creates a bar plot of neighbor counts by clade and COG letter.
#'
#' @param neighbour_count A data frame with neighbor counts.
#' @param clade_labels A named vector with clade labels.
#' @param plot_height A numeric value for plot height, or NULL for default.
#' @return A ggplot object.
#' @keywords internal
create_plot <- function(neighbour_count, clade_labels, plot_height) {
  p <- ggplot2::ggplot(neighbour_count, 
                     ggplot2::aes(x = COG_LETTER, y = n, fill = COG_LETTER)) + 
    # Add bars with identity statistic, position them side by side, and outline in black
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black", width = 1) +
    
    # Create facets based on the clade, with free y scales and labels
    ggplot2::facet_grid(~clade, scales = "free_y", 
                      labeller = ggplot2::labeller(clade = clade_labels)) +
    
    # Use a minimal theme with base font size 12, and place the legend at the right
    ggplot2::theme_bw(base_size = 12) +
    
    # Customize theme elements
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(colour = "black", size = 0.5, linetype = "solid"),
      axis.ticks.length = ggplot2::unit(.2, "cm"),
      axis.line.y = ggplot2::element_line(colour = "black", size = 0.5, linetype = "solid"),
      strip.background = ggplot2::element_blank(),
      panel.spacing.x = ggplot2::unit(-0.1, "cm"),
      legend.position = "right"
    ) +
    
    # Set x axis expansion to create space between axis and bars
    ggplot2::scale_x_discrete(expand = c(0, 1)) +
    
    # Set y axis label and fill legend title
    ggplot2::labs(y = "Count per clade", fill = "Types of Neighbours") +
    
    # Control relation to data and y axis
    ggplot2::coord_cartesian(ylim = c(0, plot_height), clip = 'off') +
    
    # Add text labels above the bars
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2, 
                     colour = 'black', size = 2, angle = 90)
  
  # Set y axis breaks and limits if plot_height is provided
  if (!is.null(plot_height)) {
    p <- p + ggplot2::scale_y_continuous(
      breaks = seq(0, plot_height, plot_height / 5),
      limits = c(0, plot_height),
      expand = c(0, 0)
    )
  }
  
  return(p)
}

#' Generate color vector
#'
#' This function generates a color vector for plots.
#'
#' @param n The number of colors to generate (default is 20).
#' @return A vector of color codes.
#' @keywords internal
generate_color_vector <- function(n = 20) {
  # Get qualitative color palettes
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  
  # Generate color vector
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, 
                             rownames(qual_col_pals)))
  
  return(col_vector)
}

#' Save plot with colors
#'
#' This function saves a plot with custom colors.
#'
#' @param plot The ggplot object to save.
#' @param filename The output filename.
#' @param col_vector A vector of color codes.
#' @param width Plot width in cm (default is 30).
#' @param height Plot height in cm (default is 15).
#' @return Invisible NULL, called for side effects.
#' @keywords internal
save_plot_with_colors <- function(plot, filename, col_vector, width = 30, height = 15) {
  # Create directory if it doesn't exist
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  # Save plot with color vector
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      width = width,
      height = height,
      units = "cm"
    )
  )
  
  pn_debug("Saved plot to:", filename)
  return(invisible(NULL))
}

#' Plot Correlation Matrix
#'
#' This function plots a correlation matrix.
#'
#' @param correlation.matrix A correlation matrix to be plotted.
#' @return A ggplot object representing the correlation matrix plot.
#' @export
plot_correlation_matrix <- function(correlation.matrix) {
  pn_info("Plotting correlation matrix")
  
  # Check if input is valid
  if (is.null(correlation.matrix) || !is.matrix(correlation.matrix)) {
    pn_error("Invalid correlation matrix provided")
    return(NULL)
  }
  
  # Convert matrix to long format for ggplot
  tryCatch({
    long <- reshape2::melt(correlation.matrix)
    
    # Create the plot
    output.p <- ggplot2::ggplot(long) + 
      ggplot2::geom_tile(ggplot2::aes(x = Var1, y = Var2, fill = value)) +
      ggplot2::geom_text(ggplot2::aes(x = Var1, y = Var2, label = round(value, 2)), 
                       size = 4, col = "black") +
      ggplot2::scale_fill_gradient(low = "white", high = "#66C2A5") +
      ggplot2::ylab('... how likely is to have a CODH from Clade ...') + 
      ggplot2::xlab('If an organism has a CODH from Clade...') +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
      ) +
      ggplot2::ggtitle("Co-occurrence Probability of Clades")
    
    pn_info("Successfully created correlation matrix plot")
    return(output.p)
    
  }, error = function(e) {
    pn_error("Failed to create correlation matrix plot:", e$message)
    return(NULL)
  })
}

#' Create and Save Correlation Matrix Plot
#'
#' This function calculates a correlation matrix, plots it, and saves the plot as a PNG file.
#'
#' @param df A data frame containing the data to calculate the correlation matrix.
#' @param vector A vector of column names to be used in the correlation matrix.
#' @param supress_output Logical indicating whether to suppress the output. Default is FALSE.
#' @param output_dir Directory to save the plot (default is based on current date).
#' @param width Plot width in cm (default is 10).
#' @param height Plot height in cm (default is 10).
#' @return A ggplot object representing the correlation matrix plot if supress_output is FALSE.
#' @export
make_correlation_matrix <- function(df, vector, supress_output = FALSE, 
                                   output_dir = NULL, width = 10, height = 10) {
  pn_info("Creating and saving correlation matrix")
  
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_dir)) {
    output_dir <- file.path('output', current_date)
  }
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate correlation matrix
  correlation.matrix <- calculate_correlation(df, vector)
  
  if (is.null(correlation.matrix)) {
    pn_info("Failed to calculate correlation matrix")
    return(NULL)
  }
  
  # Create plot
  plot <- plot_correlation_matrix(correlation.matrix)
  
  if (is.null(plot)) {
    pn_info("Failed to create correlation matrix plot")
    return(NULL)
  }
  
  # Save plot
  output_file <- file.path(output_dir, 'correlation_matrix.png')
  ggplot2::ggsave(output_file, plot, width = width, height = height, units = "cm")
  pn_info("Saved correlation matrix plot to:", output_file)
  
  # Save matrix as CSV
  output_file_csv <- file.path(output_dir, 'correlation_matrix.csv')
  readr::write_csv(as.data.frame(correlation.matrix), output_file_csv)
  pn_info("Saved correlation matrix data to:", output_file_csv)
  
  # Return plot if not suppressed
  if (supress_output == FALSE) {
    return(plot)
  }
  
  return(invisible(NULL))
}

#' Create Clade Histograms
#'
#' This function creates histograms for each clade showing the distribution of
#' proteins per organism.
#'
#' @param fasta_data A data frame containing data with columns for organism and clade.
#' @param clade_colors A vector of colors to be used for the clades. Default is a predefined set of colors.
#' @param output_dir Directory to save the plot (default is based on current date).
#' @param width Plot width in cm (default is 10).
#' @param height Plot height in cm (default is 10).
#' @return A ggplot object representing the combined histograms.
#' @export
create_clade_histograms2 <- function(fasta_data, 
                                   clade_colors = c("#FFD92F", "#A6D854", "#FC8D62", 
                                                   "#E78AC3", "#8DA0CB", "#66C2A5", 
                                                   "#56B4E9", "#E5C494", "#B3B3B3"),
                                   output_dir = NULL, width = 10, height = 10) {
  pn_info("Creating clade histograms")
  
  # Set output directory
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if (is.null(output_dir)) {
    output_dir <- file.path('output', current_date)
  }
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate clade distribution per assembly
  tryCatch({
    fasta_data.matrix <- how_many_clades_per_assembly(fasta_data)
    
    if (is.null(fasta_data.matrix)) {
      pn_info("Failed to calculate clades per assembly")
      return(NULL)
    }
    
    # Melt the data and filter rows where value > 0
    filtered_data <- suppressMessages(suppressWarnings(reshape2::melt(fasta_data.matrix))) %>%
      dplyr::filter(value > 0)
    
    # Determine the number of clades dynamically
    num_clades <- length(unique(fasta_data$clade))
    
    # Create a histogram for each clade
    clade_histogram <- ggplot2::ggplot(filtered_data, 
                                     ggplot2::aes(x = value, fill = variable)) + 
      ggplot2::geom_histogram(binwidth = 1) +
      ggplot2::stat_count(geom = 'text',
                        size = 2,
                        ggplot2::aes(label = ggplot2::after_stat(count)),
                        position = ggplot2::position_stack(vjust = 1)) +
      ggplot2::scale_fill_manual(values = clade_colors) +
      ggplot2::facet_wrap(~variable, nrow = 2) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 12),
        legend.position = "none",
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
      ) +
      ggplot2::labs(
        x = "Number of CODH in one organism",
        y = "Count of Organisms",
        title = "Distribution of Proteins per Organism by Clade"
      )
    
    # Save the plot
    output_file <- file.path(output_dir, 'clade_histogram.png')
    ggplot2::ggsave(output_file, clade_histogram, width = width, height = height, units = "cm")
    pn_info("Saved clade histogram to:", output_file)
    
    # Save the data
    output_file_csv <- file.path(output_dir, 'clade_histogram_data.csv')
    readr::write_csv(filtered_data, output_file_csv)
    pn_info("Saved clade histogram data to:", output_file_csv)
    
    return(clade_histogram)
    
  }, error = function(e) {
    pn_error("Failed to create clade histograms:", e$message)
    return(NULL)
  })
}
