#' Plot neighbors per clade
#'
#' @param combined_data A combined data frame with neighbors and clade information.
#' @param exclude_unknown_clade Logical indicating whether to exclude "unknown" clades. Default is FALSE.
#' @param exclude_unknown_cog Logical indicating whether to exclude "unknown" COG_LETTER. Default is FALSE.
#' @param output_path A string representing the output path. Default is NULL.
#' @param plot_count_codh Logical indicating whether to plot the count of CODH per clade or count of Neighbours per clade. Default is FALSE.
#' @export
plot_neighbours_per_clade <- function(combined_data, exclude_unknown_clade = FALSE, exclude_unknown_cog = FALSE, output_path = NULL, plot_count_codh = FALSE) {
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  if(is.null(output_path)) {
    output_path <- file.path("output", current_date)
  }else{
    output_path <- file.path("output", current_date, output_path)
  }
  if(exclude_unknown_clade || exclude_unknown_cog) {
    output_path <- file.path(output_path, "exclude_unknowns")
  }
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  if(plot_count_codh) {
    combined_data <- combined_data %>% mutate(ID = PIGI)
    output_path <- file.path(output_path, "count_codh")
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)}
    
  neighbour_count_per_clade <- prepare_neighbour_count(combined_data, exclude_unknown_clade, exclude_unknown_cog)
  representatives_per_clade <- prepare_representatives(combined_data, exclude_unknown_clade)
  total_neighbours_per_clade <- prepare_total_neighbours(combined_data, exclude_unknown_clade)
  how_many_neighbours_per_protein <- prepare_neighbour_count_per_protein(combined_data, exclude_unknown_clade, exclude_unknown_cog)

  write_csv(neighbour_count_per_clade, file.path(output_path, "neighbour_count_per_clade.csv"))
  write_csv(representatives_per_clade, file.path(output_path, "representatives_per_clade.csv"))
  write_csv(total_neighbours_per_clade, file.path(output_path, "total_neighbours_per_clade.csv"))
  write_csv(how_many_neighbours_per_protein, file.path(output_path, "how_many_neighbours_per_protein.csv"))

  clade_labels <- create_clade_labels(representatives_per_clade)
  plot_height <- calculate_plot_height(neighbour_count_per_clade)
  
  p <- create_plot(neighbour_count_per_clade, clade_labels, plot_height)
  print(p)
  
  # Generate color vector
  col_vector <- generate_color_vector()
  
  # Save plots with color vector
  save_plot_with_colors(p, file.path(output_path, "amount_of_neighbour_per_cladenounknown.png"), col_vector)
  save_plot_with_colors(p, file.path(output_path, "amount_of_neighbour_per_cladenounknown.svg"), col_vector)
}

prepare_neighbour_count <- function(data, exclude_unknown_clade, exclude_unknown_cog) {
  if (exclude_unknown_clade) {
    data <- data %>% filter(clade != "unknown")
  }
  
  if (exclude_unknown_cog) {
    data <- data %>% filter(COG_LETTER != "unknown")
  }
  
  # Count proteins with distinct COG_LETTER when is.neighbour is TRUE
neighbours_count <- data %>%
    filter(is.neighbour == TRUE)%>%
    count(clade, COG_LETTER, .drop = FALSE) %>%
    replace_na(list(clade = "unknown", COG_LETTER = "unknown")) 

# Count proteins when is.neighbour is FALSE
protein_count <- data %>%
    filter(is.neighbour == FALSE) %>%
    count(clade, .drop = FALSE) %>%
    replace_na(list(clade = "unknown"))

# Subtract the counts of proteins with neighbours from the counts of proteins without neighbours
no_neighbours_count <- protein_count %>%
    left_join(neighbours_count %>% group_by(clade) %>% summarise(total_neighbours = sum(n)), by = "clade") %>%
    mutate(n = n - coalesce(total_neighbours, 0)) %>%
    select(clade, n) %>%
    mutate(COG_LETTER = "no neighbours")

# Combine the two dataframes
combined_count <- bind_rows(neighbours_count, no_neighbours_count)

# Select and arrange the columns as required
final_output <- combined_count %>%
    select(clade, COG_LETTER, n) %>%
    arrange(clade, COG_LETTER)%>%
    mutate(COG_LETTER = as.factor(COG_LETTER))

return(final_output)
}

prepare_representatives <- function(data, exclude_unknown_clade) {
  if (exclude_unknown_clade) {
    data <- data %>% filter(clade != "unknown")
  }

  data <- data %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    replace_na(list(clade = "unknown")) %>%
    distinct(PIGI, clade)

  data %>%
    count(clade, .drop = FALSE) %>%
    mutate(clade = as.factor(clade))
}

prepare_total_neighbours <- function(data, exclude_unknown_clade) {
 
  if (exclude_unknown_clade) {
    data <- data %>% filter(clade != "unknown")
  }

  data <- data %>%
    filter(is.neighbour == TRUE) %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    distinct(COG_LETTER, clade)
  
  data %>%
    count(clade, .drop = FALSE) %>%
    mutate(across(where(is.factor), as.character)) %>%
    replace_na(list(clade = "unknown"))
}

prepare_neighbour_count_per_protein <- function(data, exclude_unknown_clade, exclude_unknown_cog) {
  if (exclude_unknown_clade) {
    data <- data %>% filter(clade != "unknown")
  }
  
  if (exclude_unknown_cog) {
    data <- data %>% filter(COG_LETTER != "unknown")
  }
  
  data %>%
    count(PIGI, clade, .drop = FALSE) %>%
    replace_na(list(clade = "unknown")) %>%
    mutate(n = n - 1)
}

create_clade_labels <- function(representatives) {
  clade_labels <- paste(c(LETTERS[1:6], 'hcp', 'unknown clade'), ", n=", representatives$n, sep="")
  names(clade_labels) <- c(LETTERS[1:7], 'unknown')
  clade_labels
}

calculate_plot_height <- function(neighbour_count) {
  plot_height <- round_any(max(neighbour_count$n), 100, f = ceiling)
  if (!is.finite(plot_height) || plot_height <= 0) {
    warning("Calculated plot height is not a finite positive number. Using ggplot2 default.")
    plot_height <- NULL
  }
  plot_height
}

create_plot <- function(neighbour_count, clade_labels, plot_height) {
  p <- ggplot(neighbour_count, aes(x= COG_LETTER, y= n, fill=COG_LETTER)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 1) +
    facet_grid(~clade, scales = "free_y", labeller = labeller(clade = clade_labels)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          axis.ticks.length = unit(.2, "cm"),
          axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          strip.background = element_blank(),
          panel.spacing.x = unit(-0.1, "cm"),
          legend.position = "right") +
    scale_x_discrete(expand = c(0, 1)) +
    labs(y = "Count per clade", fill = "Types of Neighbours") +
    coord_cartesian(ylim = c(0, plot_height), clip = 'off') +
    geom_text(aes(label = n), hjust = -0.2, colour = 'black', size = 2, angle = 90)
  
  if (!is.null(plot_height)) {
    p <- p + scale_y_continuous(breaks = seq(0, plot_height, plot_height / 5), 
                                limits = c(0, plot_height), 
                                expand = c(0, 0))
  }
  p
}

generate_color_vector <- function() {
  n <- 20
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector
}

save_plot_with_colors <- function(plot, filename, col_vector) {
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave(filename = filename, plot = plot, width = 30, height = 15, units = "cm")
  )
}

####################

#' Calculate Correlation Matrix
#'
#' This function calculates a correlation matrix based on the co-occurrence of elements in a data frame.
#'
#' @param df A data frame containing the data to calculate the correlation matrix.
#' @param vector A vector of column names to be used in the correlation matrix.
#' @return A correlation matrix.
#' @import dplyr
#' @export
calculate_correlation <- function(df, vector) {
  matrix <- as.data.frame.matrix(table(na.omit(df)))
  correlation.matrix <- matrix(nrow = length(vector), ncol= length(vector))
  colnames(correlation.matrix) <- vector
  rownames(correlation.matrix) <- vector
  
  for (i in seq_along(vector)) {
    for (l in seq_along(vector)) {
      var1 <- nrow(matrix %>% filter(matrix[i] >= 1 & matrix[l] >= 1)) / nrow(matrix %>% filter(matrix[i] >= 1)) 
      correlation.matrix[i,l] <- var1
    }
  }
  
  return(correlation.matrix)
}

#' Plot Correlation Matrix
#'
#' This function plots a correlation matrix.
#'
#' @param correlation.matrix A correlation matrix to be plotted from out of the calculate_correlation function.
#' @return A ggplot object representing the correlation matrix plot.
#' @import ggplot2
#' @import reshape2
#' @export
plot_correlation_matrix <- function(correlation.matrix) {
  long <- reshape2::melt(correlation.matrix)
  
  output.p <- ggplot2::ggplot(long) + 
    ggplot2::geom_tile(aes(x=Var1, y=Var2, fill=value)) +
    ggplot2::geom_text(aes(x=Var1, y=Var2, label=round(value, 2)), size=4, col = "black") +
    ggplot2::scale_fill_gradient(low = "white", high = "#66C2A5") +
    ggplot2::ylab('... how likely is to have a CODH from Clade ...') + 
    ggplot2::xlab('If an organism has a CODH from Clade...') +
    ggplot2::scale_x_discrete(position = "top") +
    theme(legend.position="none")
  
  return(output.p)
}

#' Create and Save Correlation Matrix Plot
#'
#' This function calculates a correlation matrix, plots it, and saves the plot as a PNG file.
#'
#' @param df A data frame containing the data to calculate the correlation matrix.
#' @param vector A vector of column names to be used in the correlation matrix.
#' @param supress_output Logical indicating whether to suppress the output. Default is FALSE.
#' @return A ggplot object representing the correlation matrix plot.
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @export
make_correlation_matrix <- function(df, vector, supress_output = FALSE) {
  correlation.matrix <- calculate_correlation(df, vector)
  plot <- plot_correlation_matrix(correlation.matrix)
  ggsave(paste('output/', Sys.Date(), '/correlation_matrix.png', sep=''), plot, width = 10, height = 10, units = "cm")
  if(supress_output == FALSE) {
     return(plot)}
}


#########################

#' Calculate the Number of Clades per Organism
#'
#' This function calculates the number of clades per organism from the provided FASTA data.
#'
#' @param fasta_data A data frame containing FASTA data with columns for organism and clade.
#' @return A data frame with the number of clades per organism.
#' @import dplyr
#' @export
how_many_clades_per_assembly <- function(fasta_data) {
  # Create a data frame with assembly and clade
  assembly.and.clade <- data.frame(fasta_data$assembly, fasta_data$clade)
  
  # Create a matrix from the data frame
  assembly.and.clade.matrix <- as.data.frame.matrix(table(assembly.and.clade))
  
  # Calculate the total for each clade
  assembly.and.clade.matrix$total <- rowSums(assembly.and.clade.matrix)
  
  return(assembly.and.clade.matrix)
}

#' Create Clade Histograms
#'
#' This function creates histograms for each clade and combines them into a single figure.
#'
#' @param fasta_data A data frame containing FASTA data with columns for organism and clade.
#' @param clade_colors A vector of colors to be used for the clades. Default is a predefined set of colors.
#' @return A plotly object representing the combined histograms for each clade.
#' @import ggplot2
#' @import reshape2
#' @import plotly
#' @import dplyr
#' @export
create_clade_histograms <- function(fasta_data, clade_colors = c("#FFD92F","#A6D854","#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#56B4E9","#E5C494","#B3B3B3")) {
  fasta_data.matrix <- how_many_clades_per_assembly(fasta_data)
  
  # Melt the data and filter rows where value > 0
  filtered_data <- reshape2::melt(fasta_data.matrix) %>%
    filter(value > 0)
  
  # Initialize a list to store all plots
  clade_histograms <- list()
  
  # Determine the number of clades dynamically
  num_clades <- length(unique(fasta_data$clade))
  
  # Loop over each clade
  for(i in seq_len(num_clades)){
    # Subset the data for the current clade
    clade_data <- subset(filtered_data, as.character(variable) == as.character(unique(fasta_data$clade)[i]))
    
    # Create a histogram for the current clade
    clade_histograms[[i]] <- ggplot(clade_data, aes(x=value, fill = variable)) + 
      geom_histogram(col = 'white', binwidth = 1, fill = clade_colors[i]) +
      xlim(0,8)
  }
  
  # Combine all clade histograms into a single figure
  combined_histogram <- subplot(clade_histograms, nrows = 2) %>%
    layout(title = 'Multiples of one clade distribution')
  
  # Initialize a list to store the annotations
  annotations <- list()
  
  # Loop over each clade to generate the annotations
  for(i in seq_len(num_clades)){
    annotations[[i]] <- list(
      x = (i-1)/num_clades + 1/(2*num_clades), 
      y = ifelse(i %% 2 == 0, 0.35, 0.9), 
      text = paste("Clade", as.character(unique(fasta_data$clade)[i])), 
      showarrow = F, 
      xref='paper', 
      yref='paper', 
      xanchor = "center",
      showarrow = FALSE
    )
  }
  
  # Add annotations to the combined histogram
  annotated_histogram <- combined_histogram %>% layout(annotations = annotations)
  return(annotated_histogram)
}

#' Create Clade Histograms (Alternative)
#'
#' This function creates histograms for each clade and combines them into a single figure using ggplot2.
#'
#' @param fasta_data A data frame containing FASTA data with columns for organism and clade.
#' @param clade_colors A vector of colors to be used for the clades. Default is a predefined set of colors.
#' @return A ggplot object representing the combined histograms for each clade.
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @export
create_clade_histograms2 <- function(fasta_data, clade_colors = c("#FFD92F","#A6D854","#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#56B4E9","#E5C494","#B3B3B3")) {
  fasta_data.matrix <- how_many_clades_per_assembly(fasta_data)
  
  # Melt the data and filter rows where value > 0
  filtered_data <- suppressMessages(suppressWarnings(reshape2::melt(fasta_data.matrix))) %>%
    filter(value > 0)
  
  # Determine the number of clades dynamically
  num_clades <- length(unique(fasta_data$clade))
  
  # Create a histogram for each clade
  clade_histogram <- ggplot(filtered_data, aes(x=value, fill = variable)) + 
    geom_histogram(binwidth = 1) +
    stat_count(geom = 'text',
               size = 2,
               aes(label = after_stat(count)),
               position = position_stack(vjust = 1)) +
    scale_fill_manual(values = clade_colors) +
    facet_wrap(~variable, nrow = 2) +
    theme(strip.text = element_text(size = 12), legend.position = "none") +
    labs(x ="Number of CODH in one organism", y = "Count of Organism")
  
  ggsave(paste('output/', Sys.Date(), '/clade_histogram.png', sep=''), clade_histogram, width = 10, height = 10, units = "cm")
  return(clade_histogram)
}