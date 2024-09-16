#' Rounding function
#'
#' @param x The number to be rounded.
#' @param accuracy The accuracy to which the number should be rounded.
#' @param f The rounding function to use (default is round).
#' @return The rounded number.
#' @export
round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

#' Extract a specific field from the attribute column of a GFF file
#'
#' @param x The attribute column as a character vector.
#' @param field The field to extract.
#' @param attrsep The separator used in the attribute column (default is ";").
#' @return A character vector with the extracted field values.
#' @export
getAttributeField <- function (x, field, attrsep = ";") {
  # Split the attributes
  s = strsplit(x, split = attrsep, fixed = TRUE)
  
  # Extract the desired field
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

#' Find the protein alias
#'
#' @param protein.id The protein ID to search for.
#' @param alias The alias database (default is PROTEIN_ALIAS).
#' @param verbose Whether to print additional information (default is FALSE).
#' @param identi Whether to include identity information (default is FALSE).
#' @return A data frame with the corresponding PIGI and alias.
#' @export
protein.alias <- function(protein.id, alias = PROTEIN_ALIAS, verbose = FALSE, identi = FALSE){
  # Check if the protein ID is in the alias database
  if(length(alias$alias == protein.id) > 0){
    # If it is, find the corresponding PIGI
    result <- alias %>%
      distinct(PIGI, alias) %>%
      filter(alias == protein.id) %>%
      select(PIGI)
    # Also find the corresponding identity
    if(identi){
      percent <- alias %>%
        filter(alias == protein.id) %>%
        select(identity)
    }
    # Print a message indicating the protein was found, if verbose is TRUE
    if(verbose) {
      if(identi) {
        print(paste("Protein (", protein.id, ") was found in alias database, with an identity of ", percent, " matching ", result, ".", sep = ""))
      } else {
        print(paste("Protein (", protein.id, ") was found in alias database matching ", result, ".", sep = ""))
      }
    }
    # Return the PIGI as a character string
    return(result)
  } else {
    # If the protein ID is not in the database, print a message indicating it was not found, if verbose is TRUE
    if(verbose) {
      print("Protein was not found in alias database.")
    }
    # Return the original protein ID
    return(protein.id)
  }
}

#' Get neighboring proteins
#'
#' @param gff.df The GFF data frame.
#' @param protein.id The protein ID to search for neighbors.
#' @param bp Number of base pairs to consider for neighbors (default is 300).
#' @param n Number of neighbors to find (default is 15).
#' @return A data frame with neighboring proteins.
#' @export
getNeighborProteins <- function(gff.df, protein.id, bp = 300, n = 15){
  # Initialize empty data frame for neighbors
  neighbors <- data.frame()
  
  # Get start, end, strand, and seqid for the given protein id
  protein_info <- gff.df %>% filter(ID == protein.id)
  start <- protein_info$start
  end <- protein_info$end
  strand <- protein_info$strand
  seqid <- protein_info$seqid
  
  # Loop through n upstream and downstream neighbors
  for(i in seq_len(n)){
    # Define condition for upstream neighbors
    condition <- start-bp < gff.df$end & gff.df$end < start & gff.df$type=='CDS' & gff.df$strand== strand & gff.df$seqid== seqid
    
    # If there are upstream neighbors, add the first one to the neighbors data frame
    if(nrow(gff.df[condition, ]) > 0){
      output <- gff.df[condition, ][1,]
      start <- output$start
      neighbors = rbind(neighbors, output)
    }
    
    # Define condition for downstream neighbors
    condition <- end+bp > gff.df$start & gff.df$start > end & gff.df$type=='CDS' & gff.df$strand== strand & gff.df$seqid== seqid
    
    # If there are downstream neighbors, add the first one to the neighbors data frame
    if(nrow(gff.df[condition, ]) > 0){
      output <- gff.df[condition, ][1,]
      end <- output$end
      neighbors = rbind(neighbors, output)
    }
  }
  
  return(neighbors)
}

#' Get protein neighbors from a GFF3 file
#'
#' @param input The path to the GFF3 file.
#' @param protein.id The protein ID to search for neighbors.
#' @param basepairs Number of base pairs to consider for neighbors (default is 300).
#' @param m Number of neighbors to find (default is 15).
#' @return A data frame with neighboring proteins.
#' @export
getProteinNeighborsFromGff3 <- function(input, protein.id, basepairs = 300, m = 15){
  # Inform the user which input is being processed
  message(paste("Processing input:", input))

   # Check if the file exists
  if (!file.exists(input)) {
    warning(paste("File not found:", input, "- Skipping this input."))
    return(NULL)
  }

  # Import gff3 file
  gffData <- read.gff(input, na.strings = c(".", "?"), GFF3 = TRUE)
  
  # Get Protein ID and product separately
  gffData <- gffData %>%
    mutate(ID = gsub(".*\\-", "", getAttributeField(attributes, "ID")),
           product = getAttributeField(attributes, "product"))
  
  # Get neighboring proteins
  finalOutput <- getNeighborProteins(gffData, protein.id, basepairs, m)
  
  return(finalOutput)
}

#' Collect all neighbors for a list of proteins
#'
#' @param protein.assembly A data frame with protein and assembly information.
#' @param basepairs Number of base pairs to consider for neighbors (default is 300).
#' @param m Number of neighbors to find (default is 15).
#' @param PATH Path where to ncbi_dataset folder is stored (default = data).
#' @return A data frame with all neighbors.
#' @export
collec_all_neigbour <- function(protein.assembly, basepairs = 300, max_neighbors = 15, PATH = 'data'){
  all.neighbours <- data.frame() # Initialize empty data frame for all neighbors
  for(i in seq_len(nrow(protein.assembly))){
    # Get neighbors for each protein
    np <- getProteinNeighborsFromGff3(paste(PATH,'/ncbi_dataset/data/', protein.assembly[i,1],'/genomic.gff', sep = ""), 
                                      protein.assembly[i,2], basepairs, max_neighbors)
    # Skip the iteration if np is NULL
    if (is.null(np)) {
        next
    }

    if(nrow(np) > 0){
      # Add protein and assembly information to neighbors data frame
      np <- np %>%
        mutate(PIGI = protein.assembly[i,2],  # Protein I Gave as Input
               assembly = protein.assembly[i,1])
      
      # Add neighbors to all neighbors data frame
      all.neighbours <- rbind(all.neighbours, np)}
    }

  # Get the current date in YYYY-MM-DD format
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  # Save results of this long run
  output_file <- file.path('output',current_date, paste('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv', sep = ""))
  write_csv(all.neighbours, output_file)

  fasta_output_file <- file.path('output',current_date, paste('all_neighbours_bp', basepairs, '_n', max_neighbors, '.fasta', sep = ""))
  data.frame2fasta(all.neighbours, name=c('ID','PIGI','assembly','product'), sequence='seq', output_file = fasta_output_file)
  
  return(all.neighbours)

}

concatenate_columns <- function(df, columns_vector, new_column_name = "concatenated") {
  # Concatenate specified columns into a single string separated by a comma
  # df: data frame
  # columns_vector: vector of column names to concatenate
  # new_column_name: name of the new column with concatenated values
  
  df %>%
    unite(new_column_name, all_of(columns_vector), sep = ",", remove = FALSE)
}

data.frame2fasta <- function(df, name_vector, sequence_column="seq", output_file="fasta.fasta"){
  # Create a fasta file from a data frame
  # df: data frame
  # name: column name for the name of the sequence
  # sequence: column name for the sequence
  # return: a character vector with the fasta file
  
  # Generate new df with a name column and a sequence column
  df2 <- concatenate_columns(df, name_vector, new_column_name = "name")

  df2 <- df2%>%
    select(all_of(c("name", sequence_column))) %>%
    mutate(name = paste(">", name))

  fasta_lines <- do.call(rbind, lapply(seq(nrow(df2)), function(i) t(df2[i, ])))
  write.table(fasta_lines, row.names = FALSE, col.names = FALSE, quote = FALSE, file = output_file)
}

#' Plot neighbors
#'
#' @param all.neighbours.df A data frame with all neighbors.
#' @param protein.id The protein ID to plot.
#' @export
plot.neighbours <- function(all.neighbours.df, protein.id){
  # Get assembly for given protein id
  assembly <- all.neighbours.df$assembly[all.neighbours.df$PIGI == protein.id][1]
  
  # Read gff data
  df <- read.gff(paste('data/ncbi_dataset/data/', assembly,'/genomic.gff', sep = ""), na.strings = c(".", "?"), GFF3 = TRUE)
  
  # Get Protein ID and Name separately
  df <- df %>%
    mutate(ID = gsub(".*\\-", "", getAttributeField(attributes, "ID")),  # with a little addition to get the clean id :)
           Name = getAttributeField(attributes, "Name"))
  
  # Create data frame for plotting
  df0 <- data.frame(molecule = df$seqid[df$ID==protein.id],
                    gene = protein.id, 
                    start = df$start[df$ID==protein.id],
                    end = df$end[df$ID==protein.id],
                    strand = df$strand[df$ID==protein.id],
                    orientation = 1)
  
  df1 <- all.neighbours.df %>%
    filter(PIGI == protein.id) %>%
    select(molecule = seqid, gene = product, start, end, strand, orientation = 1)
  
  df2 <- rbind(df0, df1)
  
  # Plot neighbors
  ggplot(df2, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
    geom_gene_arrow() +
    facet_wrap(~ molecule, scales = "free", ncol = 1) +
    theme_genes()
  ggsave(paste('neighbours_', protein.id, '.png', sep = ""), 
         plot = last_plot(), device = "png", path = "output", 
         width = 20, height = 5, units = "cm")
} 

#' Combine and plot neighbors
#'
#' @param neighbours_data A data frame with neighbors data.
#' @param cog_data A data frame with COG information from analyze_proteins_cog_classifier function.
#' @param clade_assign A data frame with clade assignments.
#' @return A combined data frame.
#' @export
combine_and_plot <- function(neighbours_data, cog_data, clade_assign, neighbour_annotations){
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  
  # Merge the fasta and the clade assignment
  if (is.null(clade_assign)) {
    warning("No clade assignment found. Skipping the combination.")
    neighbours_with_clades <- neighbours_data %>%
      mutate(clade = 'unknown')
  }else{
    # Merge the clade assignment
    neighbours_with_clades <- neighbours_data %>%
      left_join(clade_assign, by='PIGI')%>%
      replace_na(list(clade = 'unknown'))
  }

  # Combine all information into one dataframe if possible
  if (is.null(cog_data)) {
    warning("No additional COG data found. Skipping the combination.")
    combined_data <- neighbours_with_clades %>%
      mutate(COG_ID = NA,
              CDD_ID = NA,
              EVALUE = NA,
              GENE_NAME = NA,
              COG_NAME = 'unkown',
              COG_LETTER = 'unkown',
              COG_DESCRIPTION = 'unkown'
              )

  }else{
    # Merge the COG data
    combined_data <- neighbours_with_clades %>%
      left_join(cog_data, by = c("ID" = "QUERY_ID")) %>%
      replace_na(list(COG_ID = NA,
                      CDD_ID = NA,
                      EVALUE = NA,
                      GENE_NAME = NA,
                      COG_NAME = 'unkown',
                      COG_LETTER = 'unkown',
                      COG_DESCRIPTION = 'unkown'
                      ))
  } 

  if(!is.null(neighbour_annotations)){
    combined_data <- combined_data %>%
      left_join(neighbour_annotations, by = c("COG_NAME" = "COG_NAME")) %>%
      replace_na(list(annotation = 'unknown'))
  }

  # Write the combined data to a csv file
  write_csv(combined_data, file.path("output", current_date,"combined_df_all_neighbours_assigned.csv"))


  return(combined_data)
}

#' Plot neighbors per clade
#'
#' @param combined_data A combined data frame with neighbors and clade information.
#' @export
plot_neighbours_per_clade <- function(combined_data){
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  # Extract the amount of one type of neighbour per clade 
  neighbour_count_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER) %>%
    distinct() %>%
    count(clade, COG_LETTER, .drop = FALSE) %>%
    replace_na(list(clade = "unknown", COG_LETTER = "unknown")) %>%
    mutate(COG_LETTER = as.factor(COG_LETTER))
  write_csv(neighbour_count_per_clade, file.path("output", current_date, "neighbour_count_per_clade.csv"))
  
  # Extract how many representatives from each clade were involved
  representatives_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    replace_na(list(clade = "unknown")) %>%
    distinct(PIGI, clade) %>%
    count(clade, .drop = FALSE) %>%
    mutate(clade = as.factor(clade))
  write_csv(representatives_per_clade, file.path("output", current_date, "representatives_per_clade.csv"))
  
  # Extract how many neighbours from each clade were found
  total_neighbours_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    distinct(COG_LETTER, clade) %>%
    count(clade, .drop = FALSE) %>%
    mutate(across(where(is.factor), as.character)) %>%
    replace_na(list(clade = "unknown"))
  write_csv(total_neighbours_per_clade, file.path("output", current_date, "total_neighbours_per_clade.csv"))
  
  # Create facet label names for clade variable
  clade_labels <- paste(c(LETTERS[1:6], 'hcp', 'unkown clade'), ", n=", representatives_per_clade$n, sep="")
  names(clade_labels) <- c(LETTERS[1:7], 'unknown')
  
  
  plot_height <- round_any(max(neighbour_count_per_clade$n), 100, f = ceiling)
  
  # Ensure plot_height is a finite number
  if (!is.finite(plot_height) || plot_height <= 0) {
    warning("Calculated plot height is not a finite positive number. Using ggplot2 default.")
    plot_height <- NULL
  }

  # Plot the amount of neighbour of one sort per clade
  p <- ggplot(neighbour_count_per_clade, aes(x= COG_LETTER, y= n, fill=COG_LETTER)) + 
    # Add bars with identity statistic, position them side by side (dodge), and outline in black
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 1) +
    
    # Create facets based on the second variable, with free y scales and labels on both sides
    facet_grid(~clade, scales = "free_y", labeller = labeller(clade = clade_labels)) +
    
    
    # Use a minimal theme with base font size 12, and place the legend at the top
    theme_bw(base_size = 12) +
    
    # Remove panel grid, x axis title, x axis text, and x axis ticks
    # Customize y axis ticks and line, and remove strip background
    # Adjust spacing between panels
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          axis.ticks.length=unit(.2, "cm"),
          axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          strip.background =element_blank(),
          panel.spacing.x= unit(-0.1, "cm"),
          legend.position="right") +
    
    
    # Set x axis expansion to create space between axis and bars
    scale_x_discrete(expand=c(0,1)) +
    
    # Set y axis label and fill legend title
    labs(y = "Count per clade", fill = "Types of Neighbours") +
    
    # control relation to data an y axis basically you zoom in but this could be handy if you data is barely above the highest tick point, but still want to inclucde it to the plot? like in base r bar plot
    coord_cartesian(ylim = c(0, plot_height),
                    clip = 'off') +
    
    #generate text above the bars with the amount of neighbours
    geom_text(aes(label=n), hjust=-0.2, colour = 'black', size = 2, angle = 90) 
  # Set y axis breaks, limits, and expansion
  if (!is.null(plot_height)) {
    p <- p + scale_y_continuous(breaks = seq(0, plot_height, plot_height/5), 
                       limits = c(0,plot_height), 
                       expand = c(0,0)) 
  }
  # Print the plot
  print(p)

  # Generate color vector
  n <- 20
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Define color vector and print png
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave(
  filename = file.path("output", current_date, "amount_of_neighbour_per_clade.png"),
  width = 30, height = 15, units = "cm"
))
  
  # Define color vector and print svg
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave(
  filename = file.path("output", current_date, "amount_of_neighbour_per_clade.svg"),
  width = 30, height = 15, units = "cm"
))
}

#' Plot neighbors per clade without unknown clade
#'
#' @param combined_data A combined data frame with neighbors and clade information.
#' @export
plot_neighbours_per_clade2 <- function(combined_data){
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  # Extract the amount of one type of neighbour per clade 
  neighbour_count_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER) %>%
    distinct() %>%
    count(clade, COG_LETTER, .drop = FALSE) %>%
    replace_na(list(clade = "unknown", COG_LETTER = "unknown")) %>%
    mutate(COG_LETTER = as.factor(COG_LETTER))
  
  # Extract how many representatives from each clade were involved
  representatives_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    replace_na(list(clade = "unknown")) %>%
    distinct(PIGI, clade) %>%
    count(clade, .drop = FALSE) %>%
    mutate(clade = as.factor(clade))
  
  # Extract how many neighbours from each clade were found
  total_neighbours_per_clade <- combined_data %>%
    select(ID, clade, COG_LETTER, PIGI) %>%
    distinct(COG_LETTER, clade) %>%
    count(clade, .drop = FALSE) %>%
    mutate(across(where(is.factor), as.character)) %>%
    replace_na(list(clade = "unknown"))
  
  
  # Create facet label names for clade variable
  clade_labels <- paste(c(LETTERS[1:6],'hcp', 'unkown clade'), ", n = ", representatives_per_clade$n, sep="")
  names(clade_labels) <- c(LETTERS[1:7], 'unknown')
  
  
  
  # Plot the amount of neighbour of one sort per clade without unknown clade
  neighbour_count_per_clade_no_unkown <- neighbour_count_per_clade %>%
    filter(clade != 'unknown') %>%
    filter(COG_LETTER != 'unknown')
  
  if(nrow(neighbour_count_per_clade_no_unkown) == 0){
    warning("There are only unknown clades and/or neighbours. Skipping the plot.")
    return()
  }
  plot_height <- round_any(max(neighbour_count_per_clade$n), 100, f = ceiling)
  
  # Ensure plot_height is a finite number
  if (!is.finite(plot_height) || plot_height <= 0) {
    warning("Calculated plot height is not a finite positive number. Using ggplot2 default.")
    plot_height <- NULL
  }
  # Plot the amount of neighbour of one sort per clade
  p<-ggplot(neighbour_count_per_clade_no_unkown, aes(x= COG_LETTER, y= n, fill=COG_LETTER)) +
    # Add bars with identity statistic, position them side by side (dodge), and outline in black
    geom_bar(stat = "identity", position = "dodge", color = "black", width = 1) +
    
    # Create facets based on the second variable, with free y scales and labels on both sides
    facet_grid(~clade, scales = "free_y", labeller = labeller(clade = clade_labels)) +
    
    # Use a Brewer color palette for the fill aesthetic
   # scale_fill_brewer(palette = "Pastel2") +
    
    # Use a minimal theme with base font size 12, and place the legend at the top
    theme_bw(base_size = 12) +
    
    # Remove panel grid, x axis title, x axis text, and x axis ticks
    # Customize y axis ticks and line, and remove strip background
    # Adjust spacing between panels
    theme(panel.grid = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          axis.ticks.length=unit(.2, "cm"),
          axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
          strip.background =element_blank(),
          panel.spacing.x= unit(-0.1, "cm"),
          legend.position="right") +
    
    # Set x axis expansion to create space between axis and bars
    scale_x_discrete(expand=c(0,1)) +
    
    # Set y axis label and fill legend title
    labs(y = "Count per clade", fill = "Types of Neighbours") +
    
    # control relation to data an y axis basically you zoom in but this could be handy if you data is barely above the highest tick point, but still want to inclucde it to the plot? like in base r bar plot
    coord_cartesian(ylim = c(0, plot_height),
                    clip = 'off') +
    
    #generate text above the bars with the amount of neighbours
    geom_text(aes(label=n), hjust=-0.2, colour = 'black', size = 2, angle = 90) 
    # Set y axis breaks, limits, and expansion

  if (!is.null(plot_height)) {
    p <- p + scale_y_continuous(breaks = seq(0, plot_height, plot_height/5), 
                       limits = c(0,plot_height), 
                       expand = c(0,0)) 
  }
  # Print the plot
  print(p)
  
  # Generate color vector
  n <- 20
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Define color vector and print png
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave(file.path("output",current_date,"amount_of_neighbour_per_cladenounknown.png"),
           width = 30, height = 15, units = "cm"))
  
  # Define color vector and print svg
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave(file.path("output",current_date,"amount_of_neighbour_per_cladenounknown.svg"),
         width = 30, height = 15, units = "cm"))
}


#' Create or Validate Output File Path
#'
#' This function prompts the user to enter a date and checks if the corresponding directory exists.
#' If the directory does not exist, it prompts the user with options to rerun, change the date, or exit the program.
#' The function returns the output file path based on the provided date or the current date.
#'
#' @param basepairs An integer representing the number of base pairs.
#' @param max_neighbors An integer representing the maximum number of neighbors.
#' @param date A character string representing the date of the data to load in format YYYY-MM-DD (default is NULL).
#' @return A character string representing the output file path.
#' @examples
#' \dontrun{
#' create_or_validate_output_file_path(100, 5, "2022-01-01")
#' }
create_or_validate_output_file_path <- function(basepairs, max_neighbors, date = NULL) {
    # Get the current date in YYYY-MM-DD format
    current_date <- format(Sys.Date(), "%Y-%m-%d")
    if (is.null(date)) {
  repeat {
    
    # Ask user if they want to reload the data from an older date that they specify
    date <- readline(prompt = "Enter the date of the data you want to load (YYYY-MM-DD): ")
    
    
    # Determine the output file path based on the provided date or current date
    if (date != "") {
      output_dir <- file.path('output', date)
      if (dir.exists(output_dir)) {
        output_file <- file.path(output_dir, paste('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv', sep = ""))
        break
      } else {
        cat("There is no valid folder for your date.\n")
        choice <- readline(prompt = "Do you want to do a new (R)un, (C)hange date, or (E)xit the program? ")
        if (toupper(choice) == "E") {
          stop("Exiting the program.")
        } else if (toupper(choice) == "C") {
          next
        } else if (toupper(choice) == "R") {
          date <- ""
        }
      }
    }
    
    if (date == "") {
      output_dir <- file.path('output', current_date)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      output_file <- file.path(output_dir, paste('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv', sep = ""))
      break
    }
  }
  }else{
    output_dir <- file.path('output', date)
    if (!dir.exists(output_dir)) {
      stop("There is no valid folder for your date supplied by main().")
    }else{
      output_file <- file.path(output_dir, paste('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv', sep = ""))
    }}
  return(output_file)
}


#' Analyze neighbors
#' 
#' This function will dowload the sequence of the neighbours an subsequently analyse them via the COGclassifier.
#' 
#' @param neighbours_data A data frame with neighbors data.
#' @param output_dir The output directory.
#' 
#' @return A data frame with the COG classification.
#' @export
#' 
analyze_proteins_cog_classifier <- function(df, column= 'ID') {
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  output_dir <- paste("output/", current_date, sep = "") 
  # Export accessions of the neighbours  
  neighbours_accesions <- data.frame(df%>%select(all_of(column)))
  write_csv(neighbours_accesions, file.path(output_dir, "neighbours_accesions.csv"), col_names = FALSE)

  # Get the fasta file of the neighbours
  bash_command <- paste("./get_fasta_from_accession.sh", file.path(output_dir, "neighbours_accesions.csv"), file.path(output_dir, "neighbours.fasta"))
  cat("Running command:", bash_command, "\n")
  system2(bash_command, wait = TRUE)
  
  # Run COGclassifier
  cog_command <- paste("COGclassifier -i", file.path(output_dir, "neighbours.fasta"), "-o", file.path(output_dir, "cogclassifier"))
  cat("Running command:", cog_command, "\n")
  system2(cog_command, wait = TRUE)

  cog_classification <- read.delim(file.path(output_dir, "cogclassifier/classifier_result.tsv"))


  return(cog_classification)
} 



