#rounding function
round_any <- function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

# Function to extract column from attribute column
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

# Define a function to find the protein alias
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

# Function to get neighboring proteins
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

# Combined function, requires ape
getProteinNeighborsFromGff3 <- function(input, protein.id, basepairs = 300, m = 15){
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

# Loop through each protein
collec.all.neigbour <- function(protein.assembly, basepairs = 300, m = 15){
  all.neighbours <- data.frame() # Initialize empty data frame for all neighbors
  for(i in seq_len(nrow(protein.assembly))){
    # Get neighbors for each protein
    np <- getProteinNeighborsFromGff3(paste('data/ncbi_dataset/data/', protein.assembly[i,1],'/genomic.gff', sep = ""), 
                                      protein.assembly[i,2], basepairs, m)
    if(nrow(np) > 0){
      # Add protein and assembly information to neighbors data frame
      np <- np %>%
        mutate(PIGI = protein.assembly[i,2],  # Protein I Gave as Input
               assembly = protein.assembly[i,1])
      
      # Add neighbors to all neighbors data frame
      all.neighbours <- rbind(all.neighbours, np)}
    # save resluts of this long run
    write_csv(all.neighbours, paste('output/all_neighbours_bp',basepairs, '_n', m,'.csv', sep = ""))
  }
}

# Function to plot neighbors
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


combine_and_plot <- function(neighbours_data, cd_data, neighbour_types, clade_assign){
  # Merge the fasta and the clade assignment
  neighbours_with_clades <- neighbours_data %>%
    left_join(clade_assign, by='PIGI')
  # Combine all information into one dataframe
  combined_data <- neighbours_with_clades %>%
    left_join(cd_data, by = 'ID') %>%
    left_join(neighbour_types, by = 'Short name')
  return(combined_data)
  # Save the combined data frame
  write_csv(combined_data, "output/combined_df_all_neighbours_assigned.csv")
}

plot_neighbours_per_clade <- function(combined_data){
  # Extract the amount of one type of neighbour per clade 
  neighbour_count_per_clade <- combined_data %>%
    select(ID, clade, type.y) %>%
    distinct() %>%
    count(clade, type.y, .drop=FALSE)%>%
    replace_na(list(clade = "unknown_clade", type.y = "unknown")) %>%
    mutate(y.type = as.factor(type.y))
  
  # Extract how many representatives from each clade were involved
  representatives_per_clade <- combined_data %>%
    select(ID, clade, type.y, PIGI) %>%
    distinct(PIGI, clade) %>%
    count(clade, .drop=FALSE)
  
  # Extract how many neighbours from each clade were found
  neighbours_per_clade <- combined_data %>%
    select(ID, clade, type.y, PIGI) %>%
    distinct(type.y, clade)  %>%
    count(clade, .drop=FALSE)%>%
    mutate_if(is.factor,as.character)%>%
    replace_na(list(clade = "unknown_clade"))
  
  
  # Create facet label names for clade variable
  clade_labels <- paste(c(LETTERS[1:6], 'unkown clade'), ", n=", representatives_per_clade$n, sep="")
  names(clade_labels) <- c(LETTERS[1:6], 'unknown_clade')
  
  
  plot_height <- round_any(max(neighbour_count_per_clade$n), 100, f = ceiling)
  # Plot the amount of neighbour of one sort per clade
  ggplot(neighbour_count_per_clade, aes(x= type.y, y= n, fill=type.y)) + 
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
          panel.spacing.x= unit(-0.1, "cm")) +
    
    # Set y axis breaks, limits, and expansion
    scale_y_continuous(breaks = seq(0, plot_height, plot_height/5), 
                       limits = c(0,plot_height), 
                       expand = c(0,0)) +
    
    # Set x axis expansion to create space between axis and bars
    scale_x_discrete(expand=c(0,1)) +
    
    # Set y axis label and fill legend title
    labs(y = "Count per clade", fill = "Types of Neighbours") +
    
    # control relation to data an y axis basically you zoom in but this could be handy if you data is barely above the highest tick point, but still want to inclucde it to the plot? like in base r bar plot
    coord_cartesian(ylim = c(0, plot_height),
                    clip = 'off') +
    
    #generate text above the bars with the amount of neighbours
    geom_text(aes(label=n), hjust=-0.2, colour = 'black', size = 2, angle = 90) 
  # Generate color vector
  n <- 20
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Define color vector and print png
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave("output/amount_of_neighbour_per_clade.png",
           width = 30, height = 15, units = "cm")) 
  
  # Define color vector and print svg
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave("output/amount_of_neighbour_per_clade.svg",
           width = 30, height = 15, units = "cm"))
}

plot_neighbours_per_clade2 <- function(combined_data){
  # Extract the amount of one type of neighbour per clade 
  neighbour_count_per_clade <- combined_data %>%
    select(ID, clade, type.y) %>%
    distinct() %>%
    count(clade, type.y, .drop=FALSE)%>%
    replace_na(list(clade = "unknown_clade", type.y = "unknown"))%>%
    mutate(y.type = as.factor(type.y))
  
  # Extract how many representatives from each clade were involved
  representatives_per_clade <- combined_data %>%
    select(ID, clade, type.y, PIGI) %>%
    distinct(PIGI, clade) %>%
    count(clade, .drop=FALSE)
  
  # Extract how many neighbours from each clade were found
  neighbours_per_clade <- combined_data %>%
    select(ID, clade, type.y, PIGI) %>%
    distinct(type.y, clade)  %>%
    count(clade, .drop=FALSE)%>%
    mutate_if(is.factor,as.character)%>%
    replace_na(list(clade = "unknown_clade"))
  
  
  # Create facet label names for clade variable
  clade_labels <- paste(c(LETTERS[1:6], 'unkown clade'), ", n = ", representatives_per_clade$n, sep="")
  names(clade_labels) <- c(LETTERS[1:6], 'unknown_clade')
  
  
  
  # Plot the amount of neighbour of one sort per clade without unknown clade
  neighbour_count_per_clade_no_unkown <- neighbour_count_per_clade %>%
    filter(clade != 'unknown_clade') %>%
    filter(type.y != 'unknown')
  
  plot_height <- round_any(max(neighbour_count_per_clade_no_unkown$n), 100, f = ceiling)
  
  # Plot the amount of neighbour of one sort per clade
  ggplot(neighbour_count_per_clade_no_unkown, aes(x= type.y, y= n, fill=type.y)) +
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
          panel.spacing.x= unit(-0.1, "cm")) +
    
    # Set y axis breaks, limits, and expansion
    scale_y_continuous(breaks = seq(0, plot_height, plot_height/5), 
                       limits = c(0,plot_height), 
                       expand = c(0,0)) +
    
    # Set x axis expansion to create space between axis and bars
    scale_x_discrete(expand=c(0,1)) +
    
    # Set y axis label and fill legend title
    labs(y = "Count per clade", fill = "Types of Neighbours") +
    
    # control relation to data an y axis basically you zoom in but this could be handy if you data is barely above the highest tick point, but still want to inclucde it to the plot? like in base r bar plot
    coord_cartesian(ylim = c(0, plot_height),
                    clip = 'off') +
    
    #generate text above the bars with the amount of neighbours
    geom_text(aes(label=n), hjust=-0.2, colour = 'black', size = 2, angle = 90) 
  # Generate color vector
  n <- 20
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # Define color vector and print png
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave("output/amount_of_neighbour_per_cladenounknown.png",
           width = 30, height = 15, units = "cm"))
  
  # Define color vector and print svg
  withr::with_options(
    list(ggplot2.discrete.fill = col_vector),
    ggsave("output/amount_of_neighbour_per_cladenounknown.svg",
         width = 30, height = 15, units = "cm"))
}

