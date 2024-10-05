# Protein for GFF analysis ####

#' Read Protein and Assembly Data
#'
#' This function reads the protein and assembly data from specified files.
#'
#' @param protein_file The path to the protein file (default is 'data/proteins.txt'). This file should contain all proteins IDs that you are interested in, separated by a new line.
#' @param assembly_file The path to the assembly file (default is 'data/assm_accs.txt') This file should contain all genome assembly numbers that you want to include in your analysis, separeted by a new line.
#' @param protein_assembly_file The path to the protein assembly file (default is 'data/assm_accs_proteins.txt').
#' @return A list containing the protein, assembly, and protein assembly data frames.
#' @export
read_protein_assembly_data <- function(protein_file = 'proteins.csv', 
                                       assembly_file = 'assm_accs.csv', 
                                       protein_assembly_file = 'assm_accs_protein.csv',
                                       PATH = 'data'){
  protein_of_interest <- readline(prompt = "Enter the Accession Number of protein of interest: ")
  protein <- read.csv(file.path(PATH,protein_file), header = FALSE)
  assembly <- read.csv(file.path(PATH,assembly_file), header = FALSE)
  protein_assembly <- read.csv(file.path(PATH,protein_assembly_file), header = FALSE)
  
  return(list(protein_of_interest = protein_of_interest, 
              protein = protein, 
              assembly = assembly, 
              protein_assembly = protein_assembly))
}
#' Read clades
#'
#' This function reads clade information from text files and processes them.
#'
#' @param PATH The base path to the data directory (default is 'data').
#' @param clade_dir The directory containing clade files (default is 'clades').
#' @param pattern The pattern to match clade files (default is "Clade*.txt").
#' @return A data frame with clade assignments, or an empty data frame if no clade files are found.
#' @export
read_clades <- function(PATH = 'data', clade_dir = 'clades', pattern = "^[C]") {
  # Load clade information from text files
  PROTEIN_ALIAS <- read_representatives(PATH = PATH)
  clade_path <- file.path(PATH, clade_dir)
  clade_files <- list.files(clade_path, pattern = pattern, full.names = TRUE)
  
  # Check if any clade files are found
  if (length(clade_files) == 0) {
    warning("No clade files found. Returning an empty data frame.")
    return(data.frame(protein.id = character(), clade = character(), PIGI = character(), stringsAsFactors = FALSE))
  }
    
  # Apply the function to each file and clade, then bind rows
  clade_assign <- map2_df(clade_files, LETTERS[1:length(clade_files)], process_clade_file)
  clade_assign$PIGI <- unlist(flatten(lapply(clade_assign$protein.id, function(i) { protein.alias(i, alias = PROTEIN_ALIAS, verbose = FALSE)[1, 1] })))
  
  return(clade_assign)
}

#' Process clade file
#'
#' This function reads and processes a single clade file.
#'
#' @param file The path to the clade file.
#' @param clade The clade identifier.
#' @return A data frame with processed clade information.
#' @export
process_clade_file <- function(file, clade) {
  read_delim(file, col_names = c('protein.id','V2','V3','V4'), delim="/", show_col_types = FALSE) %>%
    select(protein.id) %>%
    mutate(clade = clade)
}


#' Get the Types of Neighbours and Their Counts
#'
#' This function calculates the types of neighbours and their counts, and generates a plot and CSV file.
#' If no cluster domain assignments are found, the function uses protein product information for NCBI annotation.
#'
#' @param cog_data Dataframe with cog_data.
#' @export
amount_of_neighbours <- function(cog_data) {
  current_date <- format(Sys.Date(), "%Y-%m-%d")

    types_of_neighbours <- cog_data %>%
    select(COG_NAME, COG_LETTER) %>%
    add_count(COG_NAME) %>%
    arrange(desc(n))%>%
    unique()
  # Plot the types of neighbours in a range
  ggplot(types_of_neighbours[c(1:100),], 
         aes(x = reorder(COG_LETTER, -n), y = n)) +
    geom_bar(stat = "identity")
  ggsave(file.path('output',current_date,'types_of_neighbours.png'))
  write_csv(types_of_neighbours, file.path('output',current_date,'types_of_neighbours.csv'))
 
  print('Please open the output folder to see the types of neighbours and annotated them as per the instructions in the README.md file (optional)')
}

#' Read Neighbour Annotations
#'
#' This function re-imports the types of neighbours after manual annotation.
#'
#' @param current_date A string representing the current date, used to locate the annotated CSV file.
#' @param empty_cells A string representing the value to replace empty cells with. Default is "unknown".
#' @return A data frame containing the annotated neighbours if the file exists, otherwise NULL.
#' @importFrom readr read_csv
#' @importFrom dplyr mutate_all rename
#' @importFrom tidyr replace_na
#' @export
read_annotations <- function(current_date, empty_cells = "unknown") {
  # Re-import the types of neighbours after manual annotation
  if(file.exists(file.path("output", current_date, "types_of_neighbours_annotated.csv"))){
    annotated_neighbours <- read_csv(file.path("output", current_date, "types_of_neighbours_annotated.csv"), show_col_types = FALSE)
    
    # Replace empty cells with NA
    annotated_neighbours <- annotated_neighbours %>%
      mutate_all(~ifelse(. == "", NA, .))
    
    # Replace NA with the specified value
    if (!is.null(empty_cells)) {
      annotated_neighbours <- annotated_neighbours %>%
        replace_na(list(COG_NAME = empty_cells, COG_LETTER= empty_cells, N = empty_cells, ANNOTATION = empty_cells))
    }
    
    # Rename columns
    colnames(annotated_neighbours)[1:4] <- c("COG_NAME", "COG_LETTER","N", "ANNOTATION")
  } else {
    annotated_neighbours <- NULL
  }
  return(annotated_neighbours)
}


# Protein Alias ####

#' Read Representatives
#'
#' This function reads the IPG, PDB, and cluster alias data from text files and combines them.
#'
#' @param PATH The base path to the data directory (default is 'data').
#' @param path The subdirectory containing the alias files (default is 'representatives').
#' @param ipg_file The name of the IPG alias file (default is 'ipg_representative.txt').
#' @param pdb_file The name of the PDB alias file (default is 'pdb_representative.txt').
#' @param cluster_file The name of the cluster alias file (default is 'cluster_representative.txt').
#' @return A data frame with combined alias data, or an empty data frame if no alias files are found.
#' @export
read_representatives <- function(PATH = 'data',
                                 path = 'representatives', 
                                 ipg_file = 'ipg_representative.txt', 
                                 pdb_file = 'pdb_representative.txt', 
                                 cluster_file = 'cluster_representative.txt') {
    
  path <- file.path(PATH, path)
  
  # Helper function to read alias files
  read_alias_file <- function(file_name, col_names) {
    file_path <- file.path(path, file_name)
    if (file.exists(file_path)) {
      read_delim(file_path, col_names = col_names, show_col_types = FALSE)
    } else {
      warning(paste("File not found:", file_path))
      return(NULL)
    }
  }
  
  # Read the IPG alias data
  ipg_alias <- read_alias_file(ipg_file, c('PIGI', 'alias'))
  if (!is.null(ipg_alias)) ipg_alias$identity <- 1
  
  # Read the PDB alias data
  pdb_alias <- read_alias_file(pdb_file, c('alias', 'PIGI'))
  if (!is.null(pdb_alias)) pdb_alias$identity <- 1
  
  # Read the cluster alias data
  cluster_alias <- read_alias_file(cluster_file, c( 'PIGI', 'alias', 'identity'))
  if (!is.null(cluster_alias)) {
    cluster_alias <- cluster_alias %>%
      mutate(identity = as.numeric(sub("%", "", identity)) / 100)
  }
  
  # Combine the alias data
  combined_alias <- bind_rows(ipg_alias, pdb_alias, cluster_alias)
  
  # Check if combined_alias is empty
  if (nrow(combined_alias) == 0) {
    warning("No alias data found. Returning an empty data frame.")
    return(data.frame(alias = character(), PIGI = character(), stringsAsFactors = FALSE))
  }
  
  # Replace PIGI values based on PDB alias replacements
  if (!is.null(pdb_alias)) {
    replacements <- setNames(pdb_alias$PIGI, pdb_alias$alias)
    combined_alias$PIGI <- paste0(sub("\\..*", "", str_replace_all(combined_alias$PIGI, replacements)), ".1")
  }
  
  # Remove duplicate entries
  combined_alias <- combined_alias %>%
    distinct(alias, PIGI)
  
  return(combined_alias)
}

