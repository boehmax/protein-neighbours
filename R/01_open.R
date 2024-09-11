# Protein for GFF analysis ####

#' Read Protein and Assembly Data
#'
#' This function reads the protein and assembly data from specified files.
#'
#' @param protein_file The path to the protein file (default is 'data/proteins.txt').
#' @param assembly_file The path to the assembly file (default is 'data/assm_accs.txt').
#' @param protein_assembly_file The path to the protein assembly file (default is 'data/assm_accs_proteins.txt').
#' @return A list containing the protein, assembly, and protein assembly data frames.
#' @export
read_protein_assembly_data <- function(protein_file = 'proteins.txt', 
                                       assembly_file = 'assm_accs.txt', 
                                       protein_assembly_file = 'assm_accs_proteins.txt',
                                       PATH = 'data'){
  protein_of_interest <- readline(prompt = "Enter the Accession Number of protein of interest: ")
  protein <- read.table(protein_file, header = FALSE)
  assembly <- read.table(assembly_file, header = FALSE)
  protein_assembly <- read.table(protein_assembly_file, header = FALSE)
  
  return(list(protein_of_interest = protein_of_interest, 
              protein = protein, 
              assembly = assembly, 
              protein_assembly = protein_assembly))
}
#' Read clades
#'
#' This function reads clade information from text files and processes them.
#'
#' @param clade_dir The directory containing clade files.
#' @param pattern The pattern to match clade files.
#' @return A data frame with clade assignments.
#' @export
read_clades <- function(PATH = 'data', clade_dir = 'clades', pattern = "Clade.*_20231130.txt"){
  # Load clade information from text files
  clade_path <- file.path(PATH, clade_dir)
  clade_files <- list.files(clade_path, pattern = "Clade.*_20231130.txt", full.names = TRUE)
  
  # Define a function to read and process each file
  process_clade_file <- function(file, clade) {
    read_delim(file, col_names = c('protein.id','V2','V3','V4'), delim="/", show_col_types = FALSE) %>%
      select(protein.id) %>%
      mutate(clade = clade)
  }
  
  # Apply the function to each file and clade, then bind rows
  clade_assign <- map2_df(clade_files, LETTERS[1:length(clade_files)], process_clade_file)
  clade_assign$PIGI <- unlist(flatten(lapply(clade_assign$protein.id, function(i){protein.alias(i, verbose = FALSE)[1,1]})))
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
#' Read Cluster Domain Information
#'
#' This function reads and processes cluster domain information from text files.
#'
#' @param classification_path The path to the classification directory (default is 'data/classification/COG').
#' @return A data frame with cluster domain assignments.
#' @export
read_cluster_domain <- function(PATH = 'data', classification_path = file.path('classification', 'COG')){
  # Load Cluster Domain information from text files
  classification_path <- file.path(PATH, classification_path)
  cd_files <- list.files(classification_path, pattern = "cog_.*hitdata.txt", full.names = TRUE)
  
  # Define a function to read and process each file
  process_cd_files <- function(file) {
    read_delim(file, col_names = c('Query',	'Hit type',	'PSSM-ID',	'From',	'To',	'E-Value',	'Bitscore',	'Accession',	'Short name',	'Incomplete',	'Superfamily'), skip= 8, delim="\t", show_col_types = FALSE) %>%
      mutate(ID = gsub(".*-","",Query)) %>%
      mutate(across(c('Superfamily','Short name','Accession'),as.factor)) %>%
      mutate(`Short name` = sub(" superfamily","",`Short name`))%>%
      mutate(ID = sub(" ","",ID))
  }
  
  # Apply the function to each file and clade, then bind rows
  cd_assign <- map_df(cd_files, process_cd_files)
  return(cd_assign)
}

#' Get the Types of Neighbours and Their Counts
#'
#' This function calculates the types of neighbours and their counts, and generates a plot and CSV file.
#'
#' @param cd_assign A data frame with cluster domain assignments.
#' @export
amount_of_neighbours <- function(cd_assign){
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  types_of_neighbours <- cd_assign %>%
    select(`Short name`) %>%
    add_count(`Short name`)  %>%
    arrange(desc(n))%>%
    unique()
  # Plot the types of neighbours in a range
  ggplot(types_of_neighbours[c(1:100),], 
         aes(x = reorder(`Short name`, -n), y = n)) +
    geom_bar(stat = "identity")
  ggsave(file.path('output',current_date,'types_of_neighbours.png'))
  write_csv(types_of_neighbours, file.path('output',current_date,'types_of_neighbours.csv'))
  print('Please open the output folder to see the types of neighbours and annotated them as per the instructions in the README.md file')
}

# Protein Alias ####

#' Read Representatives
#'
#' This function reads the IPG, PDB, and cluster alias data from text files and combines them.
#'
#' @param path The base path to the data directory (default is 'data/representatives').
#' @param ipg_file The name of the IPG alias file (default is 'ipg_representative.txt').
#' @param pdb_file The name of the PDB alias file (default is 'pdb_representative.txt').
#' @param cluster_file The name of the cluster alias file (default is 'cluster_representative.txt').
#' @return A data frame with combined alias data.
#' @export
read_representatives <- function(PATH = 'data',
                 path = 'representatives', 
                 ipg_file = 'ipg_representative.txt', 
                 pdb_file = 'pdb_representative.txt', 
                 cluster_file = 'cluster_representative.txt') {
    
  path <- file.path(PATH, path)
  # Helper function to read alias files
  read_alias_file <- function(file_name, col_names) {
  read_delim(file.path(path, file_name), col_names = col_names, show_col_types = FALSE)
  }
  
  # Read the IPG alias data
  ipg_alias <- read_alias_file(ipg_file, c('alias', 'PIGI'))
  ipg_alias$identity <- 1
  
  # Read the PDB alias data
  pdb_alias <- read_alias_file(pdb_file, c('alias', 'PIGI'))
  pdb_alias$identity <- 1
  
  # Read the cluster alias data
  cluster_alias <- read_alias_file(cluster_file, c('alias', 'PIGI', 'identity'))
  cluster_alias <- cluster_alias %>%
  mutate(identity = as.numeric(sub("%", "", identity)) / 100)
  
  # Combine the alias data
  combined_alias <- bind_rows(ipg_alias, pdb_alias, cluster_alias)
  
  # Replace PIGI values based on PDB alias replacements
  replacements <- setNames(pdb_alias$PIGI, pdb_alias$alias)
  combined_alias$PIGI <- paste0(sub("\\..*", "", str_replace_all(combined_alias$PIGI, replacements)), ".1")
  
  # Remove duplicate entries
  combined_alias <- combined_alias %>%
  distinct(alias, PIGI)
  
  return(combined_alias)
}

