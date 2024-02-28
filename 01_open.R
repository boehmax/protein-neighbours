# Protein for GFF analysis ####
# Read protein and assembly data
protein_of_interest <- readline(prompt="Enter the Accesion Number of protein of interest: ")
protein <- read.table('data/proteins.txt', header = FALSE)
assembly <- read.table('data/assm_accs.txt', header = FALSE)
protein.assembly <- read.table('data/assm_accs_proteins.txt', header = FALSE)

# Read files for neighbour classification and clades####
read_clades <- function(){
  # Load clade information from text files
  clade_files <- list.files('data/clades', pattern = "Clade.*_20231130.txt", full.names = TRUE)
  
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

# Define a function to read and process each file
process_clade_file <- function(file, clade) {
  read_delim(file, col_names = c('protein.id','V2','V3','V4'), delim="/", show_col_types = FALSE) %>%
    select(protein.id) %>%
    mutate(clade = clade)
}

classification_path <- file.path('data', 'classification', 'COG')
read_cluster_domain <- function(){
  # Load Cluster Domain information from text files
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

# Get the types of neighbours and their counts
amount_of_neighbours <- function(cd_assign){
  types_of_neighbours <- cd_assign %>%
    select(`Short name`) %>%
    add_count(`Short name`)  %>%
    arrange(desc(n))%>%
    unique()
  # Plot the types of neighbours in a range
  ggplot(types_of_neighbours[c(1:100),], 
         aes(x = reorder(`Short name`, -n), y = n)) +
    geom_bar(stat = "identity")
  ggsave("output/types_of_neighbours.png")
  write_csv(types_of_neighbours, "output/types_of_neighbours.csv")
  print('Please open the output folder to see the types of neighbours and annotated them as per the instructions in the README.md file')
}

# Protein Alias ####
# Read the IPG alias data from the text file
read_represatatives <- function(){
  ipg.alias <- read_delim('data/representatives/ipg_representative.txt', col_names = c('alias','PIGI'))
  # Set the identity column to '1'
  ipg.alias$identity <- '1'
  
  # Read the PDB alias data from the text file
  pdb.alias <- read_delim('data/representatives/pdb_representative.txt', col_names = c('alias','PIGI'))
  replacements <- setNames(pdb.alias$PIGI, pdb.alias$alias)
  # Set the identity column to '1'
  pdb.alias$identity <- '1'
  
  # Read the cluster alias data from the text file
  cluster.alias <- read_delim('data/representatives/cluster_representative.txt', col_names = c('alias','PIGI','identity'))
  # Convert the identity column to numeric values
  cluster.alias <- cluster.alias %>%
    mutate(identity = as.numeric(sub("%","",identity))/100)
  
  # Combine the IPG and cluster alias data
  PIGI.alias <- rbind(ipg.alias, cluster.alias,pdb.alias)
  PIGI.alias$PIGI <- paste(sub("\\..*","",str_replace_all(PIGI.alias$PIGI, replacements)),".1", sep="") 
  
  PIGI.alias <- PIGI.alias %>%
    distinct(alias, PIGI)
  
  return(PIGI.alias)
  }