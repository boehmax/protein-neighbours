#' Input/Output Functions for Protein Neighbor Analysis
#'
#' This file contains functions for reading and writing data files,
#' including protein and assembly data, clades, and representative proteins.
#'
#' @importFrom dplyr %>%
#' @author Maximilian Böhm
#' @modified Improved version with enhanced documentation and error handling

#' Read Protein and Assembly Data
#'
#' This function reads the protein and assembly data from specified files
#' and prompts for a protein of interest.
#'
#' @param protein_file The path to the protein file (CSV with protein IDs).
#' @param assembly_file The path to the assembly file (CSV with assembly IDs).
#' @param protein_assembly_file The path to the protein assembly file (maps proteins to assemblies).
#' @param PATH The base path to the data directory.
#' @param interactive Whether to prompt the user for input (default TRUE).
#' @param protein_of_interest Optional protein ID to use (overrides prompt).
#' @return A list containing the protein, assembly, and protein assembly data frames.
#' @export
read_protein_assembly_data <- function(protein_file = "proteins.csv", 
                                       assembly_file = "assm_accs.csv", 
                                       protein_assembly_file = "assm_accs_protein.csv",
                                       PATH = "data",
                                       interactive = TRUE,
                                       protein_of_interest = NULL) {
  # Log reading of data
  pn_info("Reading protein and assembly data from:", PATH)

  # Construct full paths
  protein_path <- file.path(PATH, protein_file)
  assembly_path <- file.path(PATH, assembly_file)
  protein_assembly_path <- file.path(PATH, protein_assembly_file)

  # Check if files exist
  if (!file.exists(protein_path)) {
    pn_error("Protein file not found:", protein_path)
    stop(paste("Protein file not found:", protein_path))
  }

  if (!file.exists(assembly_path)) {
    pn_error("Assembly file not found:", assembly_path)
    stop(paste("Assembly file not found:", assembly_path))
  }

  if (!file.exists(protein_assembly_path)) {
    pn_error("Protein-assembly mapping file not found:", protein_assembly_path)
    stop(paste("Protein-assembly mapping file not found:", protein_assembly_path))
  }

  # Read files
  tryCatch({
    protein <- read.csv(protein_path, header = FALSE, stringsAsFactors = FALSE)
    pn_info("Read", nrow(protein), "proteins from", protein_path)
  }, error = function(e) {
    pn_error("Failed to read protein file:", e$message)
    stop(paste("Failed to read protein file:", e$message))
  })

  tryCatch({
    assembly <- read.csv(assembly_path, header = FALSE, stringsAsFactors = FALSE)
    pn_info("Read", nrow(assembly), "assemblies from", assembly_path)
  }, error = function(e) {
    pn_error("Failed to read assembly file:", e$message)
    stop(paste("Failed to read assembly file:", e$message))
  })

  tryCatch({
    protein_assembly <- read.csv(protein_assembly_path, header = FALSE, stringsAsFactors = FALSE)
    pn_info("Read", nrow(protein_assembly), "protein-assembly mappings from", protein_assembly_path)
  }, error = function(e) {
    pn_error("Failed to read protein-assembly file:", e$message)
    stop(paste("Failed to read protein-assembly file:", e$message))
  })

  # Get protein of interest
  if (is.null(protein_of_interest) && interactive) {
    protein_of_interest <- readline(prompt = "Enter the Accession Number of protein of interest: ")
  } else if (is.null(protein_of_interest)) {
    protein_of_interest <- ""
    pn_info("No protein of interest specified")
  } else {
    pn_info("Using specified protein of interest:", protein_of_interest)
  }

  return(list(protein_of_interest = protein_of_interest,
              protein = protein,
              assembly = assembly,
              protein_assembly = protein_assembly))
}

#' Read clades
#'
#' This function reads clade information from text files and processes them.
#'
#' @param PATH The base path to the data directory.
#' @param clade_dir The directory containing clade files.
#' @param pattern The pattern to match clade files.
#' @return A data frame with clade assignments, or an empty data frame if no clade files are found.
#' @export
read_clades <- function(PATH = "data", clade_dir = "clades", pattern = "^[C]") {
  # Log start of clade reading
  pn_info("Reading clade information from:", file.path(PATH, clade_dir))
  pn_info("Using pattern:", pattern)

  # Load protein alias data
  PROTEIN_ALIAS <- read_representatives(PATH = PATH)
  clade_path <- file.path(PATH, clade_dir)

  # Check if clade directory exists
  if (!dir.exists(clade_path)) {
    pn_warn("Clade directory not found:", clade_path)
    return(data.frame(protein.id = character(), clade = character(), PIGI = character(), stringsAsFactors = FALSE))
  }

  # Get list of clade files
  clade_files <- list.files(clade_path, pattern = pattern, full.names = TRUE)

  # Check if any clade files are found
  if (length(clade_files) == 0) {
    pn_warn("No clade files found matching pattern:", pattern)
    return(data.frame(protein.id = character(), clade = character(), PIGI = character(), stringsAsFactors = FALSE))
  }

  # Log the files found
  pn_info("Found", length(clade_files), "clade files")
  for (file in clade_files) {
    pn_debug("Clade file:", basename(file))
  }

  # Apply the function to each file and clade, then bind rows
  tryCatch({
    clade_assign <- purrr::map2_df(clade_files, LETTERS[seq_along(clade_files)], process_clade_file)
    pn_info("Processed", nrow(clade_assign), "entries from clade files")

    # Add PIGI information
    clade_assign$PIGI <- unlist(purrr::flatten(lapply(clade_assign$protein.id, function(i) { 
      protein.alias(i, alias = PROTEIN_ALIAS, verbose = FALSE)[1, 1] 
    })))

    pn_info("Added PIGI information to clade assignments")

  }, error = function(e) {
    pn_error("Failed to process clade files:", e$message)
    return(data.frame(protein.id = character(), clade = character(), PIGI = character(), stringsAsFactors = FALSE))
  })

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
  pn_debug("Processing clade file:", file, "with clade ID:", clade)

  tryCatch({
    # Read the file and extract protein IDs
    data <- readr::read_delim(file, col_names = c("protein.id", "V2", "V3", "V4"),
                              delim = "/", show_col_types = FALSE)
    data <- data%>%
      dplyr::select(protein.id) %>%
      dplyr::mutate(clade = clade)

    pn_debug("Extracted", nrow(data), "proteins from clade", clade)
    return(data)

  }, error = function(e) {
    pn_error("Failed to process clade file:", file, "-", e$message)
    return(data.frame(protein.id = character(), clade = character(), stringsAsFactors = FALSE))
  })
}

#' Read Representatives
#'
#' This function reads the IPG, PDB, and cluster alias data from text files and combines them.
#'
#' @param PATH The base path to the data directory.
#' @param path The subdirectory containing the alias files.
#' @param ipg_file The name of the IPG alias file.
#' @param pdb_file The name of the PDB alias file.
#' @param cluster_file The name of the cluster alias file.
#' @return A data frame with combined alias data, or an empty data frame if no alias files are found.
#' @export
read_representatives <- function(PATH = 'data',
                                 path = 'representatives', 
                                 ipg_file = 'ipg_representative.txt', 
                                 pdb_file = 'pdb_representative.txt', 
                                 cluster_file = 'cluster_representative.txt') {

  pn_info("Reading protein representatives from:", file.path(PATH, path))
  rep_path <- file.path(PATH, path)

  # Check if representatives directory exists
  if (!dir.exists(rep_path)) {
    pn_warn("Representatives directory not found:", rep_path)
    return(data.frame(alias = character(), PIGI = character(), identity = numeric(), stringsAsFactors = FALSE))
  }

  # Helper function to read alias files
  read_alias_file <- function(file_name, col_names) {
    file_path <- file.path(rep_path, file_name)
    if (file.exists(file_path)) {
      pn_debug("Reading alias file:", file_path)
      tryCatch({
        data <- readr::read_delim(file_path, col_names = col_names, show_col_types = FALSE)
        pn_debug("Read", nrow(data), "entries from", file_name)
        return(data)
      }, error = function(e) {
        pn_error("Failed to read alias file:", file_path, "-", e$message)
        return(NULL)
      })
    } else {
      pn_warn("Alias file not found:", file_path)
      return(NULL)
    }
  }

  # Read the IPG alias data
  ipg_alias <- read_alias_file(ipg_file, c('PIGI', 'alias'))
  if (!is.null(ipg_alias)) {
    ipg_alias$identity <- 1
    pn_info("Read IPG alias data:", nrow(ipg_alias), "entries")
  }
  
  # Read the PDB alias data
  pdb_alias <- read_alias_file(pdb_file, c('alias', 'PIGI'))
  if (!is.null(pdb_alias)) {
    pdb_alias$identity <- 1
    pn_info("Read PDB alias data:", nrow(pdb_alias), "entries")
  }
  
  # Read the cluster alias data
  cluster_alias <- read_alias_file(cluster_file, c('PIGI', 'alias', 'identity'))
  if (!is.null(cluster_alias)) {
    cluster_alias <- cluster_alias %>%
      dplyr::mutate(identity = as.numeric(sub("%", "", identity)) / 100)
    pn_info("Read cluster alias data:", nrow(cluster_alias), "entries")
  }
  
  # Combine the alias data
  combined_alias <- dplyr::bind_rows(ipg_alias, pdb_alias, cluster_alias)
  
  # Check if combined_alias is empty
  if (nrow(combined_alias) == 0) {
    pn_warn("No alias data found. Returning an empty data frame.")
    return(data.frame(alias = character(), PIGI = character(), identity = numeric(), stringsAsFactors = FALSE))
  }
  
  # Replace PIGI values based on PDB alias replacements
  if (!is.null(pdb_alias) && nrow(pdb_alias) > 0) {
    replacements <- setNames(pdb_alias$PIGI, pdb_alias$alias)
    combined_alias$PIGI <- paste0(sub("\\..*", "", stringr::str_replace_all(combined_alias$PIGI, replacements)), ".1")
    pn_info("Applied PDB replacements to PIGI values")
  }
  
  # Remove duplicate entries
  combined_alias <- combined_alias %>%
    dplyr::distinct(alias, PIGI)
  
  pn_info("Final combined representatives data has", nrow(combined_alias), "entries")
  return(combined_alias)
}

#' Create or Validate Output File Path
#'
#' This function checks if the output directory exists and creates it if necessary.
#' It returns the output file path for the neighbors file based on the provided parameters.
#'
#' @param basepairs Number of base pairs used for neighbor detection.
#' @param max_neighbors Maximum number of neighbors considered.
#' @param date Date string for the output directory (default: current date).
#' @param config Configuration list (optional).
#' @param interactive Whether to prompt the user for input (default: TRUE).
#' @return The path to the output file.
#' @export
create_or_validate_output_file_path <- function(basepairs, max_neighbors, date = NULL, 
                                                config = NULL, interactive = TRUE) {
  # Get the current date in YYYY-MM-DD format if not provided
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  
  # If config is provided, use its output directory
  if (!is.null(config)) {
    output_base_dir <- config$paths$output_dir
  } else {
    output_base_dir <- "output"
  }
  
  # If date is provided, use it
  if (!is.null(date)) {
    output_dir <- file.path(output_base_dir, date)
    
    # Check if the directory exists
    if (!dir.exists(output_dir)) {
      if (interactive) {
        choice <- readline(prompt = paste0("Directory ", output_dir, 
                                         " does not exist. Create it? (Y/n): "))
        if (tolower(choice) %in% c("", "y", "yes")) {
          dir.create(output_dir, recursive = TRUE)
          pn_info("Created output directory:", output_dir)
        } else {
          pn_info("Using current date directory instead")
          date <- current_date
          output_dir <- file.path(output_base_dir, date)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
      } else {
        # In non-interactive mode, create the directory
        dir.create(output_dir, recursive = TRUE)
        pn_info("Created output directory:", output_dir)
      }
    }
  } else if (interactive) {
    # If no date provided and interactive, ask user for date
    repeat {
      date <- readline(prompt = "Enter the date of the data you want to load (YYYY-MM-DD) or press Enter for today: ")
      
      if (date == "") {
        date <- current_date
        output_dir <- file.path(output_base_dir, date)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        pn_info("Using current date directory:", output_dir)
        break
      } else {
        # Validate date format
        if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) {
          cat("Invalid date format. Please use YYYY-MM-DD format.\n")
          next
        }
        
        output_dir <- file.path(output_base_dir, date)
        if (dir.exists(output_dir)) {
          pn_info("Using existing directory:", output_dir)
          break
        } else {
          cat("There is no valid folder for your date.\n")
          choice <- readline(prompt = "Do you want to (C)reate this directory, use (T)oday's date, or (E)xit? ")
          if (tolower(choice) == "c") {
            dir.create(output_dir, recursive = TRUE)
            pn_info("Created output directory:", output_dir)
            break
          } else if (tolower(choice) == "t") {
            date <- current_date
            output_dir <- file.path(output_base_dir, date)
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
            pn_info("Using current date directory:", output_dir)
            break
          } else if (tolower(choice) == "e") {
            stop("Exiting the program.")
          }
        }
      }
    }
  } else {
    # Non-interactive mode with no date provided
    date <- current_date
    output_dir <- file.path(output_base_dir, date)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pn_info("Using current date directory:", output_dir)
  }
  
  # Construct the output file path
  output_file <- file.path(output_dir, paste0('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv'))
  pn_info("Output file path:", output_file)
  
  return(output_file)
}

#' Read Neighbour Annotations
#'
#' This function re-imports the types of neighbours after manual annotation.
#'
#' @param current_date A string representing the current date, used to locate the annotated CSV file.
#' @param output_dir Alternative output directory (overrides date-based directory).
#' @param empty_cells A string representing the value to replace empty cells with. Default is "unknown".
#' @return A data frame containing the annotated neighbours if the file exists, otherwise NULL.
#' @export
read_annotations <- function(current_date, output_dir = NULL, empty_cells = "unknown") {
  # Determine the file path
  if (is.null(output_dir)) {
    file_path <- file.path("output", current_date, "types_of_neighbours_annotated.csv")
  } else {
    file_path <- file.path(output_dir, "types_of_neighbours_annotated.csv")
  }
  
  pn_info("Checking for annotated neighbours file:", file_path)
  
  # Re-import the types of neighbours after manual annotation
  if(file.exists(file_path)) {
    pn_info("Found annotated neighbours file")
    
    tryCatch({
      annotated_neighbours <- readr::read_csv(file_path, show_col_types = FALSE)
      
      # Replace empty cells with NA
      annotated_neighbours <- annotated_neighbours %>%
        dplyr::mutate_all(~ifelse(. == "", NA, .))
      
      # Replace NA with the specified value
      if (!is.null(empty_cells)) {
        annotated_neighbours <- annotated_neighbours %>%
          tidyr::replace_na(list(COG_NAME = empty_cells, COG_LETTER = empty_cells, 
                            N = empty_cells, ANNOTATION = empty_cells))
      }
      
      # Rename columns
      colnames(annotated_neighbours)[1:4] <- c("COG_NAME", "COG_LETTER", "N", "ANNOTATION")
      
      pn_info("Read", nrow(annotated_neighbours), "annotated neighbour entries")
      return(annotated_neighbours)
      
    }, error = function(e) {
      pn_error("Failed to read annotated neighbours file:", e$message)
      return(NULL)
    })
  } else {
    pn_info("No annotated neighbours file found")
    return(NULL)
  }
}
