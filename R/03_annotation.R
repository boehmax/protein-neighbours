#' Annotation Module
#'
#' Functions for annotating proteins using eggNOG-mapper or alternative methods.
#'
#' @author Protein Neighbors Package

#' Clean up eggNOG output files for a fresh run
#'
#' This function removes all eggNOG-mapper output files to allow for a clean run.
#'
#' @param eggnog_dir Directory with eggNOG-mapper outputs
#' @param confirm Whether to confirm before deleting (default: TRUE)
#' @return TRUE if successful, FALSE otherwise
#' @export
cleanup_eggnog_results <- function(eggnog_dir, confirm = TRUE) {
  # Check if the directory exists
  if (!dir.exists(eggnog_dir)) {
    pn_info("No eggNOG results directory found at:", eggnog_dir)
    return(TRUE)
  }
  
  # List all files that match eggNOG output patterns
  eggnog_files <- list.files(
    eggnog_dir, 
    pattern = "eggnog_results\\.emapper\\..*", 
    full.names = TRUE
  )
  
  if (length(eggnog_files) == 0) {
    pn_info("No eggNOG output files found in:", eggnog_dir)
    return(TRUE)
  }
  
  # Check if we need confirmation
  if (confirm) {
    message("Found ", length(eggnog_files), " eggNOG output files in: ", eggnog_dir)
    message("Files to be deleted:")
    for (file in eggnog_files) {
      message("  ", basename(file))
    }
    
    answer <- readline(prompt = "Delete these files? (y/n): ")
    if (tolower(answer) != "y") {
      message("Cleanup cancelled.")
      return(FALSE)
    }
  }
  
  # Delete the files
  pn_info("Removing", length(eggnog_files), "eggNOG output files")
  deleted <- file.remove(eggnog_files)
  
  # Check if all files were deleted
  if (all(deleted)) {
    pn_info("Successfully removed all eggNOG output files")
    return(TRUE)
  } else {
    failed <- eggnog_files[!deleted]
    pn_warn("Failed to remove some eggNOG output files:")
    for (file in failed) {
      pn_warn("  ", basename(file))
    }
    return(FALSE)
  }
}

#' Analyze proteins using eggNOG or alternative methods
#'
#' This function annotates proteins using either eggNOG-mapper, COG classifier, 
#' or generates dummy annotations for testing.
#'
#' @param df A data frame with protein information.
#' @param column The column name containing protein IDs.
#' @param config The configuration list.
#' @param force Logical, whether to force recomputation of results (default: FALSE)
#' @return A data frame with the annotation results.
#' @export
analyze_proteins <- function(df, column = 'ID', config, force = FALSE) {
  # Create output directories
  current_date <- config$analysis$date
  output_dir <- file.path(config$paths$output_dir, current_date)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Log the start of protein annotation
  pn_info("Starting protein annotation with", config$annotation$tool)
  pn_info(paste("Annotating", nrow(df), "proteins"))
  
  # Export accessions of the proteins
  protein_accessions <- data.frame(df %>% dplyr::select(all_of(column)))
  accessions_file <- file.path(output_dir, "protein_accessions.csv")
  readr::write_csv(protein_accessions, accessions_file, col_names = FALSE)
  
  # Get the FASTA file of the proteins
  fasta_file <- file.path(output_dir, "proteins.fasta")
  
  # Check if the FASTA file already exists and has content
  fasta_exists <- FALSE
  if (!force && file.exists(fasta_file) && file.size(fasta_file) > 0) {
    # Check if it contains actual FASTA sequences (lines starting with >)
    tryCatch({
      fasta_content <- readLines(fasta_file, n = 100)  # Read first 100 lines
      if (any(grepl("^>", fasta_content))) {
        # Count sequences for verification
        seq_count <- sum(grepl("^>", fasta_content))
        
        # If we found at least a few sequences, consider it valid
        if (seq_count >= 1) {
          pn_info("Found existing FASTA file with", seq_count, "sequences (first 100 lines)")
          
          # Read the whole file to get a more accurate count if there are many sequences
          if (seq_count >= 50) {  # If we found a lot in the first 100 lines, there may be more
            all_content <- readLines(fasta_file)
            total_seq_count <- sum(grepl("^>", all_content))
            pn_info("Total sequences in existing FASTA file:", total_seq_count)
          }
          
          fasta_exists <- TRUE
          pn_info("Using existing FASTA file, skipping download")
        }
      }
    }, error = function(e) {
      pn_warn("Error checking existing FASTA file:", e$message)
    })
  }
  
  if (!fasta_exists) {
    pn_info("Retrieving FASTA sequences for proteins")
    
    # Use rentrez to get sequences
    if (requireNamespace("rentrez", quietly = TRUE)) {
      pn_info("Using rentrez to retrieve sequences")
      success <- retrieve_fasta_sequences(accessions_file, fasta_file)
      if (!success) {
        pn_error("Failed to retrieve FASTA sequences using rentrez")
        return(NULL)
      }
    } else {
      # Try to use efetch directly
      pn_info("rentrez not available, trying to use efetch directly")
      cmd <- paste("efetch -db protein -id $(cat", shQuote(accessions_file), 
                   "| tr '\\n' ',') -format fasta >", shQuote(fasta_file))
      status <- system(cmd)
      if (status != 0 || !file.exists(fasta_file) || file.size(fasta_file) == 0) {
        pn_error("Failed to retrieve FASTA sequences using efetch")
        
        # Last resort: generate dummy sequences for testing
        pn_warn("Generating dummy sequences for testing purposes")
        generate_dummy_sequences(protein_accessions[[1]], fasta_file)
      }
    }
  }
  
  # Check if FASTA file exists and has content (after potential download)
  if (!file.exists(fasta_file) || file.size(fasta_file) == 0) {
    pn_error("No valid FASTA sequences available for annotation")
    return(NULL)
  }
  
  # Run the annotation tool
  if (config$annotation$tool == "eggnog") {
    # Set up eggNOG output directory
    eggnog_dir <- file.path(output_dir, "eggnog")
    dir.create(eggnog_dir, recursive = TRUE, showWarnings = FALSE)
    
    # If force is TRUE, clean up existing results
    if (force) {
      pn_info("Force flag is set - cleaning up existing results")
      cleanup_eggnog_results(eggnog_dir, confirm = FALSE)
    }
    
    # Check for existing eggNOG output files
    annotations_file <- file.path(eggnog_dir, "eggnog_results.emapper.annotations")
    classifier_file <- file.path(eggnog_dir, "classifier_result.tsv")
    
    # Check if we already have processed results
    if (!force && file.exists(classifier_file) && file.size(classifier_file) > 0) {
      # We already have processed results, just load them
      pn_info("Found existing eggNOG annotation results, loading from:", classifier_file)
      try({
        cog_classification <- read.delim(classifier_file)
        pn_info("Loaded", nrow(cog_classification), "existing annotations")
        return(cog_classification)
      }, silent = TRUE)
      # If loading fails, continue with reprocessing
    }
    
    # Check if we have raw eggNOG annotations that need processing
    if (!force && file.exists(annotations_file) && file.size(annotations_file) > 0) {
      # We have eggNOG results but need to process them
      pn_info("Found existing eggNOG annotations, processing:", annotations_file)
      annotation_result <- process_eggnog_output(annotations_file, eggnog_dir)
      if (!is.null(annotation_result)) {
        pn_info("Successfully processed existing eggNOG annotations")
        return(annotation_result)
      }
      # If processing fails, continue with rerunning eggNOG
    }
    
    # First check if eggNOG-mapper is available
    emapper_available <- system("which emapper.py >/dev/null 2>&1") == 0
    
    if (emapper_available) {
      pn_info("eggNOG-mapper found. Using it for annotation.")
      
      # Extract eggNOG parameters from config
      eggnog_config <- config$annotation$eggnog
      
      # Make sure the database directory exists and is absolute path
      db_dir <- normalizePath(eggnog_config$db_dir, mustWork = FALSE)
      
      # Check if database directory exists and has files
      if (dir.exists(db_dir) && length(list.files(db_dir)) > 0) {
        # Create the tmp directory
        tmp_dir <- file.path(eggnog_dir, "tmp")
        dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
        file_to_resume  <- normalizePath(file.path(eggnog_dir, "eggnog_results.emapper.hits"), mustWork = FALSE)
        
        # Build the command with resume flag
        cmd <- paste(
          "cd", shQuote(eggnog_dir), "&&",
          "emapper.py",
          "-i", shQuote(normalizePath(fasta_file, mustWork = FALSE)),
          "--output eggnog_results",
          "--data_dir", shQuote(normalizePath(db_dir, mustWork = FALSE)),
          "--cpu", eggnog_config$cpu,
          "--temp_dir", shQuote(normalizePath(tmp_dir, mustWork = FALSE)),
          "--tax_scope", eggnog_config$tax_scope,
          "--go_evidence non-electronic",
          "--target_orthologs all",
          "--seed_ortholog_evalue", eggnog_config$seed_ortholog_evalue,
          "--seed_ortholog_score", eggnog_config$seed_ortholog_score,
          # Add resume flag to continue from existing results if available or override to start fresh
          ifelse(force, "--override", ifelse(file.exists(file_to_resume) ,"--resume",""))
        )
        
        # Run eggNOG-mapper
        pn_info("Running eggNOG-mapper command:", cmd)
        status <- system(cmd)
        
        if (status == 0) {
          # Check if the output file exists
          if (file.exists(annotations_file)) {
            pn_info("eggNOG-mapper completed successfully")
            # Process the results
            annotation_result <- process_eggnog_output(annotations_file, eggnog_dir)
            if (!is.null(annotation_result)) {
              return(annotation_result)
            }
          }
        }
        
        pn_error("eggNOG-mapper failed with status:", status)
      } else {
        pn_warn("eggNOG database directory doesn't exist or is empty:", db_dir)
      }
    } else {
      pn_warn("eggNOG-mapper not found in PATH.")
    }
    
    # Fallback to dummy annotations
    pn_warn("Falling back to dummy annotations")
    dummy_dir <- file.path(output_dir, "dummy_annotations")
    annotation_result <- create_dummy_annotations(unique(df[[column]]), dummy_dir)
    return(annotation_result)
    
  } else if (config$annotation$tool == "cog") {
    # Run COG classifier
    annotation_result <- run_cog_classifier(fasta_file, output_dir, config)
    return(annotation_result)
    
  } else if (config$annotation$tool == "dummy") {
    # Use dummy annotations
    pn_info("Using dummy annotation as specified in config")
    dummy_dir <- file.path(output_dir, "dummy_annotations")
    annotation_result <- create_dummy_annotations(unique(df[[column]]), dummy_dir)
    return(annotation_result)
    
  } else {
    pn_error("Unknown annotation tool:", config$annotation$tool)
    return(NULL)
  }
}

#' Retrieve FASTA sequences from NCBI using R
#'
#' @param accession_file Path to the file containing accession numbers
#' @param output_file Path where to save the FASTA sequences
#' @param batch_size Number of accessions to retrieve at once (default: 50)
#' @param sleep_time Time to sleep between batches in seconds (default: 0.5)
#' @return TRUE if successful, FALSE otherwise
#' @keywords internal
retrieve_fasta_sequences <- function(accession_file, output_file, batch_size = 50, sleep_time = 0.5) {
  # Check if the rentrez package is available
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    pn_error("Package 'rentrez' is needed for retrieving sequences. Please install it.")
    return(FALSE)
  }
  
  # Read accession numbers
  tryCatch({
    accessions <- readLines(accession_file)
    accessions <- accessions[accessions != ""]  # Remove empty lines
  }, error = function(e) {
    pn_error("Failed to read accession file:", e$message)
    return(FALSE)
  })
  
  pn_info("Retrieved", length(accessions), "accession numbers")
  
  # Create or clear the output file
  file.create(output_file)
  
  # Process in batches
  total_batches <- ceiling(length(accessions) / batch_size)
  successful_count <- 0
  
  for (batch_num in 1:total_batches) {
    # Get accessions for this batch
    start_idx <- (batch_num - 1) * batch_size + 1
    end_idx <- min(batch_num * batch_size, length(accessions))
    batch_accessions <- accessions[start_idx:end_idx]
    
    # Print progress
    pn_info(sprintf("Processing batch %d of %d (accessions %d-%d)", 
                    batch_num, total_batches, start_idx, end_idx))
    
    # Try to fetch sequences with retries
    max_retries <- 3
    success <- FALSE
    
    for (attempt in 1:max_retries) {
      tryCatch({
        # Fetch sequences
        sequences <- rentrez::entrez_fetch(
          db = "protein",
          id = batch_accessions,
          rettype = "fasta",
          retmode = "text"
        )
        
        # Append to output file
        cat(sequences, file = output_file, append = TRUE)
        
        # Count sequences
        seq_count <- length(gregexpr(">", sequences)[[1]])
        if (seq_count > 0) {
          successful_count <- successful_count + seq_count
          pn_info(sprintf("Retrieved %d sequences in this batch", seq_count))
          success <- TRUE
          break
        } else {
          pn_warn("No sequences retrieved in this batch, retrying...")
        }
      }, error = function(e) {
        pn_warn(sprintf("Attempt %d failed: %s", attempt, e$message))
        # Wait longer between retries
        Sys.sleep(sleep_time * attempt)
      })
    }
    
    if (!success) {
      pn_error(sprintf("Failed to retrieve sequences for batch %d after %d attempts", 
                       batch_num, max_retries))
    }
    
    # Sleep between batches to avoid overwhelming the NCBI server
    Sys.sleep(sleep_time)
  }
  
  # Check if we got any sequences
  if (successful_count == 0) {
    pn_error("Failed to retrieve any sequences")
    return(FALSE)
  }
  
  pn_info(sprintf("Successfully retrieved %d sequences out of %d accessions", 
                  successful_count, length(accessions)))
  return(TRUE)
}

#' Generate dummy sequences for testing
#'
#' @param protein_ids Vector of protein IDs
#' @param output_file Path to save the FASTA file
#' @return TRUE if successful
#' @keywords internal
generate_dummy_sequences <- function(protein_ids, output_file) {
  # Create dummy sequences
  fasta_content <- character()
  
  for (id in protein_ids) {
    # Generate a random sequence of 300-500 amino acids
    seq_length <- sample(300:500, 1)
    amino_acids <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
                     "M", "F", "P", "S", "T", "W", "Y", "V")
    sequence <- paste(sample(amino_acids, seq_length, replace = TRUE), collapse = "")
    
    # Add to FASTA content
    fasta_content <- c(fasta_content, 
                       paste0(">", id),
                       sequence)
  }
  
  # Write to file
  writeLines(fasta_content, output_file)
  
  pn_warn("Generated dummy sequences for testing purposes. DO NOT USE FOR REAL ANALYSIS.")
  return(TRUE)
}

# Replace the process_eggnog_output function in R/03_annotation.R

#' Process eggNOG-mapper output into COG classifier format
#'
#' @param annotations_file Path to the eggNOG-mapper annotations file
#' @param output_dir Directory to save processed files
#' @return A data frame with COG classifications or NULL if failed
#' @keywords internal
process_eggnog_output <- function(annotations_file, output_dir) {
  # Create a file with COG category information
  cog_categories_file <- file.path(output_dir, "cog_categories.tsv")
  cog_categories <- '
COG_LETTER	COG_CATEGORY	DESCRIPTION
J	INFORMATION STORAGE AND PROCESSING	Translation, ribosomal structure and biogenesis
A	INFORMATION STORAGE AND PROCESSING	RNA processing and modification
K	INFORMATION STORAGE AND PROCESSING	Transcription
L	INFORMATION STORAGE AND PROCESSING	Replication, recombination and repair
B	INFORMATION STORAGE AND PROCESSING	Chromatin structure and dynamics
D	CELLULAR PROCESSES AND SIGNALING	Cell cycle control, cell division, chromosome partitioning
Y	CELLULAR PROCESSES AND SIGNALING	Nuclear structure
V	CELLULAR PROCESSES AND SIGNALING	Defense mechanisms
T	CELLULAR PROCESSES AND SIGNALING	Signal transduction mechanisms
M	CELLULAR PROCESSES AND SIGNALING	Cell wall/membrane/envelope biogenesis
N	CELLULAR PROCESSES AND SIGNALING	Cell motility
Z	CELLULAR PROCESSES AND SIGNALING	Cytoskeleton
W	CELLULAR PROCESSES AND SIGNALING	Extracellular structures
U	CELLULAR PROCESSES AND SIGNALING	Intracellular trafficking, secretion, and vesicular transport
O	CELLULAR PROCESSES AND SIGNALING	Posttranslational modification, protein turnover, chaperones
X	MOBILOME	Mobilome: prophages, transposons
C	METABOLISM	Energy production and conversion
G	METABOLISM	Carbohydrate transport and metabolism
E	METABOLISM	Amino acid transport and metabolism
F	METABOLISM	Nucleotide transport and metabolism
H	METABOLISM	Coenzyme transport and metabolism
I	METABOLISM	Lipid transport and metabolism
P	METABOLISM	Inorganic ion transport and metabolism
Q	METABOLISM	Secondary metabolites biosynthesis, transport and metabolism
R	POORLY CHARACTERIZED	General function prediction only
S	POORLY CHARACTERIZED	Function unknown
'
  # Write COG categories to file
  writeLines(trimws(cog_categories), cog_categories_file)
  
  # First, examine the raw eggNOG output file directly to debug
  pn_info("Examining raw eggNOG annotations file:", annotations_file)
  if (!file.exists(annotations_file)) {
    pn_error("eggNOG annotations file not found:", annotations_file)
    return(NULL)
  }
  
  # Debug: Check file size and content
  file_size <- file.size(annotations_file)
  pn_info("eggNOG annotations file size:", file_size, "bytes")
  
  if (file_size == 0) {
    pn_error("eggNOG annotations file is empty")
    return(NULL)
  }
  
  # Debug: Check first few lines
  first_lines <- readLines(annotations_file, n = 20)
  pn_info("First few lines of eggNOG annotations file:")
  for (i in 1:min(10, length(first_lines))) {
    pn_info(first_lines[i])
  }
  
  # Count comment lines and find header line
  comment_count <- sum(grepl("^##", first_lines))
  header_line <- which(grepl("^#[^#]", first_lines))[1]
  
  pn_info("Found", comment_count, "comment lines and header at line", header_line)
  
  # Try to read the annotations file with proper handling of comment lines
  tryCatch({
    if (is.na(header_line)) {
      pn_error("Could not find header line in eggNOG output")
      return(NULL)
    }
    
    # Get the header line and remove the leading #
    header <- sub("^#", "", first_lines[header_line])
    header_fields <- strsplit(header, "\t")[[1]]
    
    pn_info("Header fields:", paste(header_fields, collapse=", "))
    
    # Read the file with specified column names
    annotations <- read.delim(
      annotations_file, 
      skip = header_line,
      header = FALSE,
      col.names = header_fields,
      stringsAsFactors = FALSE,
      comment.char = "#",
      quote = "",
      na.strings = c("-", "")
    )
    
    # Debug info
    pn_info("Successfully read", nrow(annotations), "rows from eggNOG annotations")
    pn_info("Columns found:", paste(colnames(annotations), collapse=", "))
    
    # Check if we have COG categories
    has_cog_categories <- FALSE
    if ("COG_category" %in% colnames(annotations)) {
      pn_info("Found COG_category column")
      has_cog_categories <- TRUE
    } else {
      pn_warn("No COG_category column found in eggNOG annotations")
      # Try alternative column names that might contain COG categories
      possible_cols <- c("COG_categories", "cog_category", "cog_categories")
      for (col in possible_cols) {
        if (col %in% colnames(annotations)) {
          pn_info(paste("Using", col, "column instead of COG_category"))
          colnames(annotations)[colnames(annotations) == col] <- "COG_category"
          has_cog_categories <- TRUE
          break
        }
      }
    }
    
    # Check if we have Description
    has_description <- FALSE
    if ("Description" %in% colnames(annotations)) {
      pn_info("Found Description column")
      has_description <- TRUE
    } else {
      pn_warn("No Description column found in eggNOG annotations")
      # Try alternative column names that might contain descriptions
      possible_cols <- c("description", "desc", "protein_description", "product")
      for (col in possible_cols) {
        if (col %in% colnames(annotations)) {
          pn_info(paste("Using", col, "column instead of Description"))
          colnames(annotations)[colnames(annotations) == col] <- "Description"
          has_description <- TRUE
          break
        }
      }
    }
#    
#    # Initialize classifier results
#    classifier_results <- data.frame(
#      QUERY_ID = character(),
#      COG_ID = character(),
#      CDD_ID = character(),
#      EVALUE = character(),
#      GENE_NAME = character(),
#      COG_NAME = character(),
#      COG_LETTER = character(),
#      COG_DESCRIPTION = character(),
#      stringsAsFactors = FALSE
#    )
#    
#    # Read COG categories for reference
#    cog_cats <- read.delim(cog_categories_file, stringsAsFactors = FALSE)
#    
#    # Process each annotation
#    for (i in 1:nrow(annotations)) {
#      # Get the query ID
#      query_id <- annotations$query[i]
#      
#      # Get COG letters from COG_category column
#      cog_letters <- c()
#      if (has_cog_categories && !is.na(annotations$COG_category[i]) && annotations$COG_category[i] != "-") {
#        cog_letters <- strsplit(annotations$COG_category[i], "")[[1]]
#      }
#      
#      # If no COG category was found, default to S (unknown)
#      if (length(cog_letters) == 0) {
#        cog_letters <- "S"
#      }
#      
#      # Get COG ID from eggNOG_OGs
#      cog_id <- NULL
#      if ("eggNOG_OGs" %in% colnames(annotations) && 
#          !is.na(annotations$eggNOG_OGs[i]) && 
#          annotations$eggNOG_OGs[i] != "-") {
#        ogs <- strsplit(annotations$eggNOG_OGs[i], ",")[[1]]
#        cog_matches <- grep("^COG[0-9]+@", ogs, value = TRUE)
#        if (length(cog_matches) > 0) {
#          cog_id <- sub("@.*$", "", cog_matches[1])
#        }
#      }
#      
#      # If no COG ID was found place as NA
#      if (is.null(cog_id)) {
#        cog_id <- NA
#      }
#      
#      # Get description
#      description <- "Function unknown"
#      if ("Description" %in% colnames(annotations) && 
#          !is.na(annotations$Description[i]) && 
#          annotations$Description[i] != "-") {
#        description <- annotations$Description[i]
#      }
#      
#      # Get e-value
#      evalue <- "-"
#      if ("evalue" %in% colnames(annotations) && !is.na(annotations$evalue[i])) {
#        evalue <- annotations$evalue[i]
#      }
#      
#      # Get preferred name
#      gene_name <- "-"
#      if ("Preferred_name" %in% colnames(annotations) && 
#          !is.na(annotations$Preferred_name[i]) && 
#          annotations$Preferred_name[i] != "-") {
#        gene_name <- annotations$Preferred_name[i]
#      }
#      
#      # Create an entry for each COG letter
#      for (letter in cog_letters) {
#        # Find category description
#        cat_idx <- which(cog_cats$COG_LETTER == letter)
#        if (length(cat_idx) > 0) {
#          category <- cog_cats$COG_CATEGORY[cat_idx]
#          
#          # Add to results
#          classifier_results <- rbind(classifier_results, data.frame(
#            QUERY_ID = query_id,
#            COG_ID = cog_id,
#            CDD_ID = "-",
#            EVALUE = evalue,
#            GENE_NAME = gene_name,
#            COG_NAME = category,
#            COG_LETTER = letter,
#            COG_DESCRIPTION = description,
#            stringsAsFactors = FALSE
#          ))
#        }
#      }
#    }
#    
#    # If no annotations were found, create at least one dummy entry
#    if (nrow(classifier_results) == 0) {
#      pn_warn("No COG annotations could be extracted. Creating a diverse set of dummy annotations.")
#      
#      # Create diverse dummy annotations instead of all S
#      protein_ids <- unique(annotations$query)
#      dummy_results <- create_diverse_dummy_annotations(protein_ids)
#      classifier_results <- dummy_results
#    } else {
#      pn_info("Successfully extracted", nrow(classifier_results), "COG annotations")
#      
#      # Check for diversity in COG letters
#      cog_letter_counts <- table(classifier_results$COG_LETTER)
#      pn_info("COG letter distribution in processed results:")
#      for (letter in names(cog_letter_counts)) {
#        pn_info("  ", letter, ":", cog_letter_counts[letter])
#      }
#      
#      # If only S is present, create more diverse annotations
#      if (length(cog_letter_counts) == 1 && names(cog_letter_counts)[1] == "S") {
#        pn_warn("Only 'S' (Function unknown) annotations found. Creating more diverse annotations.")
#        protein_ids <- unique(annotations$query)
#        classifier_results <- create_diverse_dummy_annotations(protein_ids)
#      }
#    }
#    
#    # Write results to file
#    classifier_file <- file.path(output_dir, "classifier_result.tsv")
#    write.table(classifier_results, classifier_file, sep = "\t", row.names = FALSE, quote = FALSE)
#    
#    pn_info("Successfully processed eggNOG annotations into", nrow(classifier_results), "COG classifier entries")
#    
#  
#    
#    return(classifier_results)
#    
#  }, error = function(e) {
#    pn_error("Failed to process eggNOG-mapper output:", e$message)
#    pn_error("Error traceback:", conditionCall(e))
#    
#    # Create more diverse dummy annotations as fallback
#    pn_warn("Creating diverse dummy annotations as fallback")
#    
#    # Try to extract protein IDs from the annotations file
#    protein_ids <- tryCatch({
#      # Get protein IDs from header lines
#      lines <- readLines(annotations_file)
#      data_lines <- lines[!grepl("^#", lines)]
#      
#      if (length(data_lines) > 0) {
#        # Extract first column (query ID)
#        unique(sapply(strsplit(data_lines, "\t"), function(x) x[1]))
#      } else {
#        c("dummy1", "dummy2", "dummy3")  # Fallback
#      }
#    }, error = function(e2) {
#      c("dummy1", "dummy2", "dummy3")  # Fallback
#    })
#    
#    
    # Write to file
    classifier_results <- annotations
    classifier_file <- file.path(output_dir, "classifier_result.tsv")
    write.table(classifier_results, classifier_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    return(classifier_results)
  })
}

#' Create diverse dummy COG annotations for testing
#'
#' @param protein_ids Vector of protein IDs to annotate
#' @return A data frame with diverse dummy COG annotations
#' @keywords internal
create_diverse_dummy_annotations <- function(protein_ids) {
  # Define COG categories and their descriptions
  cog_categories <- data.frame(
    COG_LETTER = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", 
                   "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
    COG_NAME = c(rep("INFORMATION STORAGE AND PROCESSING", 5),
                 rep("CELLULAR PROCESSES AND SIGNALING", 10),
                 "MOBILOME",
                 rep("METABOLISM", 8),
                 rep("POORLY CHARACTERIZED", 2)),
    COG_DESCRIPTION = c(
      "Translation, ribosomal structure and biogenesis",
      "RNA processing and modification",
      "Transcription",
      "Replication, recombination and repair",
      "Chromatin structure and dynamics",
      "Cell cycle control, cell division, chromosome partitioning",
      "Nuclear structure",
      "Defense mechanisms",
      "Signal transduction mechanisms",
      "Cell wall/membrane/envelope biogenesis",
      "Cell motility",
      "Cytoskeleton",
      "Extracellular structures",
      "Intracellular trafficking, secretion, and vesicular transport",
      "Posttranslational modification, protein turnover, chaperones",
      "Mobilome: prophages, transposons",
      "Energy production and conversion",
      "Carbohydrate transport and metabolism",
      "Amino acid transport and metabolism",
      "Nucleotide transport and metabolism",
      "Coenzyme transport and metabolism",
      "Lipid transport and metabolism",
      "Inorganic ion transport and metabolism",
      "Secondary metabolites biosynthesis, transport and metabolism",
      "General function prediction only",
      "Function unknown"
    ),
    stringsAsFactors = FALSE
  )
  
  # Ensure a wide range of COG categories, not just S
  # To make it more realistic, weight the probability of common categories
  common_categories <- c("K", "E", "T", "M", "C", "G", "P", "O", "L", "J", "S", "R")
  weights <- rep(1, nrow(cog_categories))
  weights[match(common_categories, cog_categories$COG_LETTER)] <- 3  # 3x more likely
  
  # Create dummy results with diverse COG assignments
  results <- data.frame()
  
  for (id in protein_ids) {
    # Assign 1-3 COG categories to each protein
    num_cogs <- sample(1:3, 1)
    
    # For each protein, select categories without replacement
    selected_indices <- sample(1:nrow(cog_categories), num_cogs, prob = weights, replace = FALSE)
    
    for (idx in selected_indices) {
      # Create a random COG ID
      cog_id <- paste0("COG", sprintf("%04d", sample(1:9999, 1)))
      
      # Add to results
      results <- rbind(results, data.frame(
        QUERY_ID = id,
        COG_ID = cog_id,
        CDD_ID = "-",
        EVALUE = paste0(runif(1, min = 1e-100, max = 1e-10)),
        GENE_NAME = "-",
        COG_NAME = cog_categories$COG_NAME[idx],
        COG_LETTER = cog_categories$COG_LETTER[idx],
        COG_DESCRIPTION = cog_categories$COG_DESCRIPTION[idx],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  pn_warn("Created", nrow(results), "diverse dummy annotations for", length(protein_ids), "proteins")
  return(results)
}

#' Create dummy COG annotations for testing
#'
#' @param protein_ids Vector of protein IDs to annotate
#' @param output_dir Directory to save the dummy results
#' @return A data frame with dummy COG annotations
#' @export
create_dummy_annotations <- function(protein_ids, output_dir) {
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Define COG categories and their descriptions
  cog_categories <- data.frame(
    COG_LETTER = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", 
                   "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
    COG_NAME = c(rep("INFORMATION STORAGE AND PROCESSING", 5),
                 rep("CELLULAR PROCESSES AND SIGNALING", 10),
                 "MOBILOME",
                 rep("METABOLISM", 8),
                 rep("POORLY CHARACTERIZED", 2)),
    COG_DESCRIPTION = c(
      "Translation, ribosomal structure and biogenesis",
      "RNA processing and modification",
      "Transcription",
      "Replication, recombination and repair",
      "Chromatin structure and dynamics",
      "Cell cycle control, cell division, chromosome partitioning",
      "Nuclear structure",
      "Defense mechanisms",
      "Signal transduction mechanisms",
      "Cell wall/membrane/envelope biogenesis",
      "Cell motility",
      "Cytoskeleton",
      "Extracellular structures",
      "Intracellular trafficking, secretion, and vesicular transport",
      "Posttranslational modification, protein turnover, chaperones",
      "Mobilome: prophages, transposons",
      "Energy production and conversion",
      "Carbohydrate transport and metabolism",
      "Amino acid transport and metabolism",
      "Nucleotide transport and metabolism",
      "Coenzyme transport and metabolism",
      "Lipid transport and metabolism",
      "Inorganic ion transport and metabolism",
      "Secondary metabolites biosynthesis, transport and metabolism",
      "General function prediction only",
      "Function unknown"
    ),
    stringsAsFactors = FALSE
  )
  
  # Create dummy results with random COG assignments
  results <- data.frame()
  
  for (id in protein_ids) {
    # Randomly assign 1-3 COG categories to each protein
    num_cogs <- sample(1:3, 1)
    
    for (i in 1:num_cogs) {
      # Randomly select a COG category
      cog_idx <- sample(1:nrow(cog_categories), 1)
      
      # Create a random COG ID
      cog_id <- paste0("COG", sprintf("%04d", sample(1:9999, 1)))
      
      # Add to results
      results <- rbind(results, data.frame(
        QUERY_ID = id,
        COG_ID = cog_id,
        CDD_ID = "-",
        EVALUE = paste0(runif(1, min = 1e-100, max = 1e-10)),
        GENE_NAME = "-",
        COG_NAME = cog_categories$COG_NAME[cog_idx],
        COG_LETTER = cog_categories$COG_LETTER[cog_idx],
        COG_DESCRIPTION = cog_categories$COG_DESCRIPTION[cog_idx],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save to file
  output_file <- file.path(output_dir, "classifier_result.tsv")
  write.table(results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  pn_warn("Using DUMMY annotations. This is only for testing purposes!")
  pn_warn("Install eggNOG-mapper for real functional annotations.")
  
  return(results)
}

#' Run COG classifier for protein annotation
#'
#' This function runs COG classifier on a FASTA file to annotate proteins.
#'
#' @param fasta_file Path to the FASTA file.
#' @param output_dir Output directory for results.
#' @param config The configuration list.
#' @return A data frame with the annotation results.
#' @keywords internal
run_cog_classifier <- function(fasta_file, output_dir, config) {
  pn_info("Running COG classifier for protein annotation")
  
  # Create COG output directory
  cog_dir <- file.path(output_dir, "cogclassifier")
  dir.create(cog_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if COGclassifier is available
  cog_cmd <- Sys.which("COGclassifier")
  if (cog_cmd == "") {
    pn_error("COGclassifier not found in PATH. Falling back to dummy annotations.")
    
    # Read FASTA file to get protein IDs
    fasta_content <- readLines(fasta_file)
    header_lines <- grep("^>", fasta_content)
    protein_ids <- gsub("^>", "", fasta_content[header_lines])
    
    # Create dummy annotations
    return(create_dummy_annotations(protein_ids, cog_dir))
  }
  
  # Run COGclassifier
  cog_command <- paste("COGclassifier -i", shQuote(fasta_file), "-o", shQuote(cog_dir))
  pn_info("Running command:", cog_command)
  
  system_result <- system(cog_command)
  if (system_result != 0) {
    pn_error("Failed to run COG classifier")
    
    # Read FASTA file to get protein IDs
    fasta_content <- readLines(fasta_file)
    header_lines <- grep("^>", fasta_content)
    protein_ids <- gsub("^>", "", fasta_content[header_lines])
    
    # Create dummy annotations
    return(create_dummy_annotations(protein_ids, cog_dir))
  }
  
  # Read the results
  result_file <- file.path(cog_dir, "classifier_result.tsv")
  if (!file.exists(result_file)) {
    pn_error("COG classifier result file not found:", result_file)
    
    # Read FASTA file to get protein IDs
    fasta_content <- readLines(fasta_file)
    header_lines <- grep("^>", fasta_content)
    protein_ids <- gsub("^>", "", fasta_content[header_lines])
    
    # Create dummy annotations
    return(create_dummy_annotations(protein_ids, cog_dir))
  }
  
  cog_classification <- read.delim(result_file)
  
  # Log success
  pn_info("Successfully annotated proteins with COG classifier")
  pn_info(paste("Found", nrow(cog_classification), "annotations"))
  
  return(cog_classification)
}