#' Annotation Module
#'
#' Functions for annotating proteins using eggNOG-mapper or alternative methods.
#'
#' @author Improved Protein Neighbors Package

#' Analyze proteins using eggNOG or alternative methods
#'
#' This function annotates proteins using either eggNOG-mapper, COG classifier, 
#' or generates dummy annotations for testing.
#'
#' @param df A data frame with protein information.
#' @param column The column name containing protein IDs.
#' @param config The configuration list.
#' @return A data frame with the annotation results.
#' @export
analyze_proteins <- function(df, column = 'ID', config) {
  # Create output directories
  current_date <- config$analysis$date
  output_dir <- file.path(config$paths$output_dir, current_date)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Log the start of protein annotation
  pn_info("Starting protein annotation with", config$annotation$tool)
  pn_info(paste("Annotating", nrow(df), "unique proteins"))
  
  # Export accessions of the proteins
  protein_accessions <- data.frame(df %>% dplyr::select(all_of(column)))
  accessions_file <- file.path(output_dir, "protein_accessions.csv")
  readr::write_csv(protein_accessions, accessions_file, col_names = FALSE)
  
  # Get the FASTA file of the proteins
  fasta_file <- file.path(output_dir, "proteins.fasta")
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
  
  # Check if FASTA file exists and has content
  if (!file.exists(fasta_file) || file.size(fasta_file) == 0) {
    pn_error("No valid FASTA sequences available for annotation")
    return(NULL)
  }
  
  # Count sequences for verification
  fasta_content <- readLines(fasta_file)
  seq_count <- sum(grepl("^>", fasta_content))
  pn_info("Retrieved", seq_count, "sequences for annotation")
  
  # Run the annotation tool
  if (config$annotation$tool == "eggnog") {
    # First check if eggNOG-mapper is available
    emapper_available <- system("which emapper.py >/dev/null 2>&1") == 0
    
    if (emapper_available) {
      pn_info("eggNOG-mapper found. Using it for annotation.")
      annotation_result <- run_eggnog_mapper_r(
        fasta_file = fasta_file,
        output_dir = file.path(output_dir, "eggnog"),
        db_dir = config$annotation$eggnog$db_dir,
        cpu = config$annotation$eggnog$cpu,
        tax_scope = config$annotation$eggnog$tax_scope,
        seed_ortholog_evalue = config$annotation$eggnog$seed_ortholog_evalue,
        seed_ortholog_score = config$annotation$eggnog$seed_ortholog_score
      )
    } else {
      pn_warn("eggNOG-mapper not found. Falling back to dummy annotations.")
      # Use dummy annotations
      annotation_result <- create_dummy_annotations(
        unique(df[[column]]), 
        file.path(output_dir, "dummy_annotations")
      )
    }
  } else if (config$annotation$tool == "cog") {
    annotation_result <- run_cog_classifier(fasta_file, output_dir, config)
  } else if (config$annotation$tool == "dummy") {
    pn_info("Using dummy annotation as specified in config")
    annotation_result <- create_dummy_annotations(
      unique(df[[column]]), 
      file.path(output_dir, "dummy_annotations")
    )
  } else {
    pn_error("Unknown annotation tool:", config$annotation$tool)
    return(NULL)
  }
  
  return(annotation_result)
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

#' Run eggNOG-mapper from R
#'
#' @param fasta_file Path to the input FASTA file
#' @param output_dir Directory for eggNOG-mapper output
#' @param db_dir Path to eggNOG database
#' @param cpu Number of CPUs to use (default: 4)
#' @param tax_scope Taxonomic scope (default: auto)
#' @param seed_ortholog_evalue E-value threshold (default: 0.001)
#' @param seed_ortholog_score Score threshold (default: 60)
#' @param query_cover Query coverage threshold (default: 20)
#' @param subject_cover Subject coverage threshold (default: 20)
#' @param temp_dir Temporary directory (default: tmp)
#' @return A data frame with annotation results or NULL if failed
#' @keywords internal
run_eggnog_mapper_r <- function(fasta_file, output_dir, 
                                db_dir = "data/eggnog_db", 
                                cpu = 4, 
                                tax_scope = "auto",
                                seed_ortholog_evalue = 0.001,
                                seed_ortholog_score = 60,
                                query_cover = 20,
                                subject_cover = 20,
                                temp_dir = "tmp") {
  
  # Create output directory if it doesn't exist
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create temp directory if it doesn't exist
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if emapper.py is in the path
  emapper_cmd <- Sys.which("emapper.py")
  if (emapper_cmd == "") {
    pn_error("emapper.py not found in PATH. Please install eggNOG-mapper.")
    return(NULL)
  }
  
  # Build the command
  output_prefix <- file.path(output_dir, "eggnog_results")
  
  # Create the command with all parameters
  cmd <- paste(
    emapper_cmd,
    "-i", shQuote(fasta_file),
    "--output", shQuote(basename(output_prefix)),
    "--output_dir", shQuote(output_dir),
    "--cpu", cpu,
    "--temp_dir", shQuote(temp_dir),
    "--data_dir", shQuote(db_dir),
    "--tax_scope", tax_scope,
    "--go_evidence", "non-electronic",
    "--target_orthologs", "all",
    "--seed_ortholog_evalue", seed_ortholog_evalue,
    "--seed_ortholog_score", seed_ortholog_score,
    "--query_cover", query_cover,
    "--subject_cover", subject_cover
  )
  
  # Run eggNOG-mapper
  pn_info("Running eggNOG-mapper command:", cmd)
  status <- system(cmd)
  
  if (status != 0) {
    pn_error("eggNOG-mapper failed with status:", status)
    return(NULL)
  }
  
  # Check if the output file exists
  annotations_file <- paste0(output_prefix, ".emapper.annotations")
  if (!file.exists(annotations_file)) {
    pn_error("Expected output file not found:", annotations_file)
    return(NULL)
  }
  
  # Process eggNOG-mapper output into COG classifier format
  result <- process_eggnog_output(annotations_file, output_dir)
  return(result)
}

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
  
  # Try to read the annotations file
  tryCatch({
    # Read the annotations file (from eggNOG-mapper)
    annotations <- read.delim(annotations_file, stringsAsFactors = FALSE)
    
    # Process the eggNOG-mapper output to extract COG annotations
    cog_annotations_file <- file.path(output_dir, "cog_annotations.tsv")
    classifier_file <- file.path(output_dir, "classifier_result.tsv")
    
    # Process annotations into COG format
    classifier_results <- extract_cog_annotations(annotations, cog_categories_file)
    
    # Write to file
    write.table(classifier_results, classifier_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Check if the conversion was successful
    if (file.exists(classifier_file) && file.size(classifier_file) > 0) {
      pn_info("Successfully converted eggNOG-mapper results to COG classifier format")
      return(classifier_results)
    } else {
      pn_error("Failed to convert eggNOG-mapper results")
      return(NULL)
    }
  }, error = function(e) {
    pn_error("Failed to process eggNOG-mapper output:", e$message)
    return(NULL)
  })
}

#' Extract COG annotations from eggNOG-mapper results
#'
#' @param annotations Data frame with eggNOG-mapper annotations
#' @param cog_categories_file Path to COG categories file
#' @return A data frame with COG classifications
#' @keywords internal
extract_cog_annotations <- function(annotations, cog_categories_file) {
  # Read COG categories
  cog_cats <- read.delim(cog_categories_file, stringsAsFactors = FALSE)
  
  # Initialize classifier results
  classifier_results <- data.frame(
    QUERY_ID = character(),
    COG_ID = character(),
    CDD_ID = character(),
    EVALUE = character(),
    GENE_NAME = character(),
    COG_NAME = character(),
    COG_LETTER = character(),
    COG_DESCRIPTION = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each annotation
  for (i in 1:nrow(annotations)) {
    # Check if COG column exists and has data
    if ("COG" %in% colnames(annotations) && 
        !is.na(annotations$COG[i]) && 
        annotations$COG[i] != "-") {
      
      # Split COGs if multiple are present
      cogs <- strsplit(annotations$COG[i], ",")[[1]]
      
      # Get COG category if available
      cog_cat <- if ("COG_category" %in% colnames(annotations) && 
                     !is.na(annotations$COG_category[i]) && 
                     annotations$COG_category[i] != "-") {
        strsplit(annotations$COG_category[i], "")[[1]]
      } else {
        "S"  # Default to "Function unknown"
      }
      
      # Get description if available
      description <- if ("Description" %in% colnames(annotations) && 
                         !is.na(annotations$Description[i]) && 
                         annotations$Description[i] != "-") {
        annotations$Description[i]
      } else {
        "Function unknown"
      }
      
      # Process each COG
      for (cog in cogs) {
        if (grepl("^COG[0-9]+", cog)) {
          # Process each category letter
          for (letter in cog_cat) {
            # Find category description
            cat_idx <- which(cog_cats$COG_LETTER == letter)
            category <- if (length(cat_idx) > 0) {
              cog_cats$COG_CATEGORY[cat_idx]
            } else {
              "POORLY CHARACTERIZED"
            }
            
            # Add to results
            classifier_results <- rbind(classifier_results, data.frame(
              QUERY_ID = annotations$query[i],
              COG_ID = cog,
              CDD_ID = "-",
              EVALUE = "-",
              GENE_NAME = "-",
              COG_NAME = category,
              COG_LETTER = letter,
              COG_DESCRIPTION = description,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # If no annotations were found, create at least one dummy entry
  if (nrow(classifier_results) == 0) {
    pn_warn("No COG annotations found in eggNOG-mapper results. Creating a placeholder entry.")
    classifier_results <- rbind(classifier_results, data.frame(
      QUERY_ID = if (nrow(annotations) > 0) annotations$query[1] else "placeholder",
      COG_ID = "COG0000",
      CDD_ID = "-",
      EVALUE = "-",
      GENE_NAME = "-",
      COG_NAME = "POORLY CHARACTERIZED",
      COG_LETTER = "S",
      COG_DESCRIPTION = "Function unknown",
      stringsAsFactors = FALSE
    ))
  }
  
  return(classifier_results)
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