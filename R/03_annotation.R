#' Analyze proteins using eggNOG or COG classifier
#'
#' This function annotates proteins using either eggNOG-mapper or COGclassifier.
#' It extracts the FASTA sequences for the proteins and runs the selected annotation tool.
#'
#' @param df A data frame with protein information.
#' @param column The column name containing protein IDs.
#' @param config The configuration list.
#' @return A data frame with the annotation results.
#' @export
analyze_proteins <- function(df, column = 'ID', config) {
  current_date <- config$analysis$date
  output_dir <- file.path(config$paths$output_dir, current_date)
  
  # Log the start of protein annotation
  pn_info("Starting protein annotation with", config$annotation$tool)
  pn_info(paste("Annotating", nrow(df), "proteins"))
  
  # Export accessions of the proteins
  protein_accessions <- data.frame(df %>% select(all_of(column)))
  accessions_file <- file.path(output_dir, "protein_accessions.csv")
  write_csv(protein_accessions, accessions_file, col_names = FALSE)
  
  # Get the FASTA file of the proteins
  fasta_file <- file.path(output_dir, "proteins.fasta")
  pn_info("Retrieving FASTA sequences for proteins")
  
  # Run the FASTA retrieval script
  fasta_script <- system.file("scripts", "get_fasta_from_accession.sh", package = "proteinNeighbours")
  if (!file.exists(fasta_script)) {
    fasta_script <- "get_fasta_from_accession.sh"  # Fallback to current directory
  }
  
  bash_command <- paste(fasta_script, accessions_file, fasta_file)
  pn_info("Running command:", bash_command)
  
  system_result <- system2(bash_command, wait = TRUE)
  if (system_result != 0) {
    pn_error("Failed to retrieve FASTA sequences")
    return(NULL)
  }
  
  # Run the annotation tool
  if (config$annotation$tool == "eggnog") {
    annotation_result <- run_eggnog_mapper(fasta_file, output_dir, config)
  } else if (config$annotation$tool == "cog") {
    annotation_result <- run_cog_classifier(fasta_file, output_dir, config)
  } else {
    pn_error("Unknown annotation tool:", config$annotation$tool)
    return(NULL)
  }
  
  return(annotation_result)
}

#' Run eggNOG-mapper for protein annotation
#'
#' This function runs eggNOG-mapper on a FASTA file to annotate proteins.
#'
#' @param fasta_file Path to the FASTA file.
#' @param output_dir Output directory for results.
#' @param config The configuration list.
#' @return A data frame with the annotation results.
#' @keywords internal
run_eggnog_mapper <- function(fasta_file, output_dir, config) {
  pn_info("Running eggNOG-mapper for protein annotation")
  
  # Create eggNOG output directory
  eggnog_dir <- file.path(output_dir, "eggnog")
  dir.create(eggnog_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build command for eggNOG-mapper
  eggnog_script <- system.file("scripts", "run_eggnog_mapper.sh", package = "proteinNeighbours")
  if (!file.exists(eggnog_script)) {
    eggnog_script <- "inst/scripts/run_eggnog_mapper.sh"  # Fallback to current directory
  }
  
  # Extract eggNOG parameters from config
  eggnog_config <- config$annotation$eggnog
  
  # Build the command
  cmd_parts <- c(
    eggnog_script,
    "-i", fasta_file,
    "-o", eggnog_dir,
    "-d", eggnog_config$db_dir,
    "-c", eggnog_config$cpu,
    "-t", eggnog_config$tax_scope,
    "-e", eggnog_config$seed_ortholog_evalue,
    "-s", eggnog_config$seed_ortholog_score
  )
  
  cmd <- paste(cmd_parts, collapse = " ")
  pn_info("Running command:", cmd)
  
  # Run eggNOG-mapper
  system_result <- system2(cmd, wait = TRUE)
  if (system_result != 0) {
    pn_error("Failed to run eggNOG-mapper")
    return(NULL)
  }
  
  # Read the results
  result_file <- file.path(eggnog_dir, "classifier_result.tsv")
  if (!file.exists(result_file)) {
    pn_error("eggNOG-mapper result file not found:", result_file)
    return(NULL)
  }
  
  # Read the annotation results
  cog_classification <- read.delim(result_file)
  
  # Log success
  pn_info("Successfully annotated proteins with eggNOG-mapper")
  pn_info(paste("Found", nrow(cog_classification), "annotations"))
  
  return(cog_classification)
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
  
  # Run COGclassifier
  cog_command <- paste("COGclassifier -i", fasta_file, "-o", cog_dir)
  pn_info("Running command:", cog_command)
  
  system_result <- system2(cog_command, wait = TRUE)
  if (system_result != 0) {
    pn_error("Failed to run COG classifier")
    return(NULL)
  }
  
  # Read the results
  result_file <- file.path(cog_dir, "classifier_result.tsv")
  if (!file.exists(result_file)) {
    pn_error("COG classifier result file not found:", result_file)
    return(NULL)
  }
  
  cog_classification <- read.delim(result_file)
  
  # Log success
  pn_info("Successfully annotated proteins with COG classifier")
  pn_info(paste("Found", nrow(cog_classification), "annotations"))
  
  return(cog_classification)
}


