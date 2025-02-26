#' Neighbor Identification Functions
#'
#' This file contains functions for identifying neighboring proteins
#' in genomic data, parsing GFF files, and extracting attributes.
#'
#' @author Maximilian Böhm
#' @modified Improved version with enhanced documentation and error handling

#' Extract a specific field from the attribute column of a GFF file
#'
#' This function parses the GFF attribute field and extracts a specific attribute.
#'
#' @param x The attribute column as a character vector.
#' @param field The field to extract.
#' @param attrsep The separator used in the attribute column (default is ";").
#' @return A character vector with the extracted field values.
#' @export
getAttributeField <- function(x, field, attrsep = ";") {
  # Split the attributes
  s <- strsplit(x, split = attrsep, fixed = TRUE)
  
  # Extract the desired field
  result <- sapply(s, function(atts) {
    a <- strsplit(atts, split = "=", fixed = TRUE)
    m <- match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv <- a[[m]][2]
    } else {
      rv <- as.character(NA)
    }
    return(rv)
  })
  
  return(result)
}

#' Find the protein alias
#'
#' This function looks up a protein ID in the alias database and returns its PIGI.
#'
#' @param protein.id The protein ID to search for.
#' @param alias The alias database.
#' @param verbose Whether to print additional information (default is FALSE).
#' @param identi Whether to include identity information (default is FALSE).
#' @return A data frame with the corresponding PIGI and alias.
#' @export
protein.alias <- function(protein.id, alias, verbose = FALSE, identi = FALSE) {
  # Check if the protein ID is in the alias database
  if (protein.id %in% alias$alias) {
    # If it is, find the corresponding PIGI
    result <- alias %>%
      dplyr::distinct(PIGI, alias) %>%
      dplyr::filter(alias == protein.id) %>%
      dplyr::select(PIGI)
    
    # Also find the corresponding identity
    if (identi) {
      percent <- alias %>%
        dplyr::filter(alias == protein.id) %>%
        dplyr::select(identity)
    }
    
    # Print a message indicating the protein was found, if verbose is TRUE
    if (verbose) {
      if (identi) {
        log_info(paste0("Protein (", protein.id, ") was found in alias database, with an identity of ", 
                      percent, " matching ", result, "."))
      } else {
        log_info(paste0("Protein (", protein.id, ") was found in alias database matching ", result, "."))
      }
    }
    
    # Return the PIGI as a data frame
    return(result)
  } else {
    # If the protein ID is not in the database, print a message indicating it was not found, if verbose is TRUE
    if (verbose) {
      log_warn("Protein was not found in alias database:", protein.id)
    }
    
    # Return the original protein ID as a data frame
    return(data.frame(PIGI = protein.id, stringsAsFactors = FALSE))
  }
}

#' Get protein information from GFF file
#' 
#' This function extracts information about a specific protein from a GFF file.
#'
#' @param input The path to the GFF file.
#' @param protein.id The protein ID to search for.
#' @return A data frame with the protein information.
#' @export
getProteinInfoFromGff <- function(input, protein.id) {
  # Inform the user which input is being processed
  log_debug(paste("Processing input:", input))
  
  # Check if the file exists
  if (!file.exists(input)) {
    log_warn(paste("File not found:", input, "- Skipping this input."))
    return(NULL)
  }
  
  # Import GFF file
  tryCatch({
    gffData <- ape::read.gff(input, na.strings = c(".", "?"), GFF3 = TRUE)
    log_debug("Successfully read GFF file:", input)
  }, error = function(e) {
    log_error("Failed to read GFF file:", input, "-", e$message)
    return(NULL)
  })
  
  # Get Protein ID and product separately
  gffData <- gffData %>%
    dplyr::mutate(ID = gsub(".*\\-", "", getAttributeField(attributes, "ID")),
           product = getAttributeField(attributes, "product"))
  
  # Get protein information
  protein_info <- gffData %>% dplyr::filter(ID == protein.id)
  
  if (nrow(protein_info) == 0) {
    log_warn("Protein ID not found in GFF file:", protein.id)
  } else {
    log_debug("Found protein information for:", protein.id)
  }
  
  return(protein_info)
}

#' Collect all protein information
#' 
#' This function collects information about all proteins in the protein.assembly data frame.
#'
#' @param protein.assembly A data frame with protein and assembly information.
#' @param PATH Path where the ncbi_dataset folder is stored (default = data).
#' @return A data frame with all protein information.
#' @export
collect_all_protein_info <- function(protein.assembly, PATH = 'data') {
  log_info("Collecting information for", nrow(protein.assembly), "proteins")
  
  # Initialize empty data frame for all protein information
  all.protein.info <- data.frame()
  
  # For progress tracking
  total_proteins <- nrow(protein.assembly)
  progress_step <- max(1, floor(total_proteins / 20))  # Show progress roughly every 5%
  
  # Process each protein
  for (i in seq_len(total_proteins)) {
    # Show progress
    if (i %% progress_step == 0 || i == total_proteins) {
      progress_pct <- round(i / total_proteins * 100)
      log_info(paste0("Processing protein ", i, " of ", total_proteins, " (", progress_pct, "%)"))
    }
    
    # Construct GFF file path
    gff_file <- file.path(PATH, 'ncbi_dataset/data', protein.assembly[i, 1], 'genomic.gff')
    
    # Get protein information
    pi <- getProteinInfoFromGff(gff_file, protein.assembly[i, 2])
    
    # Skip the iteration if pi is NULL
    if (is.null(pi)) {
      log_warn("Skipping protein", protein.assembly[i, 2], "with assembly", protein.assembly[i, 1])
      next
    }
    
    # Add protein and assembly information to all protein information data frame
    pi <- pi %>%
      dplyr::mutate(PIGI = protein.assembly[i, 2],  # Protein I Gave as Input
                    assembly = protein.assembly[i, 1])
    
    # Add protein information to all protein information data frame
    all.protein.info <- rbind(all.protein.info, pi)
  }
  
  # Mark these as non-neighbor proteins
  all.protein.info$is.neighbour <- FALSE
  
  # Get the current date in YYYY-MM-DD format
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  
  # Save results of this long run
  output_file <- file.path('output', current_date, 'all_protein_info.csv')
  readr::write_csv(all.protein.info, output_file)
  log_info("Saved protein information to:", output_file)
  
  return(all.protein.info)
}

#' Get neighboring proteins
#'
#' This function identifies neighboring proteins in a GFF data frame.
#'
#' @param gff.df The GFF data frame.
#' @param protein.id The protein ID to search for neighbors.
#' @param bp Number of base pairs to consider for neighbors (default is 300).
#' @param n Number of neighbors to find (default is 15).
#' @return A data frame with neighboring proteins.
#' @export
getNeighborProteins <- function(gff.df, protein.id, bp = 300, n = 15) {
  log_debug("Finding neighbors for protein:", protein.id)
  
  # Initialize empty data frame for neighbors
  neighbors <- data.frame()
  
  # Get start, end, strand, and seqid for the given protein id
  protein_info <- gff.df %>% dplyr::filter(ID == protein.id)
  
  if (nrow(protein_info) == 0) {
    log_warn("Protein ID not found in GFF data:", protein.id)
    return(neighbors)
  }
  
  start <- protein_info$start
  end <- protein_info$end
  strand <- protein_info$strand
  seqid <- protein_info$seqid
  
  log_debug(paste("Protein location: seqid =", seqid, ", start =", start, 
                ", end =", end, ", strand =", strand))
  
  # Loop through n upstream and downstream neighbors
  for (i in seq_len(n)) {
    # Define condition for upstream neighbors
    condition <- start - bp < gff.df$end & 
                gff.df$end < start & 
                gff.df$type == 'CDS' & 
                gff.df$strand == strand & 
                gff.df$seqid == seqid
    
    # If there are upstream neighbors, add the first one to the neighbors data frame
    if (sum(condition) > 0) {
      output <- gff.df[condition, ][1, ]
      start <- output$start
      neighbors <- rbind(neighbors, output)
      log_debug(paste("Found upstream neighbor", i, ":", output$ID))
    }
    
    # Define condition for downstream neighbors
    condition <- end + bp > gff.df$start & 
                gff.df$start > end & 
                gff.df$type == 'CDS' & 
                gff.df$strand == strand & 
                gff.df$seqid == seqid
    
    # If there are downstream neighbors, add the first one to the neighbors data frame
    if (sum(condition) > 0) {
      output <- gff.df[condition, ][1, ]
      end <- output$end
      neighbors <- rbind(neighbors, output)
      log_debug(paste("Found downstream neighbor", i, ":", output$ID))
    }
  }
  
  log_info("Found", nrow(neighbors), "neighbors for protein", protein.id)
  return(neighbors)
}

#' Get protein neighbors from a GFF3 file
#'
#' This function reads a GFF3 file and identifies neighboring proteins.
#'
#' @param input The path to the GFF3 file.
#' @param protein.id The protein ID to search for neighbors.
#' @param basepairs Number of base pairs to consider for neighbors (default is 300).
#' @param m Number of neighbors to find (default is 15).
#' @return A data frame with neighboring proteins.
#' @export
getProteinNeighborsFromGff3 <- function(input, protein.id, basepairs = 300, m = 15) {
  # Inform the user which input is being processed
  log_info(paste("Processing input:", input))
  
  # Check if the file exists
  if (!file.exists(input)) {
    log_warn(paste("File not found:", input, "- Skipping this input."))
    return(NULL)
  }
  
  # Import GFF3 file
  tryCatch({
    gffData <- ape::read.gff(input, na.strings = c(".", "?"), GFF3 = TRUE)
    log_debug("Successfully read GFF file:", input)
  }, error = function(e) {
    log_error("Failed to read GFF file:", input, "-", e$message)
    return(NULL)
  })
  
  # Get Protein ID and product separately
  gffData <- gffData %>%
    dplyr::mutate(ID = gsub(".*\\-", "", getAttributeField(attributes, "ID")),
                 product = getAttributeField(attributes, "product"))
  
  # Get neighboring proteins
  finalOutput <- getNeighborProteins(gffData, protein.id, basepairs, m)
  
  return(finalOutput)
}

#' Collect all neighbors for a list of proteins
#'
#' This function identifies neighbors for all proteins in the protein.assembly data frame.
#'
#' @param protein.assembly A data frame with protein and assembly information.
#' @param basepairs Number of base pairs to consider for neighbors (default is 300).
#' @param max_neighbors Number of neighbors to find (default is 15).
#' @param PATH Path where the ncbi_dataset folder is stored (default = data).
#' @return A data frame with all neighbors.
#' @export
collect_all_neighbours <- function(protein.assembly, basepairs = 300, max_neighbors = 15, PATH = 'data') {
  log_info("Collecting neighbors for", nrow(protein.assembly), "proteins")
  log_info(paste("Parameters: basepairs =", basepairs, ", max_neighbors =", max_neighbors))
  
  # Initialize empty data frame for all neighbors
  all.neighbours <- data.frame()
  
  # For progress tracking
  total_proteins <- nrow(protein.assembly)
  progress_step <- max(1, floor(total_proteins / 20))  # Show progress roughly every 5%
  
  # Process each protein
  for (i in seq_len(total_proteins)) {
    # Show progress
    if (i %% progress_step == 0 || i == total_proteins) {
      progress_pct <- round(i / total_proteins * 100)
      log_info(paste0("Processing protein ", i, " of ", total_proteins, " (", progress_pct, "%)"))
    }
    
    # Construct GFF file path
    gff_file <- file.path(PATH, 'ncbi_dataset/data', protein.assembly[i, 1], 'genomic.gff')
    
    # Get neighbors for this protein
    np <- getProteinNeighborsFromGff3(gff_file, protein.assembly[i, 2], basepairs, max_neighbors)
    
    # Skip the iteration if np is NULL
    if (is.null(np)) {
      log_warn("No neighbors found for protein", protein.assembly[i, 2], "with assembly", protein.assembly[i, 1])
      next
    }
    
    # Add neighbors to the output if found
    if (nrow(np) > 0) {
      # Add protein and assembly information to neighbors data frame
      np <- np %>%
        dplyr::mutate(PIGI = protein.assembly[i, 2],  # Protein I Gave as Input
                     assembly = protein.assembly[i, 1])
      
      # Add neighbors to all neighbors data frame
      all.neighbours <- rbind(all.neighbours, np)
      log_debug(paste("Added", nrow(np), "neighbors for protein", protein.assembly[i, 2]))
    } else {
      log_debug(paste("No neighbors found for protein", protein.assembly[i, 2]))
    }
  }
  
  # Mark these as neighbor proteins
  all.neighbours$is.neighbour <- TRUE
  
  # Get the current date in YYYY-MM-DD format
  current_date <- format(Sys.Date(), "%Y-%m-%d")
  
  # Save results of this long run
  output_file <- file.path('output', current_date, 
                          paste0('all_neighbours_bp', basepairs, '_n', max_neighbors, '.csv'))
  readr::write_csv(all.neighbours, output_file)
  log_info("Saved neighbor information to:", output_file)
  
  return(all.neighbours)
}

#' Concatenate columns in a data frame
#'
#' This utility function concatenates specified columns into a new column.
#'
#' @param df The data frame.
#' @param columns_vector Vector of column names to concatenate.
#' @param new_column_name Name for the new column.
#' @return Data frame with the new concatenated column.
#' @keywords internal
concatenate_columns <- function(df, columns_vector, new_column_name = "concatenated") {
  # Concatenate specified columns into a single string separated by a comma
  df %>%
    tidyr::unite(new_column_name, dplyr::all_of(columns_vector), sep = ",", remove = FALSE)
}

#' Convert data frame to FASTA format
#'
#' This function creates a FASTA file from a data frame with sequence information.
#'
#' @param df The data frame containing sequence data.
#' @param name_vector Vector of column names to use for sequence names.
#' @param sequence_column Name of the column containing sequences.
#' @param output_file Path for the output FASTA file.
#' @return Invisible NULL, called for side effects.
#' @export
data.frame2fasta <- function(df, name_vector, sequence_column = "seq", output_file = "fasta.fasta") {
  log_info("Converting data frame to FASTA format")
  log_info(paste("Output file:", output_file))
  
  # Generate new df with a name column and a sequence column
  df2 <- concatenate_columns(df, name_vector, new_column_name = "name")
  
  df2 <- df2 %>%
    dplyr::select(dplyr::all_of(c("name", sequence_column))) %>%
    dplyr::mutate(name = paste(">", name))
  
  # Create FASTA format
  tryCatch({
    fasta_lines <- do.call(rbind, lapply(seq(nrow(df2)), function(i) t(df2[i, ])))
    write.table(fasta_lines, row.names = FALSE, col.names = FALSE, quote = FALSE, file = output_file)
    log_info("Successfully wrote FASTA file with", nrow(df2), "sequences")
  }, error = function(e) {
    log_error("Failed to write FASTA file:", e$message)
  })
  
  return(invisible(NULL))
}

#' Plot neighbors of a specific protein
#'
#' This function plots the genomic neighbors of a specific protein.
#'
#' @param all.neighbours.df A data frame with all neighbors.
#' @param protein.id The protein ID to plot.
#' @param output_dir Directory to save the plot (default is "output").
#' @param width Plot width in cm (default is 20).
#' @param height Plot height in cm (default is 5).
#' @return Invisible NULL, called for side effects.
#' @export
plot_neighbours <- function(all.neighbours.df, protein.id, output_dir = "output", 
                           width = 20, height = 5) {
  log_info("Plotting neighbors for protein:", protein.id)
  
  # Check if protein_id exists in the neighbors data frame
  if (!(protein.id %in% all.neighbours.df$PIGI)) {
    log_warn("Protein ID not found in neighbors data frame:", protein.id)
    return(invisible(NULL))
  }
  
  # Get assembly for given protein id
  assembly <- all.neighbours.df$assembly[all.neighbours.df$PIGI == protein.id][1]
  log_debug("Using assembly:", assembly)
  
  # Construct GFF file path
  gff_file <- file.path('data/ncbi_dataset/data', assembly, 'genomic.gff')
  
  # Check if the GFF file exists
  if (!file.exists(gff_file)) {
    log_error("GFF file not found:", gff_file)
    return(invisible(NULL))
  }
  
  # Read GFF data
  tryCatch({
    df <- ape::read.gff(gff_file, na.strings = c(".", "?"), GFF3 = TRUE)
    log_debug("Successfully read GFF file for plotting")
  }, error = function(e) {
    log_error("Failed to read GFF file for plotting:", e$message)
    return(invisible(NULL))
  })
  
  # Get Protein ID and Name separately
  df <- df %>%
    dplyr::mutate(ID = gsub(".*\\-", "", getAttributeField(attributes, "ID")),
                 Name = getAttributeField(attributes, "Name"))
  
  # Find the target protein in the GFF data
  target_protein <- df %>% dplyr::filter(ID == protein.id)
  if (nrow(target_protein) == 0) {
    log_warn("Protein ID not found in GFF file:", protein.id)
    return(invisible(NULL))
  }
  
  # Create data frame for the target protein
  df0 <- data.frame(molecule = target_protein$seqid[1],
                    gene = protein.id, 
                    start = target_protein$start[1],
                    end = target_protein$end[1],
                    strand = target_protein$strand[1],
                    orientation = 1)
  
  # Filter neighbors for this protein
  df1 <- all.neighbours.df %>%
    dplyr::filter(PIGI == protein.id) %>%
    dplyr::select(molecule = seqid, gene = product, start, end, strand, orientation = 1)
  
  # Combine target protein and neighbors
  df2 <- rbind(df0, df1)
  
  # Plot neighbors
  tryCatch({
    p <- ggplot2::ggplot(df2, ggplot2::aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
      gggenes::geom_gene_arrow() +
      ggplot2::facet_wrap(~ molecule, scales = "free", ncol = 1) +
      gggenes::theme_genes()
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save the plot
    output_file <- file.path(output_dir, paste0('neighbours_', protein.id, '.png'))
    ggplot2::ggsave(output_file, plot = p, device = "png", width = width, height = height, units = "cm")
    log_info("Saved neighbor plot to:", output_file)
    
  }, error = function(e) {
    log_error("Failed to create or save neighbor plot:", e$message)
  })
  
  return(invisible(NULL))
}
