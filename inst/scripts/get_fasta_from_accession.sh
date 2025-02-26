#!/bin/bash
# Script to retrieve FASTA sequences from NCBI using accession numbers
# Author: Maximilian Böhm
# Modified: Improved version with better error handling, usage instructions, and performance

# Display help/usage information
display_help() {
    echo "Usage: $0 [options] input_file output_file"
    echo ""
    echo "Description:"
    echo "  This script retrieves FASTA sequences from NCBI using a list of accession numbers."
    echo ""
    echo "Arguments:"
    echo "  input_file    Text file containing one accession number per line"
    echo "  output_file   Output file to save the FASTA sequences"
    echo ""
    echo "Options:"
    echo "  -h, --help    Display this help message and exit"
    echo "  -s, --sleep   Sleep time between requests in seconds (default: 0.34)"
    echo "  -d, --db      NCBI database to query (default: protein)"
    echo "  -e, --email   Your email address for NCBI (optional but recommended)"
    echo "  -b, --batch   Batch size for retrieving multiple sequences at once (default: 50)"
    echo "  -r, --retries Number of retries for failed requests (default: 3)"
    echo ""
    echo "Examples:"
    echo "  $0 accession_list.txt sequences.fasta"
    echo "  $0 --sleep 0.5 --db nucleotide --email user@example.com accession_list.txt sequences.fasta"
    exit 0
}

# Progress bar function (from https://github.com/fearside/ProgressBar/)
function ProgressBar {
    local currentState=$1
    local totalState=$2
    local startTime=$3

    # Process data
    let _progress=(${currentState}*100/${totalState}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
    
    # Build progressbar string lengths
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

    # Calculate elapsed time
    local elapsedTime=$(( $(date +%s) - ${startTime} ))
    
    # Estimate remaining time
    local estimatedTotalTime=$(( elapsedTime * totalState / currentState ))
    local remainingTime=$(( estimatedTotalTime - elapsedTime ))

    # Format time as hours:minutes
    local hours=$(( remainingTime / 3600 ))
    local minutes=$(( (remainingTime % 3600) / 60 ))
    local formattedTime=$(printf "%02d:%02d" $hours $minutes)

    # Build progressbar strings and print the ProgressBar line
    printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%% (${formattedTime} remaining)"
}

# Set default values
SLEEP_TIME=0.34
DATABASE="protein"
EMAIL=""
BATCH_SIZE=50
RETRIES=3
ERROR_LOG="fasta_download_error.log"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            display_help
            ;;
        -s|--sleep)
            SLEEP_TIME="$2"
            shift 2
            ;;
        -d|--db)
            DATABASE="$2"
            shift 2
            ;;
        -e|--email)
            EMAIL="$2"
            shift 2
            ;;
        -b|--batch)
            BATCH_SIZE="$2"
            shift 2
            ;;
        -r|--retries)
            RETRIES="$2"
            shift 2
            ;;
        *)
            if [[ -z "$INPUT_FILE" ]]; then
                INPUT_FILE="$1"
            elif [[ -z "$OUTPUT_FILE" ]]; then
                OUTPUT_FILE="$1"
            else
                echo "Error: Unknown parameter '$1'"
                display_help
            fi
            shift
            ;;
    esac
done

# Check for correct usage
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_FILE" ]; then
  echo "Error: Input file and output file are required." >&2
  display_help
fi

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "Error: File '$INPUT_FILE' not found." >&2
  exit 1
fi

# Check if efetch command is available
if ! command -v efetch &> /dev/null; then
  echo "Error: 'efetch' command not found. Please install NCBI Entrez Direct utilities:" >&2
  echo "  https://www.ncbi.nlm.nih.gov/books/NBK179288/" >&2
  exit 1
fi

# Create or clear the output file
> "$OUTPUT_FILE"

# Clear the error log
> "$ERROR_LOG"

# Print settings
echo "FASTA Retrieval Settings:"
echo "  Input file: $INPUT_FILE"
echo "  Output file: $OUTPUT_FILE"
echo "  Database: $DATABASE"
echo "  Sleep time: $SLEEP_TIME seconds"
echo "  Batch size: $BATCH_SIZE"
echo "  Retries: $RETRIES"
if [ ! -z "$EMAIL" ]; then
  echo "  Email: $EMAIL"
fi
echo ""

# Get the total number of lines in the input file
total_lines=$(wc -l < "$INPUT_FILE")
current_line=0

# Record the start time
start_time=$(date +%s)

# Create a temporary directory for batch processing
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Split the input file into batches
split -l "$BATCH_SIZE" "$INPUT_FILE" "$TEMP_DIR/batch_"

# Process each batch
for batch_file in "$TEMP_DIR"/batch_*; do
    # Create a comma-separated list of IDs for this batch
    BATCH_IDS=$(tr '\n' ',' < "$batch_file" | sed 's/,$//')
    
    # Try to fetch the sequences, with retries
    success=false
    for ((attempt=1; attempt<=RETRIES; attempt++)); do
        if [ ! -z "$EMAIL" ]; then
            efetch -db "$DATABASE" -id "$BATCH_IDS" -format fasta -email "$EMAIL" >> "$OUTPUT_FILE" 2>> "$ERROR_LOG"
        else
            efetch -db "$DATABASE" -id "$BATCH_IDS" -format fasta >> "$OUTPUT_FILE" 2>> "$ERROR_LOG"
        fi
        
        if [ $? -eq 0 ]; then
            success=true
            break
        else
            echo "Retry $attempt/$RETRIES for batch $(basename "$batch_file")..." >&2
            sleep $(echo "$SLEEP_TIME * 2" | bc)  # Double sleep time for retries
        fi
    done
    
    if [ "$success" = false ]; then
        echo "Failed to download batch $(basename "$batch_file") after $RETRIES retries" >&2
        echo "Accessions:" >> "$ERROR_LOG"
        cat "$batch_file" >> "$ERROR_LOG"
        echo "" >> "$ERROR_LOG"
    fi
    
    # Update progress
    lines_in_batch=$(wc -l < "$batch_file")
    current_line=$((current_line + lines_in_batch))
    ProgressBar ${current_line} ${total_lines} ${start_time}
    
    # Add a delay to avoid overwhelming the server
    sleep "$SLEEP_TIME"
done

# Check the output file
if [ ! -s "$OUTPUT_FILE" ]; then
    echo -e "\nWarning: Output file is empty. No sequences were downloaded." >&2
    exit 1
fi

# Count sequences in the output file
seq_count=$(grep -c "^>" "$OUTPUT_FILE")

printf '\nFinished!\n'
echo "Downloaded $seq_count sequences out of $total_lines accessions"

# Check for errors
if [ -s "$ERROR_LOG" ]; then
    echo "Some errors occurred during download. See $ERROR_LOG for details."
    echo "Failed to download $(($total_lines - $seq_count)) sequences."
fi

exit 0
