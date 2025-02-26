#!/bin/bash
# Script to run eggNOG-mapper on a FASTA file and format the output for the
# protein neighbors analysis package
#
# Author: Maximilian Böhm
# Date: February 2025

# Display help message
display_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required options:"
    echo "  -i, --input FILE       Input FASTA file"
    echo "  -o, --output DIR       Output directory"
    echo ""
    echo "Optional parameters:"
    echo "  -d, --db_dir DIR       eggNOG database directory (default: /data/eggnog_db)"
    echo "  -c, --cpu INT          Number of CPU cores to use (default: 4)"
    echo "  -t, --tax_scope STR    Taxonomic scope (default: auto)"
    echo "  -g, --go_evidence STR  GO evidence level (default: non-electronic)"
    echo "  -s, --seed_orthologs   Seed ortholog options: all (default), one2one, many2one, one2many, many2many"
    echo "  -e, --evalue FLOAT     E-value threshold (default: 0.001)"
    echo "  -q, --query_cov FLOAT  Query coverage threshold (default: 20)"
    echo "  -j, --subject_cov FLOAT Subject coverage threshold (default: 20)"
    echo "  -h, --help             Display this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 -i proteins.fasta -o results_dir -c 8 -d /path/to/eggnog_db"
    exit 0
}

# Progress bar function for better user feedback
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

    local hours=$(( remainingTime / 3600 ))
    local minutes=$(( (remainingTime % 3600) / 60 ))
    local formattedTime=$(printf "%02d:%02d" $hours $minutes)

    # Build progressbar strings and print the ProgressBar line
    printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%% (${formattedTime} remaining)"
}

# Default values
DB_DIR="/data/eggnog_db"
CPU=4
TAX_SCOPE="auto"
GO_EVIDENCE="non-electronic"
TARGET_ORTHOLOGS="all"
SEED_ORTHOLOG_EVALUE=0.001
SEED_ORTHOLOG_SCORE=60
QUERY_COVERAGE=20
SUBJECT_COVERAGE=20
TEMP_DIR="tmp"
INPUT_FASTA=""
OUTPUT_DIR=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input)
            INPUT_FASTA="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--db_dir)
            DB_DIR="$2"
            shift 2
            ;;
        -c|--cpu)
            CPU="$2"
            shift 2
            ;;
        -t|--tax_scope)
            TAX_SCOPE="$2"
            shift 2
            ;;
        -g|--go_evidence)
            GO_EVIDENCE="$2"
            shift 2
            ;;
        -s|--seed_orthologs)
            TARGET_ORTHOLOGS="$2"
            shift 2
            ;;
        -e|--evalue)
            SEED_ORTHOLOG_EVALUE="$2"
            shift 2
            ;;
        -q|--query_cov)
            QUERY_COVERAGE="$2"
            shift 2
            ;;
        -j|--subject_cov)
            SUBJECT_COVERAGE="$2"
            shift 2
            ;;
        -h|--help)
            display_help
            ;;
        *)
            echo "Unknown option: $1"
            display_help
            ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Input FASTA and output directory are required."
    display_help
fi

# Check if the input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create temp directory if it doesn't exist
mkdir -p "$TEMP_DIR"

# Record the start time
start_time=$(date +%s)

echo "Running eggNOG-mapper on $INPUT_FASTA with output to $OUTPUT_DIR"
echo "Parameters:"
echo "  Database directory: $DB_DIR"
echo "  CPUs: $CPU"
echo "  Taxonomic scope: $TAX_SCOPE"
echo "  E-value: $SEED_ORTHOLOG_EVALUE"
echo "  Score: $SEED_ORTHOLOG_SCORE"
echo "  Query coverage: $QUERY_COVERAGE"
echo "  Subject coverage: $SUBJECT_COVERAGE"
echo "  Target orthologs: $TARGET_ORTHOLOGS"

# Check if emapper.py is installed
if ! command -v emapper.py &> /dev/null; then
    echo "Error: emapper.py is not installed or not in PATH. Please install eggNOG-mapper."
    exit 1
fi

# Run eggNOG-mapper
echo "Starting eggNOG-mapper analysis..."
emapper.py \
    -i "$INPUT_FASTA" \
    --output "$OUTPUT_DIR/eggnog_results" \
    --output_dir "$OUTPUT_DIR" \
    --cpu $CPU \
    --temp_dir "$TEMP_DIR" \
    --data_dir "$DB_DIR" \
    --tax_scope $TAX_SCOPE \
    --go_evidence $GO_EVIDENCE \
    --target_orthologs $TARGET_ORTHOLOGS \
    --seed_ortholog_evalue $SEED_ORTHOLOG_EVALUE \
    --seed_ortholog_score $SEED_ORTHOLOG_SCORE \
    --query_cover $QUERY_COVERAGE \
    --subject_cover $SUBJECT_COVERAGE

# Check if the annotation was successful
if [ $? -ne 0 ]; then
    echo "Error: eggNOG-mapper failed to run properly."
    exit 1
fi

# Process the output to extract COG annotations
echo "Processing eggNOG-mapper results to extract COG annotations..."

# Create a file with COG category information
cat > "$OUTPUT_DIR/cog_categories.tsv" << 'EOF'
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
Q	METABOLISM	Secondary metabolites biosynthesis, transport and catabolism
R	POORLY CHARACTERIZED	General function prediction only
S	POORLY CHARACTERIZED	Function unknown
EOF

# Extract COG annotations from eggNOG-mapper output
echo "Extracting COG annotations from eggNOG-mapper results..."
awk -F'\t' 'NR>1 && $8 != "-" {
    split($8, cogs, ",");
    for (i in cogs) {
        if (match(cogs[i], /^COG[0-9]+/)) {
            cog_id = cogs[i];
            cog_category = (length($9) > 0 && $9 != "-") ? $9 : "S";
            cog_desc = (length($10) > 0 && $10 != "-") ? $10 : "Function unknown";
            print $1 "\t" cog_id "\t" cog_category "\t" cog_desc;
        }
    }
}' "$OUTPUT_DIR/eggnog_results.emapper.annotations" > "$OUTPUT_DIR/cog_annotations.tsv"

# Create a COG classifier-like output for compatibility with the existing pipeline
echo "Creating COG classifier-compatible output format..."
(echo -e "QUERY_ID\tCOG_ID\tCDD_ID\tEVALUE\tGENE_NAME\tCOG_NAME\tCOG_LETTER\tCOG_DESCRIPTION"
 awk -F'\t' '{
     query = $1;
     cog_id = $2;
     cog_letter = $3;
     cog_desc = $4;
     
     # Look up the COG category description
     cmd = "grep \"^" cog_letter "\\t\" \"" ENVIRON["OUTPUT_DIR"] "/cog_categories.tsv\"";
     cmd | getline coginfo;
     split(coginfo, cogparts, "\t");
     cog_category = cogparts[2];
     
     print query "\t" cog_id "\t-\t-\t-\t" cog_category "\t" cog_letter "\t" cog_desc;
     
 }' "$OUTPUT_DIR/cog_annotations.tsv") > "$OUTPUT_DIR/classifier_result.tsv"

# Check if the conversion was successful
if [ -s "$OUTPUT_DIR/classifier_result.tsv" ]; then
    echo "Successfully converted eggNOG-mapper results to COG classifier format."
    # Count the number of annotations
    annotation_count=$(grep -v "^QUERY_ID" "$OUTPUT_DIR/classifier_result.tsv" | wc -l)
    echo "Generated $annotation_count COG annotations for the proteins."
else
    echo "Warning: Conversion to COG classifier format produced an empty file."
fi

# Calculate elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "Analysis complete in ${hours}h ${minutes}m ${seconds}s."
echo "Results are available in $OUTPUT_DIR"

# Output file information
result_size=$(du -h "$OUTPUT_DIR/classifier_result.tsv" | cut -f1)
echo "COG classifier format results ($result_size): $OUTPUT_DIR/classifier_result.tsv"
echo "Raw eggNOG-mapper results: $OUTPUT_DIR/eggnog_results.emapper.annotations"

# Generate a brief summary report
echo "Summary of findings:" > "$OUTPUT_DIR/annotation_summary.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/annotation_summary.txt"
echo "Input file: $INPUT_FASTA" >> "$OUTPUT_DIR/annotation_summary.txt"
echo "Total proteins analyzed: $(grep -c "^>" "$INPUT_FASTA")" >> "$OUTPUT_DIR/annotation_summary.txt"
echo "Proteins with COG annotations: $annotation_count" >> "$OUTPUT_DIR/annotation_summary.txt"
echo "COG category distribution:" >> "$OUTPUT_DIR/annotation_summary.txt"
awk -F'\t' 'NR>1 {count[$7]++} END {for (letter in count) printf "%s\t%d\n", letter, count[letter]}' "$OUTPUT_DIR/classifier_result.tsv" | sort -k2,2nr >> "$OUTPUT_DIR/annotation_summary.txt"

echo "Generated summary report: $OUTPUT_DIR/annotation_summary.txt"
