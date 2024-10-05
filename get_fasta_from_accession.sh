#!/bin/bash
# found online: https://github.com/fearside/ProgressBar/
# 1. Create ProgressBar function
# 1.1 Input is currentState($1), totalState($2), and startTime($3)
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

    # 1.2 Build progressbar strings and print the ProgressBar line
    # 1.2.1 Output example:                           
    # 1.2.1.1 Progress : [########################################] 100% (00:30 remaining)
    printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%% (${formattedTime} remaining)"
}

# Check for correct usage
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 input_file output_file" >&2
  exit 1
fi

input_file="$1"
output_file="$2"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "Error: File '$input_file' not found." >&2
  exit 1
fi

# Clear the output file if it exists
> "$output_file"

# Get the total number of lines in the input file
total_lines=$(wc -l < "$input_file")
current_line=0

# Record the start time
start_time=$(date +%s)

# Read each line from the input file and feed it to the command
while IFS= read -r line; do
  # Replace 'your_command' with the actual command you want to run
  # For example, if you want to use 'efetch' with each line as an ID:
  efetch -db protein -id "$line" -format fasta >> "$output_file" 2>> fasta_download_error.log
  sleep 0.34  # Add a delay to avoid overwhelming the server

  # Update the progress bar
  ((current_line++))
  ProgressBar ${current_line} ${total_lines} ${start_time}
done < "$input_file"

printf '\nFinished!\n'