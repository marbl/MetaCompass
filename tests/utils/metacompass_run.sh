#!/bin/bash

set -x

# Define a Bash function to run Metacompass
run_metacompass() {
  local NXF_OPTS="-Dleveldb.mmap=false -Xms500m -Xmx2g"
  export NXF_OPTS

  local output_read_folder="$1"
  local metacompass_path="$2"
  local reference_db="$3"
  local forward_read="$4"
  local reverse_read="$5"

  nextflow -log "${output_read_folder}/nextflow.log" run "$metacompass_path" \
    --reference_db "$reference_db" \
    --forward "$forward_read" \
    --reverse "$reverse_read" \
    --output "$output_read_folder" \
    --threads 1 \
    --trace_file_name "$output_read_folder/trace.txt" \
    -with-timeline "$output_read_folder/timeline.html" \
    -with-dag "$output_read_folder/dag.png"
}

# Example usage:
# run_metacompass
# "/path/to/output_folder"
# "/path/to/metacompass_script.nf"
# "/path/to/reference_db" "/path/to/forward_read" "/path/to/reverse_read"

# Call this function with command-line arguments
run_metacompass "$@"
