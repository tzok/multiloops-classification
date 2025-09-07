#!/bin/bash

# Define directories
input_folder="./forgi_cg"
output_folder="./forgi_graph_files"
script_path="./parse_forgi.py"

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Process each file in the input folder
for input_file in "$input_folder"/*.cg; do
    # Extract the base filename without extension
    base_name=$(basename "$input_file" .cg)
    
    # Define the output file path
    output_file="$output_folder/${base_name}.json"
    
    # Run the Python script
    python "$script_path" "$input_file" "$output_file"
done