#!/bin/bash

# Define a list of identifiers
identifiers=(
    "1NKW"
    "1S72"
    "2AW4"
    "2J01"
    "1UN6"
    "1Y26"
    "2AVY"
    "2J00"
    "2ZJR"
    "1NYI"
    "1E8O"
    "1L9A"
    "1LNG"
    "2CZJ"
    "1U8D"
    "2B57"
    "2EES"
    "2HO7"
    "2NZ4"
    "2GDI"
    "2HOJ"
    "2CKY"
    "2QBZ"
    "2OIU"
    "1U6B"
    "1X8W"
    "2A64"
    "1NBS"
    "3E5C"
    "3EGZ"
    "3KTW"
    "3MXH"
    "1FG0"
    "1GID"
    "1J5E"
    "1MFQ"
    "1MMS"
    "1QA6"
    "3D2G"
    "2EET"
)

# Directory to save PDB files
output_dir="./pdb_files/"
mkdir -p "$output_dir"

# Directory to save FORGI files
forgi_output_dir="./forgi_cg/"
mkdir -p "$forgi_output_dir"

# Process each identifier
for id in "${identifiers[@]}"; do
    file_path="${output_dir}${id}.pdb"
    
    # Download PDB file if not already present
    if [[ ! -f "$file_path" ]]; then
        wget -q "https://files.rcsb.org/download/${id}.pdb" -P "$output_dir"
    fi

    # Convert the PDB file to FORGI format
    rnaConvert.py "$file_path" -T forgi --filename "${forgi_output_dir}${id}.cg"
done