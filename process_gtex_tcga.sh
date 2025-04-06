#!/bin/bash

process_tissue_data() {
    local tissue=$1
    
    if [ -z "$tissue" ]; then
        echo "Error: No tissue type specified"
        echo "Usage: $0 <tissue_type>"
        exit 1
    fi
    
    echo "Processing $tissue GTEX-TCGA data..."
    
    # Convert tissue to uppercase for comparison labels
    local tissue_upper=$(echo "$tissue" | tr '[:lower:]' '[:upper:]')
    
    # File names based on tissue type
    local base_file="${tissue}-gtex-tcga.txt"
    local labels_file="${tissue}-gtex-tcga.labels.txt"
    
    # Check if input files exist
    if [ ! -f "$base_file" ]; then
        echo "Error: Input file $base_file not found"
        exit 1
    fi
    
    if [ ! -f "$labels_file" ]; then
        echo "Error: Labels file $labels_file not found"
        exit 1
    fi
    
    # convert RSEM decimal values to integers
    echo "Converting decimal values to integers..."
    cat "$base_file" | sed -e 's/\.[0-9]*//g' -e 's/ *$//' > "${tissue}-gtex-tcga-integer.txt"
    
    # remove duplicate gene rows
    echo "Removing duplicate gene rows..."
    cat "${tissue}-gtex-tcga-integer.txt" | awk '!a[$1]++' > "${tissue}-gtex-tcga-integer-unique.txt"
    
    # remove dashes 
    echo "Replacing dashes with underscores in gene names..."
    cat "${tissue}-gtex-tcga-integer-unique.txt" | sed 's/-/_/g' > "${tissue}-gtex-tcga-clean.txt"
    
    # make the gem a csv
    echo "Converting to CSV format..."
    cat "${tissue}-gtex-tcga-clean.txt" | sed 's/\s/,/g' > "${tissue}-gtex-tcga-clean.csv"
    
    # create a deseq2 group definition and comparison file
    echo "Creating comparison file for DESeq2..."
    
    # convert dashes to underscores in labels
    cat "$labels_file" | sed 's/-/_/g' > "${tissue}-gtex-tcga-dash.labels.txt"
    
    # add comparison label to each row
    cat "${tissue}-gtex-tcga-dash.labels.txt" | sed 's/\s/,/g' | sed "s/$/,GTEX_${tissue_upper}_TCGA/" > "${tissue}-gtex-tcga.comparison.tmp"
    
    # make header row
    cat "${tissue}-gtex-tcga.comparison.tmp" | sed 's/sample,label,GTEX_'"${tissue_upper}"'_TCGA/Sample,Group,Comparison/' > "${tissue}-gtex-tcga.comparison.csv"
    
    # Clean up temporary files
    rm "${tissue}-gtex-tcga.comparison.tmp"
    
    echo "✅ Processing complete!"
    echo "Output files:"
    echo "- ${tissue}-gtex-tcga-integer.txt"
    echo "- ${tissue}-gtex-tcga-integer-unique.txt"
    echo "- ${tissue}-gtex-tcga-clean.txt"
    echo "- ${tissue}-gtex-tcga-clean.csv"
    echo "- ${tissue}-gtex-tcga-dash.labels.txt"
    echo "- ${tissue}-gtex-tcga.comparison.csv"
}

# Call the function with command line argument or prompt for input
if [ $# -eq 0 ]; then
    echo -n "Please enter tissue type (e.g., bladder, cervix): "
    read tissue_type
else
    tissue_type="$1"
fi

process_tissue_data "$tissue_type"