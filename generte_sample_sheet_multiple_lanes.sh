#!/bin/bash

# Usage: ./generate_sample_sheet_multiple_lanes.sh --fastq-dir /path/to/fastq_files

# parse command line arguments for the FASTQ directory path
FASTQ_DIR=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fastq-dir=*) FASTQ_DIR="${1#*=}"; shift ;;
        --fastq-dir) FASTQ_DIR="$2"; shift; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
done

if [ -z "$FASTQ_DIR" ]; then
    echo "FASTQ directory path not provided. Use --fastq-dir to specify the directory."
    exit 1
fi

echo "Processing FASTQ files in directory: $FASTQ_DIR"
cd "$FASTQ_DIR" || exit # navigate to the FASTQ directory

temp_file=$(mktemp) # create a temp file to store intermediate results

find "$FASTQ_DIR" -name "*.fastq.gz" | while read -r file_name; do
    absolute_path=$(realpath "$file_name")
    
    # extract sample name, lane, and read number
    sample=$(basename "$file_name" | sed -e 's/\(.*\)_L00[0-9]*_R[12]_001.fastq.gz/\1/')
    lane=$(basename "$file_name" | sed -e 's/.*_L00\([0-9]*\)_R[12]_001.fastq.gz/\1/')
    read_number=$(basename "$file_name" | sed -e 's/.*_R\([12]\)_001.fastq.gz/\1/')
    
    echo "$sample,$lane,$read_number,$absolute_path" >> "$temp_file"
done

# sort and process the temp file to group by sample and lane, then generate the sample sheet
awk -F, '{print $1","$2","$4}' "$temp_file" | sort | awk -F, '{
    if (NR==1 || $1!=prev_sample || $2!=prev_lane) {
        if (NR>1) print prev_line;
        prev_sample=$1; prev_lane=$2; r1=$3; r2="";
    } else {
        r2=$3;
    }
    prev_line=prev_sample","r1","r2;
}
END {
    print prev_line;
}' | awk 'BEGIN {print "sample,fastq_1,fastq_2"} {print}' > sample_sheet.csv

rm "$temp_file"

echo "Sample sheet generated: ${FASTQ_DIR}sample_sheet.csv"
