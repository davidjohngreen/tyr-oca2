#!/bin/bash

# Check if one argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

input_file=$1

# Read the file line by line
while read -r line; do
    # Get the relevent VCF file field from the line
    file1=$(echo $line | awk '{print $5}')

    # Get the third field from the line (gene name), add the 'variant_files/' prefix and '_variants.txt' suffix
    file3=$(echo $line | awk '{print $4}')
    file3="${file3}_regions.bed"

    # Submit the job via bsub
    bsub -q inter -P re_gecip_hearing_and_sight -o output -e error -R "rusage[mem=1GB] span[hosts=1]" ./pipeline.sh $file3 $file1
done < "$input_file"
