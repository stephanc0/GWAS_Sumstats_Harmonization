#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Get the input file from the argument
input_file=$1

# Extract the base name (without extensions)
base_name=$(basename "$input_file" .vcf.gz)

# Generate intermediate and output file names
bed_file="${base_name}.bed"
output_bed="${base_name}_hg38.bed"
unmapped_bed="${base_name}_unmapped.bed"

# Create the BED file
zcat "$input_file" | tail -n +2 | awk '{print "chr"$1"\t"$2-1"\t"$2"\t"$1"\t"$2}' > "$bed_file"

# Run liftOver
./liftOver "$bed_file" hg19ToHg38.over.chain.gz "$output_bed" "$unmapped_bed"

# Inform the user
echo "LiftOver complete!"
echo "Output BED: $output_bed"
echo "Unmapped BED: $unmapped_bed"