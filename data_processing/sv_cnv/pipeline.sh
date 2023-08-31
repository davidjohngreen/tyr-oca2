#!/bin/bash

module load bio/BCFtools/1.9-foss-2018b
module load bio/VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load bio/PLINK/1.9b_4.1-x86_64

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>"
    exit 1
fi

# Assign input arguments to variables
file1=$1
file2=$2


# Extract the prefix from file1 (the file containing the regions) (e.g., "XXXX" from "XXXX_regions.bed")
prefix="${file1%_regions.bed}"


# Define derived filenames using the prefix
raw_vcf="outputs/${prefix}_raw.vcf"
filtered_vcf="outputs/${prefix}_filtered.vcf"
nodups_vcf="outputs/${prefix}_nodups.vcf"
AB_filtered_vcf="outputs/${prefix}_AB_filtered.vcf"
annotated_csv="outputs/${prefix}_annotated.csv"
nodups_filtered_vcf="outputs/${prefix}_nodups_filtered.csv"
transpose_out="outputs/${prefix}_plink_out"

# Define auxillary files for the transposition stage
snps_to_extract="variant_files/${prefix}_snps.txt"


bcftools view -f PASS -R variant_files/$file1 $file2 -o $raw_vcf -Ov
python script_filter_allele_balance.py $raw_vcf $AB_filtered_vcf
bcftools annotate -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' -Ov -o $annotated_csv $AB_filtered_vcf
python script_remove_dups.py $annotated_csv $nodups_vcf
vcftools --vcf $nodups_vcf --minQ 20 --minDP 10 --recode --out $nodups_filtered_vcf
python script_transpose_vcf.py $nodups_filtered_vcf.recode.vcf $snps_to_extract $transpose_out


# Remove intermediate files
rm $raw_vcf $AB_filtered_vcf $nodups_vcf $annotated_csv $nodups_filtered_vcf.recode.vcf
