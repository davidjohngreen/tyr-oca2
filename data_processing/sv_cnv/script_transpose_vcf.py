import sys

def parse_vcf(vcf_file, variants_file, output_file):
    individuals = []
    variant_names = []
    alt_allele_counts = []

    # Read the variants to include from the file
    with open(variants_file, 'r') as variants:
        variants_to_include = set(variant.strip() for variant in variants)

    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('##'):  # Skip header lines starting with '##'
                continue
            elif line.startswith('#CHROM'):  # Header line with sample names
                individuals = line.strip().split('\t')[9:]
                alt_allele_counts = [[] for _ in individuals]
            elif not line.startswith('#'):  # Variant lines
                fields = line.strip().split('\t')
                variant_name = fields[0] + '_' + fields[1] + '_' + fields[3] + '_' + fields[4]  # Variant name based on CHROM and POS fields
                if variant_name in variants_to_include:

                    variant_names.append(variant_name)
                    genotypes = fields[9:]  # Genotypes start from column 10
                    for i, genotype in enumerate(genotypes):
                        alleles = genotype.split(':')[0]
                        if alleles == './.':
                            alt_allele_counts[i].append('NA')
                        else:
                            alt_count = alleles.count('1')
                            alt_allele_counts[i].append(str(alt_count))

    # Write the transposed counts to the output file
    with open(output_file, 'w') as outfile:
        outfile.write('Sample,' + ','.join(variant_names) + '\n')
        for i, individual in enumerate(individuals):
            outfile.write(individual + ',' + ','.join(alt_allele_counts[i]) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python vcf_to_allele_counts.py <input_vcf> <variants_file> <output_file>")
    else:
        vcf_file = sys.argv[1]
        variants_file = sys.argv[2]
        output_file = sys.argv[3]
        parse_vcf(vcf_file, variants_file, output_file)
