import os
import sys

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python3 script.py <file1> <file2>")
    sys.exit(1)

# Assign input arguments to variables
file1 = sys.argv[1]
file2 = sys.argv[2]


file = open(file1,'r')
output = open(file2, 'w+')

for line in file:
    if line.startswith('##'):
        output.write(line)
    if line.startswith('#CHROM'):
        output.write(line)
    if line.startswith('chr'):
        info = line.strip().split('\t')[:9]
        genotypes = line.strip().split('\t')[9:]
        for line in genotypes:
            items = line.strip().split(':')
            genotype = items[0]
            reads = int(items[4])
            reads_ref = int(items[6].split(',')[0])
            reads_alt = int(items[6].split(',')[1])
            '''filter'''
            if genotype == '0/0':
                if float(reads_alt)/reads > 0.1:
                    items[0] = './.'
            elif genotype == '1/1':
                if float(reads_alt)/reads < 0.9:
                    items[0] = './.'
            elif genotype == '0/1' or genotype == '1/0':
                minor_count = min(reads_alt,reads_ref)
                if float(minor_count)/reads < 0.2:
                    items[0] = './.'
            else:
                items[0] = items[0]
            info.append(":".join(items))
        output.write("\t".join(info))
        output.write('\n')
