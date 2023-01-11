# need to filter the haplotype files to only include specific individuals
# for haplotype analysis, we only want those with unambiguous individuals
# for the analysis of common and rare individuals, we only want those with complete genotypes at all three sites

file = open('status-by-ID.txt', 'r')
out1 = open('status-by-ID-full-genotypes.txt', 'w+')
out2 = open('status-by-ID-full-unabiguous_haplotypes.txt', 'w+')

for line in file:
    if line.startswith('LPid'):
        continue
    items = line.strip().split('\t')[1:]

    if items[0] != './.' and items[1] != './.' and items[2] != './.':
        out1.write(line)

    num_missing = items.count('./.')
    num_het = items.count('0/1')

    if num_missing > 0 and num_het < 1:
        continue
    else:
        out2.write(line)
