# with genotype data in column format (one column per variant, one row per individual), runs through and assigns homozygous haplotype
# for example, an individual with genotypes 0/0, 1/1, and 1/1 for the promotor variant (C>T) and two missense changes (G>A and C>A) has a haplotype CAA

file = open('status-by-ID.txt', 'r')
output = open('haplotypes.txt', 'w+')
probandFile = open('proband-IDs.txt', 'r')
probands = []

for i in probandFile:
    probands.append(i.strip('\n'))

# refers to indices of the columns
var_402 = 0
var_192 = 1
modifier = 2

ref = 0
hmz = 2


IDs = {'CAG': [], 'CAA': [], 'TAG': [], 'TAA': [], \
    'TCA': [], 'CCA': [], 'TCG': [], 'CCG': [], '': []}


def genotype_generator(n):
    if n != './.':
        genotype = [int(i) for i in n.strip().split('/')]
        return sum(genotype)

# skip header, then look at columns 2, 3, and 4 (with indices 0, 1, and 2 as above after making "items" from the line without the first column (LpID)
for line in file:
    if line.startswith('LPid'):
        continue
    if line.strip().split('\t')[0] in probands:
        items = line.strip().split('\t')[1:]

        """check haplotypes"""
        if genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == ref:
            IDs['CAG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == hmz:
            IDs['CAA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == ref:
            IDs['TAG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == hmz and \
        genotype_generator(items[var_402]) == hmz:
            IDs['TAA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == hmz:
            IDs['TCA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == hmz:
            IDs['CCA'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == hmz and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == ref:
            IDs['TCG'].append(line.strip().split('\t')[0])

        elif genotype_generator(items[modifier]) == ref and \
        genotype_generator(items[var_192]) == ref and \
        genotype_generator(items[var_402]) == ref:
            IDs['CCG'].append(line.strip().split('\t')[0])
        else:
            IDs[''].append(line.strip().split('\t')[0])

# for each sample ID, print out the ID plus the homozygous haplotype            
for i in IDs:
    for j in IDs[i]:
        output.write(i + '\t' + j + '\n')
