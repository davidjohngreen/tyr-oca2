# set a dictionary to store the No. of variant alleles for each individual for each gene
patientDict = {}
out = open('rare_alleles_per_individual.txt', 'w+')

# get the genotype from format 0/0, 0/1, or 1/1
def genotype_generator(n):
    total = 0
    alleles = n.split('/')
    for i in alleles:
        total += int(i)
    return total

# list of genes linking to output of status-by-ID.pl script, which parses VCF to collect status for each individual for a list of variants
genes = open('genes.txt', 'r')
for gene in genes:
    geneName = gene.strip().split()[0]
    header = []
    with open('output/{GENE}_status-by-ID.txt'.format(GENE=geneName), 'r') as f:
        first_line = f.readline().strip().split('\t')
        header = first_line

    with open('output/{GENE}_status-by-ID.txt'.format(GENE=geneName), 'r') as file:
        for line in file:
            if line.startswith('LPid'):
                continue
            indv_total = 0
            indv_variants = {}

            items = line.strip().split('\t')
            for i in range(1,len(items)):
                variant_alleles = genotype_generator(items[i])
                if variant_alleles > 0:
                    indv_total += variant_alleles
                    indv_variants[header[i]] = variant_alleles
            
            # if sample ID already in the dictionary, add the variants for this gene to that individual's dictionary entry, else create an entry 
            # also adds a column with which variants were found, for interest             
            if line.strip().split('\t')[0] in patientDict:
                patientDict[line.strip().split('\t')[0]][0] += indv_total
                patientDict[line.strip().split('\t')[0]][1].update(indv_variants)
            else:
                patientDict[line.strip().split('\t')[0]] = [0, {}]
                patientDict[line.strip().split('\t')[0]][0] = indv_total
                patientDict[line.strip().split('\t')[0]][1] = indv_variants

# print the sample ID followed by the number of variants and which ones they are
for patient in patientDict:
    print(patient, patientDict[patient])
    out.write(patient + '\t' + str(patientDict[patient][0]) + '\t' + ' '.join(patientDict[patient][1]) + '\n')
