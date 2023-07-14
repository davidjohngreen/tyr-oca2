import os
genes = open('genes.txt', 'r')

for line in genes:
    gene_name = line.strip().split(' ')[0]

    file = open('{GENE}_filtered_GQ20_DP10_vcf.vcf.recode.vcf'.format(GENE=gene_name), 'r')
    output = open('{GENE}_filtered_GQ20_DP10_vcf.recode.vcf_AB_filtered.vcf'.format(GENE=gene_name), 'w+')

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
