# recode the haplotypes to count the numeber of deleterious alleles
# for the promotor (modifier) variant, we count REFs (as C is the lower expressing allele) and for the two missense changes we count ALTs


file = open('status-by-ID.txt')
output = open('recoded-status-by-ID.txt', 'w+')
output.write('LPid\tmodifier\tsnp_192\tsnp_402\ttotal\n')


for line in file:
    if line.startswith('LPid'):
        continue
    items = line.strip().split('\t')

    modifier = items[3].split('/')
    modifier_count = modifier.count('0')

    snp_402 = items[1].split('/')
    snp_402_count = snp_402.count('1')

    snp_192 = items[2].split('/')
    snp_192_count = snp_192.count('1')

    total = modifier_count + snp_192_count + snp_402_count

    output.write(line.split('\t')[0] + '\t' + str(modifier_count) + '\t' + str(snp_192_count) + '\t' \
    + str(snp_402_count) + '\t' + str(total) + '\n')
