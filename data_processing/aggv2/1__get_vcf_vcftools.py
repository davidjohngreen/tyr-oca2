# extracts variants from a vcf file and limits to those that have the "PASS" flag
# the "region" .bed files refer to HGMD coordinates for variants, available through the paid professional version of HGMD

import os

genes = open('genes.txt', 'r')

for gene in genes:
    items = gene.strip().split()
    geneName = items[0]
    region = '{GENE}.bed'.format(GENE=geneName)
    file = items[1]
    out = '{GENE}_vcf.vcf'.format(GENE=geneName)
    command = 'bcftools view -f PASS -R {REGION} {FILE} -o {OUT}'.format(REGION=region,FILE=file,OUT=out,GENE=geneName)
    print(command)
    os.system(command)

"""run vcftools to filter on minDP and minGQ"""
import os
genes = open('genes.txt', 'r')

for line in genes:
    gene_name = line.strip().split(' ')[0]

    command = 'vcftools --vcf {GENE}_vcf.vcf --minGQ 20 --minDP 10 --recode --out {GENE}_filtered_vcf.vcf'.format(GENE=gene_name)

    print(command)
    os.system(command)
