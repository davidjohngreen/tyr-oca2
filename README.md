# tyr-oca2
scripts used to analyse variant combinations associated with ocular albinism

scripts in the folder "data_processing" perform simple bioinformatics manipulations within the Genomics England Research Environment to extract relevent genotypic information and apply standard QC filters to the aggv2 dataset (the production of this dataset is explained at: https://re-docs.genomicsengland.co.uk/aggv2/)

the script "analysis.r" brings together the extracted data to perform simple case-control analysis using the R package "logistf"

the script "ukbiobank_analysis.r" performs simple statistics on phenotypic data extracted from the UK Biobank Research Analysis Platform
