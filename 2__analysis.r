library(dplyr)
library(tidyverse)
library(logistf)
library(janitor)

#import haplotypes for individuals with known number of common variants
haplotypes = read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/haplotype_groupings_full_genotypes.csv', sep='\t',stringsAsFactors = T, header = T)

#import the aggv2 sample stats file to remove redacted individuals and get the ancestry information and participant ID'''
sample_stats = read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/sample_stats.txt', header = T, stringsAsFactors = T, sep = '\t')
sample_stats = sample_stats %>%
    select(Participant.Id, Platekey, Participant.Phenotypic.Sex,
            Pred.African.Ancestries, Pred.South.Asian.Ancestries,
            Pred.East.Asian.Ancestries, Pred.European.Ancestries,
            Pred.American.Ancestries) %>%
    mutate(ancestry = case_when(Pred.European.Ancestries >= 0.8 ~ "eur",
            Pred.African.Ancestries >= 0.8 ~ "afr",
            Pred.South.Asian.Ancestries >= 0.8 ~ "sas",
            Pred.East.Asian.Ancestries >= 0.8 ~ "eas",
            Pred.American.Ancestries >= 0.8 ~ "amr")) %>%
    select(Participant.Id, Platekey, Participant.Phenotypic.Sex, ancestry) %>%
    rename(LPid = Platekey, participant_id=Participant.Id,sex=Participant.Phenotypic.Sex)

#use an inner join to map the haplotypes file to the sample stats, retaining only those samples found in both datasets'''
combined = inner_join(haplotypes, sample_stats, by='LPid') %>%
          rename(plate_key = LPid)

#import the number of additional HGMD-listed variants for each sample'''
additional_variants = read.csv('/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/additional.csv', header=T, stringsAsFactors=T, sep='\t')

#update the dataset to include the rare variants and mark NA as 0'''
combined = left_join(combined, additional_variants, by='plate_key') %>%
          replace_na(list(additional = 0))

#import the IDs of individuals with albinism, rename partidipant id, add to the dataset using left join and set NA to 0'''
albinism <- read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/albinism_IDs.csv', stringsAsFactors = T, header = T, sep = '\t') %>% 
          rename(participant_id = partid)
combined <- left_join(combined, albinism, by="participant_id") %>% 
            mutate(combined, albinism = ifelse(is.na(albinism), 0, albinism)) %>%
            mutate(combined, cohort = 'GeL')

#import the number of common variants per individual'''
common_count <- read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/recoded-status-by-ID-full-genotypes.txt', stringsAsFactors = T, header = T, sep = '\t') %>% 
            select('plate_key', 'common_TYR')
combined <- left_join(combined, common_count, by='plate_key')

#import the 443 and 402 variant files"""
variant_oca2 <- read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/albinism_digenic/OCA2_recodeA.raw', stringsAsFactors = T, header = T, sep = ' ') %>% 
            select('IID', 'chr15.27985101.C.T_T') %>% 
            rename(LPid = IID, snp_oca2=chr15.27985101.C.T_T)
variant_oca2 <- na.omit(variant_oca2)
variant_tyr = read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/albinism_digenic/402.txt', stringsAsFactors = T, header = T, sep = '\t')

#process the snps and set statuses in a new column"""
snps = inner_join(variant_oca2, variant_tyr, by='LPid') %>% 
            rename(plate_key='LPid')

#add the number of rare additional variants in tyr"""
additionaltyr = read.csv('/home/dgreen/re_gecip/hearing_and_sight/dgreen/TYR/regression_analysis/final_thing/accessory_files/additional_TYR.csv', header=T, stringsAsFactors=T, sep='\t')
combined = left_join(combined, additionaltyr, by='plate_key') %>% 
            replace_na(list(additional_tyr = 0))

combined = inner_join(combined, snps, by='plate_key')


#CREATE THE FULL DATASET BY MERGING BORDEAUX AND GEL AND REMOVE THE INDIVIDUAL WITH INFETERMINATE SEX"""
#change this to match the oca digenic input file"""
cohort <- read.csv(file = '/home/dgreen/re_gecip/hearing_and_sight/dgreen/albinism_digenic/cohort.csv', stringsAsFactors = T, header = T, sep = '\t')
full_dataset <- bind_rows(combined, cohort) %>% select(-c(additional_tyr)) %>% filter(sex != 'Indeterminate')

#clean up environment'''
rm(additional_variants, additionaltyr, common_count, sample_stats, variant_oca2, variant_tyr, snps, cohort, combined, albinism, haplotypes)


#create a factor containing the different statuses to compare for basic analysis"""
full_dataset <- full_dataset %>% mutate(status = case_when(snp_oca2 == 1 & snp_402 == 1 ~ 'het_het')) %>% replace_na(list(status = 'all_other'))
full_dataset <- full_dataset %>% mutate(status_ref_ref = case_when(snp_oca2 == 1 & snp_402 == 1 ~ 'het_het',
                                                                   snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref'))

#assign the new groups for gradient analysis as per panos' idea"""
full_dataset <- full_dataset %>% mutate(status_all_new = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 0 & snp_402 == 1 ~ 'one_het',
                                                           snp_oca2 == 1 & snp_402 == 0 ~ 'one_het',
                                                           snp_oca2 == 1 & snp_402 == 1 ~ 'het_het',
                                                           snp_oca2 == 1 & snp_402 == 2 ~ 'homzandhetz',
                                                           snp_oca2 == 2 & snp_402 == 1 ~ 'homzandhetz',
                                                           snp_oca2 == 2 & snp_402 == 2 ~ 'homzandhetz',
                                                           snp_oca2 == 2 & snp_402 == 0 ~ 'homandref',
                                                           snp_oca2 == 0 & snp_402 == 2 ~ 'refandhom'))
full_dataset$status_all_new <- as.factor(full_dataset$status_all_new)
full_dataset$status_all_new <- relevel(full_dataset$status_all_new, ref='ref_ref')

#assign the different gradient levels'''
#ONE_HET'''
full_dataset <- full_dataset %>% mutate(one_het = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 0 & snp_402 == 1 ~ 'one_het',
                                                           snp_oca2 == 1 & snp_402 == 0 ~ 'one_het'))

full_dataset$one_het <- as.factor(full_dataset$one_het)
full_dataset$one_het <- relevel(full_dataset$one_het, ref = 'ref_ref')

#HET_HET'''
full_dataset <- full_dataset %>% mutate(het_het = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 1 & snp_402 == 1 ~ 'het_het'))


full_dataset$het_het <- as.factor(full_dataset$het_het)
full_dataset$het_het <- relevel(full_dataset$het_het, ref = 'ref_ref')

#HOMZ_AND_HETS'''
full_dataset <- full_dataset %>% mutate(homzandhets = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 1 & snp_402 == 2 ~ 'homzandhetz',
                                                           snp_oca2 == 2 & snp_402 == 1 ~ 'homzandhetz',
                                                           snp_oca2 == 2 & snp_402 == 2 ~ 'homzandhetz'))

full_dataset$homzandhets <- as.factor(full_dataset$homzandhets)
full_dataset$homzandhets <- relevel(full_dataset$homzandhets, ref = 'ref_ref')

#REF_AND_HOM'''
full_dataset <- full_dataset %>% mutate(refandhom = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 0 & snp_402 == 2 ~ 'ref_and_hom'))

full_dataset$refandhom <- as.factor(full_dataset$refandhom)
full_dataset$refandhom <- relevel(full_dataset$refandhom, ref = 'ref_ref')

#HOM_AND_REF'''
full_dataset <- full_dataset %>% mutate(homandref = case_when(snp_oca2 == 0 & snp_402 == 0 ~ 'ref_ref',
                                                           snp_oca2 == 2 & snp_402 == 0 ~ 'hom_and_ref'))

full_dataset$homandref <- as.factor(full_dataset$homandref)
full_dataset$homandref <- relevel(full_dataset$homandref, ref = 'ref_ref')



#get other numbers for the tables'''
full_dataset <- full_dataset %>% mutate(table_numbers = case_when(snp_oca2 == 1 & snp_402 == 0 ~ 'het_ref',
                                                                  snp_oca2 == 0 & snp_402 == 1 ~ 'ref_het',
                                                                  snp_oca2 == 2 & snp_402 == 1 ~ 'hom_het',
                                                                  snp_oca2 == 1 & snp_402 == 2 ~ 'het_hom',
                                                                  snp_oca2 == 2 & snp_402 == 2 ~ 'hom_hom'))


#set columns to correct data types'''
full_dataset$additional <- as.numeric(full_dataset$additional)
full_dataset$ancestry <- as.factor(full_dataset$ancestry)
full_dataset$status <- as.factor(full_dataset$status)
full_dataset$status_ref_ref <- as.factor(full_dataset$status_ref_ref)
full_dataset$common_TYR <- as.numeric(full_dataset$common_TYR)
full_dataset$snp_402 <- as.numeric(full_dataset$snp_402)
full_dataset$snp_oca2 <- as.numeric(full_dataset$snp_oca2)
full_dataset$sex <- as.factor(full_dataset$sex)
full_dataset$status_all_new <- as.factor(full_dataset$status_all_new)


#correct the rare and common variant columns to account for 443 and 402'''
full_dataset <- full_dataset %>% mutate(additional_corrected = case_when(snp_oca2 == 1 ~ additional - snp_oca2,
                                                                         snp_oca2 == 2 ~ additional - snp_oca2,
                                                                         snp_oca2 == 0 ~ additional))

full_dataset <- full_dataset %>% mutate(common_TYR_corrected = case_when(snp_402 == 1 ~ common_TYR - snp_402,
                                                                         snp_402 == 2 ~ common_TYR - snp_402,
                                                                         snp_402 == 0 ~ common_TYR))



#PREPARE DATA FOR ANALYSIS"""
#relevel factors factors and change order"""
full_dataset$status <- relevel(full_dataset$status, ref = 'all_other')
full_dataset$status_ref_ref <- relevel(full_dataset$status_ref_ref, ref = 'ref_ref')
full_dataset$status_all_new <- relevel(full_dataset$status_all_new, ref = 'ref_ref')
full_dataset$ancestry <- relevel(full_dataset$ancestry, ref = "eur")
full_dataset$sex <- toupper(full_dataset$sex)
full_dataset <- full_dataset %>% relocate(additional_corrected, .after = snp_402)
full_dataset <- full_dataset %>% relocate(common_TYR_corrected, .after = additional_corrected)




#filter to split the dataset into the three different datasets - full, european and GeL only"""
full_dataset_european <- full_dataset %>% filter(additional_corrected == 0 & ancestry == 'eur')
gel_dataset <- full_dataset %>% filter(cohort == 'GeL' & sex != 'INDETERMINATE')
gel_european_dataset <- full_dataset_european %>% filter(cohort == 'GeL' & sex != 'INDETERMINATE')




#MAIN OVERALL ORs FOR HET-HET VS. REF-REF AND FOR HET-HET VS. ALL OTHER"""
#VS. All OTHER"""
fit1 <- logistf(data=full_dataset_european, albinism~status+common_TYR_corrected+sex, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
primary_OR_all_other <- data.frame(exp(cbind(OR = coef(fit1),confint(fit1)))) %>% mutate_if(is.numeric, round, digits = 1)

fit2 <- logistf(data=full_dataset, albinism~status+common_TYR_corrected+ancestry+additional_corrected+sex, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
secondary_OR_all_other <- data.frame(exp(cbind(OR = coef(fit2),confint(fit2)))) %>% mutate_if(is.numeric, round, digits = 1)

#VS. REF-REF"""
fit3 <- logistf(data=full_dataset_european, albinism~status_ref_ref+common_TYR_corrected+sex, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
primary_OR_ref_ref <- data.frame(exp(cbind(OR = coef(fit3),confint(fit3)))) %>% mutate_if(is.numeric, round, digits = 1)

fit4 <- logistf(data=full_dataset, albinism~status_ref_ref+common_TYR_corrected+ancestry+additional_corrected+sex, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
secondary_OR_ref_ref <- data.frame(exp(cbind(OR = coef(fit4),confint(fit4)))) %>% mutate_if(is.numeric, round, digits = 1)






#GEL and BORDEAUX PRIMARY ANALYSIS"""
table(full_dataset_european$albinism)

fit5 <- logistf(data=full_dataset_european, albinism~one_het+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_primary_one_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset_european$albinism, full_dataset_european$one_het)

fit5 <- logistf(data=full_dataset_european, albinism~het_het+common_TYR_corrected+sex, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_primary_het_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset_european$albinism, full_dataset_european$het_het)

fit5 <- logistf(data=full_dataset_european, albinism~homzandhets+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_primary_homzandhets <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset_european$albinism, full_dataset_european$homzandhets)

fit5 <- logistf(data=full_dataset_european, albinism~refandhom+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_primary_refandhom <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset_european$albinism, full_dataset_european$refandhom)

fit5 <- logistf(data=full_dataset_european, albinism~homandref+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_primary_homandref <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset_european$albinism, full_dataset_european$hom_ref)



#GEL and BORDEAUX SECONDARY ANALYSIS """
table(full_dataset$albinism)

fit5 <- logistf(data=full_dataset, albinism~one_het+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_secondary_one_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset$albinism, full_dataset$one_het)

fit5 <- logistf(data=full_dataset, albinism~het_het+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_secondary_het_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset$albinism, full_dataset$het_het)

fit5 <- logistf(data=full_dataset, albinism~homzandhets+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_secondary_homzandhets <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset$albinism, full_dataset$homzandhets)

fit5 <- logistf(data=full_dataset, albinism~refandhom+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_secondary_refandhom <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset$albinism, full_dataset$refandhom)

fit5 <- logistf(data=full_dataset, albinism~homandref+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gradient_secondary_homandref <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(full_dataset$albinism, full_dataset$homandref)


#GEL ONLY PRIMARY ANALYSIS'''
table(gel_european_dataset$albinism)

fit5 <- logistf(data=gel_european_dataset, albinism~one_het+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_euro_gradient_secondary_one_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_european_dataset$albinism, gel_european_dataset$one_het)

fit5 <- logistf(data=gel_european_dataset, albinism~het_het+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_euro_gradient_secondary_het_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_european_dataset$albinism, gel_european_dataset$het_het)

fit5 <- logistf(data=gel_european_dataset, albinism~homzandhets+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_euro_gradient_secondary_homzandhets <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_european_dataset$albinism, gel_european_dataset$homzandhets)

fit5 <- logistf(data=gel_european_dataset, albinism~refandhom+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_euro_gradient_secondary_refandhom <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_european_dataset$albinism, gel_european_dataset$refandhom)

fit5 <- logistf(data=gel_european_dataset, albinism~homandref+sex+common_TYR_corrected, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_euro_gradient_secondary_homandref <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_european_dataset$albinism, gel_european_dataset$homandref)



#GEL ONLY ANALYSIS'''
table(gel_dataset$albinism)

fit5 <- logistf(data=gel_dataset, albinism~one_het+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_gradient_secondary_one_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_dataset$albinism, gel_dataset$one_het)

fit5 <- logistf(data=gel_dataset, albinism~het_het+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_gradient_secondary_het_het <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_dataset$albinism, gel_dataset$het_het)

fit5 <- logistf(data=gel_dataset, albinism~homzandhets+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_gradient_secondary_homzandhets <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_dataset$albinism, gel_dataset$homzandhets)

fit5 <- logistf(data=gel_dataset, albinism~refandhom+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_gradient_secondary_refandhom <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_dataset$albinism, gel_dataset$refandhom)

fit5 <- logistf(data=gel_dataset, albinism~homandref+sex+common_TYR_corrected+additional_corrected+ancestry, firth=TRUE, pl=TRUE, control=logistf.control(maxit=100, maxstep=20))
gel_gradient_secondary_homandref <- data.frame(exp(cbind(OR = coef(fit5),confint(fit5)))) %>% mutate_if(is.numeric, round, digits = 1)
table(gel_dataset$albinism, gel_dataset$homandref)



