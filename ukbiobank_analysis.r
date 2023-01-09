# to extract a dataset from RAP you need to run the following:
# dx extract_dataset record-GK2J5YQJyQ6bkzJyGZjJ17kQ --fields "participant.eid,participant.p22001" -o output_name.csv
# note that non-europeans were excluded at the cohort building stage in RAP

library(dplyr)
library(ggpubr)
library(logistf)


participant_table <- read.csv(file = 'uk_biobank/phenotype_data.csv')
albinism_IDs <- read.csv(file = 'uk_biobank/albinism_IDs.csv')
snp_402 <- read.csv(file = 'TYR_402_results.txt', sep=' ')
snp_443 <- read.csv(file = 'OCA2_443_results.txt', sep='\t')
sex <- read.csv(file = 'cohort_with_sex.csv', sep=',')

# invert the direction of the 443 snp
snp_443 <- snp_443 %>% mutate(rs121918166_C_corrected = case_when(rs121918166_C == 0 ~ 2,
																  rs121918166_C == 1 ~ 1,
																  rs121918166_C == 2 ~ 0))

# inner join to get both snps together
snps <- inner_join(snp_443, snp_402, by='participant_id') %>% 
		select(c('participant_id','rs121918166_C_corrected','rs1126809_A'))

snps <- inner_join(snps, sex, by='participant_id')

# join the snps and phenotype data, mutate to calculate the mean of logmar and retinal thickness.
pheno.geno <- inner_join(snps, participant_table, by='participant_id') %>%
		rowwise() %>% 
		mutate(logmar_mean = mean(c(logmar_left,logmar_right))) %>%
		mutate(central_retinal_thickness_mean = mean(c(central_retinal_thickness_left,central_retinal_thickness_right))) %>%
		mutate(status = case_when(rs121918166_C_corrected == 0 & rs1126809_A == 0 ~ 'ref_ref',
								  rs121918166_C_corrected == 0 & rs1126809_A == 1 ~ 'one_het',
								  rs121918166_C_corrected == 1 & rs1126809_A == 0 ~ 'one_het',
								  rs121918166_C_corrected == 1 & rs1126809_A == 1 ~ 'het_het',
								  rs121918166_C_corrected == 2 & rs1126809_A == 1 ~ 'hetzandhomz',
								  rs121918166_C_corrected == 1 & rs1126809_A == 2 ~ 'hetzandhomz',
								  rs121918166_C_corrected == 2 & rs1126809_A == 2 ~ 'hetzandhomz',
								  rs121918166_C_corrected == 0 & rs1126809_A == 2 ~ 'ref443_hom402')) %>%
		mutate(status_vs_other = case_when(rs121918166_C_corrected == 1 & rs1126809_A == 1 ~ 'het_het')) %>% replace_na(list(status_vs_other = 'other'))




# produce summary statistics for the groups, only counting those with a non-NA value for the relevant field
# MEAN
group_by(pheno.geno, status) %>%
  summarise(
  	non_na_count = sum(!is.na(logmar_mean)),
    mean = mean(logmar_mean, na.rm = TRUE),
    median = median(logmar_mean, na.rm = TRUE),
    sd = sd(logmar_mean, na.rm = TRUE)
  )

group_by(pheno.geno, status) %>%
  summarise(
  	non_na_count = sum(!is.na(central_retinal_thickness_mean)),
    mean = mean(central_retinal_thickness_mean, na.rm = TRUE),
    median = median(central_retinal_thickness_mean, na.rm = TRUE),
    sd = sd(central_retinal_thickness_mean, na.rm = TRUE)
  )

# LEFT
group_by(pheno.geno, status) %>%
  summarise(
  	non_na_count = sum(!is.na(logmar_left)),
    mean = mean(logmar_left, na.rm = TRUE),
    median = median(logmar_left, na.rm = TRUE),
    sd = sd(logmar_left, na.rm = TRUE)
  )

group_by(pheno.geno, status) %>%
  summarise(
  	non_na_count = sum(!is.na(central_retinal_thickness_left)),
    mean = mean(central_retinal_thickness_left, na.rm = TRUE),
    median = median(central_retinal_thickness_left, na.rm = TRUE),
    sd = sd(central_retinal_thickness_left, na.rm = TRUE)
  )


# kruskal to compare the gorups
kruskal.test(logmar_left ~ status, data = pheno.geno)
kruskal.test(central_retinal_thickness_left ~ status, data = pheno.geno)


kruskal.test(logmar_left ~ status_vs_other, data = pheno.geno)
kruskal.test(central_retinal_thickness_left ~ status_vs_other, data = pheno.geno)


# perform pairwise comparisons of the groups using pairwise wilxocon (non-parametric)
# MEAN
pairwise.wilcox.test(pheno.geno$logmar_mean, pheno.geno$status,
                 	 p.adjust.method = "BH")

pairwise.wilcox.test(pheno.geno$central_retinal_thickness_mean, pheno.geno$status,
                 	 p.adjust.method = "BH")

# LEFT
pairwise.wilcox.test(pheno.geno$logmar_left, pheno.geno$status,
                 	 p.adjust.method = "BH")

pairwise.wilcox.test(pheno.geno$central_retinal_thickness_left, pheno.geno$status,
                 	 p.adjust.method = "BH")




# change order of factor


target <- c('one_het', 'ref_ref', 'het_het', 'hetzandhomz')
pheno.geno.plot <- pheno.geno %>% filter(status %in% target)

pheno.geno.plot$status <- factor(pheno.geno.plot$status, levels = c("hetzandhomz", "het_het", "one_het", "ref_ref"))

levels(pheno.geno.plot$status) <- list(hetzandhomz  = "Group D", het_het = "Group C", one_het = "Group B", ref_ref = "Group A")

# MAKE PLOTS LEFT
jpeg("logmar_left.jpeg", units="cm", width=10, height=10, res=300)
logmar_violin <- ggplot(pheno.geno.plot, aes(x=status, y=logmar_left, fill=status)) + 
				 geom_violin(adjust = 3, trim=FALSE, fill='#FFFFFF') + coord_flip() +
				 guides(fill = guide_legend(reverse = TRUE)) +
				 stat_boxplot(geom ='errorbar', width = 0.1) +
				 geom_boxplot(width=0.1, fill="white",outlier.shape = NA) +
				 scale_x_discrete(labels=c("Group D", "Group C", "Group B", "Group A"))
				 # stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='black', fill='#FFFFFF')
logmar_violin
dev.off()

jpeg("thickness_left.jpeg", units="cm", width=10, height=10, res=300)
thickness_violin <- ggplot(pheno.geno.plot, aes(x=status, y=central_retinal_thickness_left, fill=status)) + 
				 geom_violin(adjust = 3, trim=FALSE, fill='#FFFFFF') + coord_flip() +
				 guides(fill = guide_legend(reverse = TRUE)) +
				 stat_boxplot(geom ='errorbar', width = 0.1) +
				 geom_boxplot(width=0.1, fill="white",outlier.shape = NA) +
				 scale_x_discrete(labels=c("Group D", "Group C", "Group B", "Group A"))
				 # stat_summary(fun.y=mean, geom="point", shape=23, size=2, color='black', fill='#FFFFFF')
thickness_violin
dev.off()









