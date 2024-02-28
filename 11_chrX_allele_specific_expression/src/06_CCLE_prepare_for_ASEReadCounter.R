library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

## Samples with SNPs call
snp.dir <- here('11_chrX_allele_specific_expression/data/CCLE_SNPs_call_on_WES')

snp.files <- list.files(snp.dir, recursive=TRUE, full.names=TRUE)

snp.files.df <- snp.files %>%
  as.data.frame() %>%
  setNames('path') %>%
  mutate(file.name=basename(path)) %>%
  mutate(ModelID=basename(dirname(path))) %>%
  mutate(ModelID=sub('-', '.', ModelID)) %>%
  mutate(pass=case_when(grepl('^PASS', file.name) ~ TRUE, TRUE ~ FALSE))
write.table(snp.files.df, file=here('11_chrX_allele_specific_expression/output/06_CCLE_prepare_for_ASEReadCounter', 'WES_gSNPs_files.txt'), sep='\t', row.names=FALSE, quote=FALSE, na='')

terra.data.list <- read.delim(here('11_chrX_allele_specific_expression/data/CCLE_Terra_output', 'sample.tsv'))

snp.google.backet.url <- snp.files.df %>%
  filter(pass==TRUE) %>%
  filter(ModelID %in% sub('-', '.', terra.data.list$entity.sample_id)) %>%
  select(ModelID, file.name) %>%
  mutate(file.type=case_when(grepl('tbi$', file.name) ~ 'WES_gSNPs_VCF_index', TRUE ~ 'WES_gSNPs_VCF')) %>%
  mutate(file.name=paste0('gs://fc-secure-516f8144-b7c9-458a-934f-2cc78f09cfb3/', file.name)) %>%
  pivot_wider(names_from=file.type, values_from=file.name) %>%
  mutate(ModelID=sub('\\.', '-', ModelID))

terra.data.list.new <- terra.data.list %>%
  left_join(snp.google.backet.url, by=c('entity.sample_id'='ModelID')) %>%
  mutate(WES_gSNPs_VCF=case_when(is.na(WES_gSNPs_VCF) ~ '', TRUE ~ WES_gSNPs_VCF)) %>%
  mutate(WES_gSNPs_VCF_index=case_when(is.na(WES_gSNPs_VCF_index) ~ '', TRUE ~ WES_gSNPs_VCF_index)) %>%
  rename(`entity:sample_id`=entity.sample_id)
write.table(terra.data.list.new, file=here('11_chrX_allele_specific_expression/output/06_CCLE_prepare_for_ASEReadCounter', 'sample_table_for_Terra.tsv'), sep='\t', row.names=FALSE, quote=FALSE, na='')

snp.pass.vcf.files.to.be.uploaded <- snp.files.df %>%
  filter(pass==TRUE) %>%
  filter(ModelID %in% sub('-', '.', terra.data.list$entity.sample_id)) %>%
  select(path)
write.table(snp.pass.vcf.files.to.be.uploaded, file=here('11_chrX_allele_specific_expression/output/06_CCLE_prepare_for_ASEReadCounter', 'gSNPs_PASS_VCFs_toBeUploadedToGoogleBucket.txt'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE, na='')
