library(tidyverse)
library(here)

terra.sample.file <- here('11_chrX_allele_specific_expression/data/Terra_output', 'sample.tsv')
terra.sample <- read.delim(terra.sample.file, na.strings=c('', 'NA')) %>%
  select(c('entity.sample_id', 'tcga_sample_id', 'participant', 'sample_type', 'mRNASeq_bam_analysis_id', 'mRNASeq_bai_path', 'mRNASeq_bam_path')) %>%
  mutate(tcga_participant_id=sub('-[0-9]{2}$', '', tcga_sample_id)) %>%
  separate(col=participant, into=c('project', 'tss', 'id'), sep='-', remove=FALSE) %>%
  select(-c('tss', 'id')) %>%
  select(entity.sample_id, project, tcga_sample_id, participant, tcga_participant_id, everything())

vcf.samples <- read.delim(here('11_chrX_allele_specific_expression/output/01_TCGA_prepare_vcf_for_ASEReadCounter', 'sample.list.txt'), header=FALSE) %>%
  setNames('vcf_barcode') %>%
  mutate(vcf_tcga_sample_id=str_extract(vcf_barcode, '^TCGA-.{2}-.{4}-.{2}')) %>%
  mutate(tcga_participant_id=sub('-[0-9]{2}$', '', vcf_tcga_sample_id)) %>%
  mutate(vcf_sample_type_code=str_extract(vcf_tcga_sample_id, '[0-9]{2}$')) %>%
  mutate(vcf_sample_type=case_when(vcf_sample_type_code=='01' ~ 'TP',
                                    vcf_sample_type_code=='10' ~ 'NB',
                                    vcf_sample_type_code=='11' ~ 'NT',
                                    vcf_sample_type_code=='12' ~ 'NBC',
                                    vcf_sample_type_code=='14' ~ 'NBM'))

terra.table <- terra.sample %>%
  filter(!is.na(mRNASeq_bam_path)) %>%
  inner_join(vcf.samples, by='tcga_participant_id') %>%
  rename(`entity:sample_id`=entity.sample_id) %>%
  rename(`pfb:mRNASeq_bam_analysis_id`=mRNASeq_bam_analysis_id) %>%
  rename(`pfb:mRNASeq_bai_path`=mRNASeq_bai_path) %>%
  rename(`pfb:mRNASeq_bam_path`=mRNASeq_bam_path)
write.table(terra.table, file=here('11_chrX_allele_specific_expression/output/02_TCGA_make_input_table_for_Terra', 'sample_table_for_Terra.tsv'), sep='\t', row.names=FALSE, quote=FALSE, na='')

sample.set.table <- terra.table %>%
  select(project, `entity:sample_id`) %>%
  rename(`membership:sample_set_id`=project) %>%
  rename(sample=`entity:sample_id`)
write.table(sample.set.table, file=here('11_chrX_allele_specific_expression/output/02_TCGA_make_input_table_for_Terra', 'sample_set_table_for_Terra.tsv'), sep='\t', row.names=FALSE, quote=FALSE, na='')
