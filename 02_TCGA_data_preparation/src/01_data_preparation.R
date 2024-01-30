# dir <- '/xchip/beroukhimlab/kei/project/Tangent/20220725_Tangent_sex_TCGA_WES'
dir <- '/xchip/beroukhimlab/kei/project/TangentXY_scripts/02_TCGA_data_preparation'
setwd(dir)

library(useful)
library(tidyverse)
library(here)

doc.n <- readRDS(file=here('../../Tangent/20220725_Tangent_sex_TCGA_WES/output/DOC/merged', 'TCGA_WES_hg19_N.RData'))
doc.t <- readRDS(file=here('../../Tangent/20220725_Tangent_sex_TCGA_WES/output/DOC/merged', 'TCGA_WES_hg19_T.RData'))


#Make annotation for probes
pos <- doc.n %>%
  rownames() %>%
  as.data.frame() %>%
  separate(data=., col=., into=c('Chr', 'pos'), sep=':') %>%
  separate(data=., col=pos, into=c('Start', 'Stop'), sep='-') %>%
  mutate(Start=as.numeric(Start), Stop=as.numeric(Stop))

## Whole exome capture kit (developed by Agilent and customized by Broad Institute) can be downloaded from URL below
## https://bitbucket.org/cghub/cghub-capture-kit-info/src/b33c1cb33593/BI/vendor/Agilent/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed?at=master
exome.capture.kit.file <- here('data', 'whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed')
exome.capture.kit.bed <- read.delim(exome.capture.kit.file, header=F) %>%
  setNames(c('Chr', 'Start', 'Stop', 'Target', 'Strand')) %>%
  mutate(Start=Start + 1)

## Add annotation about PAR
par.file <- here('data', 'PAR.txt')
par <- read_delim(par.file)
par.hg19 <- par %>%
  filter(Version=='hg19')

hg19 <- rCGH::hg19 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))

pos.annotated <- pos %>%
  left_join(exome.capture.kit.bed, by=c('Chr', 'Start', 'Stop')) %>%
  mutate(Par=case_when(Chr=='X' & Stop < par.hg19 %>% filter(Chr=='X' & Name=='par1') %>% pull(Stop) ~ 'par1',
                        Chr=='X' & Start > par.hg19 %>% filter(Chr=='X' & Name=='par2') %>% pull(Start) ~ 'par2',
                        Chr=='Y' & Stop < par.hg19 %>% filter(Chr=='Y' & Name=='par1') %>% pull(Stop) ~ 'par1',
                        Chr=='Y' & Start > par.hg19 %>% filter(Chr=='Y' & Name=='par2') %>% pull(Start) ~ 'par2')) %>%
  mutate(Masked=case_when(Chr=='Y' & !is.na(Par) ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(probe=paste0('Probe', 1:n())) %>%
  left_join(hg19, by=c('Chr'='chrom')) %>%
  mutate(arm=case_when(Stop < centromerStart ~ 'p',
                        Start > centromerEnd ~ 'q')) %>%
  mutate(chrom=paste0('chr', Chr))
saveRDS(pos.annotated, file=here('output', 'probes.hg19.annotated.rds'), compress=FALSE)
