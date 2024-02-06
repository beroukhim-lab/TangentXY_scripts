library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled (female <-> male)
## 3. Remove outliers (Replace outlier signal with merginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by median of each sample 
## 6. Removal of common germline CNVs


## 6. Removal of common germline CNVs
dat.scaled <- readRDS(file=here('02_TCGA_data_preparation/tmp/02_5_DOC_Preprocessing_scaleByMedian', 'dat.scaled.rds'))

probes <- dat.scaled %>%
  rownames() %>%
  as.data.frame() %>%
  setNames('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':', remove=FALSE) %>%
  separate(col=pos, into=c('Start', 'Stop'), sep='-') %>%
  mutate(Start=as.numeric(Start), Stop=as.numeric(Stop)) %>%
  mutate(probe=paste0('probe', 1:n()))

common.germ.cnv.file <- here('02_TCGA_data_preparation/data', 'CNV.hg19.bypos.111213.txt')
common.germ.cnv <- read_delim(common.germ.cnv.file) %>%
  mutate(Chromosome=as.character(Chromosome)) %>%
  mutate(Chromosome=case_when(Chromosome==23 ~ 'X', TRUE ~ Chromosome)) %>%
  mutate(diff=End-Start, diff2=Flankingend-Flankingstart)

probes.with.flag <- probes %>%
  left_join(common.germ.cnv[, c('Chromosome', 'Flankingstart', 'Flankingend')], by=c('Chr'='Chromosome')) %>%
  mutate(overlap=case_when(Stop < Flankingstart | Start > Flankingend ~ FALSE, TRUE ~ TRUE)) %>%
  mutate(overlap=case_when(Chr=='Y' ~ FALSE, TRUE ~ overlap)) %>%
  select(-c('Flankingstart', 'Flankingend')) %>%
  distinct() %>%
  mutate(overlap=factor(.$overlap, levels=c(TRUE, FALSE))) %>%
  arrange(probe, overlap) %>%
  filter(!duplicated(probe)) %>%
  mutate(Chr=factor(.$Chr, levels=c(.$Chr %>% unique() %>% gtools::mixedsort()))) %>%
  arrange(Chr, Start)

probes.to.be.used <- probes.with.flag %>%
  filter(overlap==FALSE) %>%
  select(-overlap)

dat.gcnv.flt <- dat.scaled[probes.to.be.used$locus, ]

saveRDS(dat.gcnv.flt, file=here('02_TCGA_data_preparation/tmp/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'dat.gcnv.flt.RData'), compress=FALSE)

## Separate normal samples and tumor samples for downstream analyses
dat.for.analysis.n <- dat.gcnv.flt %>%
  as.data.frame() %>%
  select(ends_with('NB') | ends_with('NBC') | ends_with('NT'))
dat.for.analysis.t <- dat.gcnv.flt %>%
  as.data.frame() %>%
  select(ends_with('TAP') | ends_with('TB') | ends_with('TM') | ends_with('TP'))

saveRDS(dat.for.analysis.n, file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'), compress=FALSE)
saveRDS(dat.for.analysis.t, file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'), compress=FALSE)
