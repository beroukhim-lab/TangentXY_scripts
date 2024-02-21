library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled by biological gender (female <-> male)
## 3. Remove outliers (Replace outlier signal with merginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by median of each sample 
## 6. Removal of common germline CNVs


## 6. Removal of common germline CNVs
dat.scaled <- readRDS(file=here('05_CCLE_data_preparation/output/03_5_DOC_Preprocessing_scalingByMedian', 'dat.scaled.rds'))

hg38.interval <- read.delim(here('05_CCLE_data_preparation/output/01_lifgOver_Hg19_to_GRCh38/exomeAgilentCaptureKit_GRCh38_data.txt'))

probes.to.be.used <- hg38.interval %>%
  mutate(Start=Start+1) %>%
  mutate(locus=paste(Chr, paste(Start, Stop, sep='-'), sep=':')) %>%
  filter(locus %in% rownames(dat.scaled)) %>%
  filter(freqcnv==FALSE)

dat.gcnv.flt <- dat.scaled[probes.to.be.used$locus, ]

saveRDS(dat.gcnv.flt, file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'dat.gcnv.flt.RData'), compress=FALSE)

## Separate non-cancerous cell lines and cancer cell lines for Tangent
sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))

normal.samples <- sif %>%
  filter(DepMap_ID %in% colnames(dat.gcnv.flt)) %>%
  filter(primary_disease=='Non-Cancerous') %>%
  filter(!is.na(sex)) %>%
  pull(DepMap_ID)

dat.for.analysis.n <- dat.gcnv.flt %>%
  as.data.frame() %>%
  select(normal.samples)

dat.for.analysis.t <- dat.gcnv.flt %>%
  as.data.frame() %>%
  select(!normal.samples)

saveRDS(dat.for.analysis.n, file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_N_QCed_commonCNVremoved.rds'), compress=FALSE)
saveRDS(dat.for.analysis.t, file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_T_QCed_commonCNVremoved.rds'), compress=FALSE)
