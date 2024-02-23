library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled by biological gender (female <-> male)
## 3. Remove outliers (Replace outlier signal with marginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by median of each sample 
## 6. Removal of common germline CNVs


## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s

## Options
zero.row.thresh <- 0.4 # Threshold value specifed as double within [0,1].  
#       Warning is issued for markers for which the corresponding row in data frame contains
#       more than zero.row.thresh x (number of samples) zeros. Default value for exome: 0.4
zero.col.thresh <- 0.05 # Threshold value specified as double within [0,1].
#       Samples whose corresponding column contains more than zero.col.thresh x (number of markers) zeros are removed.
#       Default value for exome: 0.05.


load(file=here('05_CCLE_data_preparation/data', 'DOC_CCLE474CellLines.RData')) ## object: DOC

exome.capture.kit.hg38 <- read.delim(file=here('05_CCLE_data_preparation/output/01_lifgOver_Hg19_to_GRCh38', 'exomeAgilentCaptureKit_GRCh38_data.txt')) %>%
  mutate(Probe=paste0('Probe', 1:n())) %>%
  mutate(Start=Start + 1)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

Mahmoud.supp.file <- here('05_CCLE_data_preparation/data/CCLE_Mahmoud2019Nature', '41586_2019_1186_MOESM4_ESM.xlsx')
annotations <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Cell Line Annotations')
datasets <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Datasets')

sample.with.wes <- datasets %>%
  filter(WES_CCLE=='TRUE')

sample.with.wes.depmapid <- annotations %>%
  filter(CCLE_ID %in% sample.with.wes$CCLE_ID) %>%
  pull(depMapID) %>%
  sub('-', '.', .)

## Make a depthOfCoverage with only samples deposited to SRA (PRJNA523380)
doc <- DOC %>%
  unite(col=pos, c('Start', 'Stop'), sep='-') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus') %>%
  select(sample.with.wes.depmapid)

pos.used <- exome.capture.kit.hg38 %>%
  filter(Masked==FALSE) %>%
  mutate(locus=paste(Chr, paste(Start, Stop, sep='-'), sep=':'))

dat <- doc[pos.used$locus, ]

nas <- sapply(dat[!grepl('chrX|chrY', rownames(dat)),], function(x) sum(length(which(is.na(x)))))
nas %>% range() ## 0 0

## There are no NAs in DOC data except Masked==TRUE loci in chrY PAR that have been already removed. Treat 0 as NA instead.
## 1-1. Samples (columns)
count.zeros <- function(list) {
  n.zero <- length(list[list==0])
  return(n.zero)
}

number.of.zeros <- lapply(dat[!grepl('chrX|chrY', rownames(dat)), ], count.zeros)
qc.df <- data.frame(num.zero=unlist(number.of.zeros)) %>%
  rownames_to_column('ModelID') %>%
  mutate(zero.ratio=num.zero/nrow(dat[!grepl('chrX|chrY', rownames(dat)), ])) %>%
  mutate(too.many.zeros=ifelse(zero.ratio >= zero.col.thresh, TRUE, FALSE)) %>%
  left_join(sif %>% select(ModelID, CellLineName, StrippedCellLineName, Sex, OncotreePrimaryDisease, OncotreeLineage), by='ModelID')
saveRDS(qc.df, file=here('05_CCLE_data_preparation/output/03_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'qc.df.rds'), compress=FALSE)

## There were no bad samples
dat.sample.flt <- dat

## 1-2. Probes (rows)
probes.with.bad.marker <- pos.used %>%
  mutate(zero.num=length(dat.sample.flt) - apply(dat.sample.flt, 1, Matrix::nnzero)) %>%
  mutate(zero.ratio=zero.num / length(dat.sample.flt)) %>%
  mutate(bad.marker=case_when(zero.ratio >= zero.row.thresh & Chr!='chrY' ~ TRUE, TRUE ~ FALSE)) # Basically ChrY should be 0 in female samples, so bad markers on ChrY cannot be determined in this way.

dat.sample.flt.male.chrY <- dat.sample.flt[grepl('chrY', rownames(dat.sample.flt)), qc.df %>% filter(ModelID %in% colnames(dat.sample.flt)) %>% filter(Sex=='Male') %>% pull(ModelID)]

probes.with.bad.marker.chrY <- pos.used %>%
  filter(Chr=='chrY') %>%
  mutate(zero.num=length(dat.sample.flt.male.chrY) - apply(dat.sample.flt.male.chrY, 1, Matrix::nnzero)) %>%
  mutate(zero.ratio=zero.num / length(dat.sample.flt)) %>%
  mutate(bad.marker=case_when(zero.ratio >= zero.row.thresh ~ TRUE, TRUE ~ FALSE))

dat.sample.probe.flt <- dat.sample.flt[probes.with.bad.marker %>% filter(bad.marker==FALSE) %>% pull(locus), ]
saveRDS(dat.sample.probe.flt, file=here('05_CCLE_data_preparation/output/03_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'dat.sample.probe.flt.rds'), compress=FALSE)

