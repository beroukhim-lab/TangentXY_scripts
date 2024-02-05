library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled (female <-> male)
## 3. Remove outliers (Replace outlier signal with marginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by mean/median of each sample 
## 6. Removal of common germline CNVs
## 7. Remove mislabeled samples (NB <-> TP)

## Options
floor.val.frac <- 0.001 # Threshold value specified as double within [0,1].
#       CN values will be floored at floor.val.frac x (data mean)
zero.row.thresh <- 0.4 # Threshold value specifed as double within [0,1].  
#       Warning is issued for markers for which the corresponding row in data frame contains
#       more than zero.row.thresh x (number of samples) zeros. Default value for exome: 0.4
zero.col.thresh <- 0.05 # Threshold value specified as double within [0,1].
#       Samples whose corresponding column contains more than zero.col.thresh x (number of markers) zeros are removed.
#       Default value for exome: 0.05.
scale.method <- 'median' # Function used to scale each column. Can be mean or median.


## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.n <- readRDS(file=here('02_TCGA_data_preparation/data', 'TCGA_WES_hg19_N.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/data', 'TCGA_WES_hg19_T.rds'))
doc <- bind_cols(doc.n, doc.t)

pos <- readRDS(file=here('02_TCGA_data_preparation/output/01_WES_probe_annotation', 'probes.hg19.annotated.rds'))

pos.used <- pos %>%
  filter(Masked==FALSE) %>%
  mutate(locus=paste(Chr, paste(Start, Stop, sep='-'), sep=':'))

dat <- doc[pos.used$locus, ]


## (There are no NAs in DOC data except Masked==TRUE loci in chrY PAR that have been removed. Treat 0 as NA instead.)
## 1-1. Samples (columns)
count.na <- function(list) {
  n.zero <- length(list[list==0])
  return(n.zero)
}

number.of.zeros <- lapply(dat[!grepl('X|Y', rownames(dat)), ], count.na)
qc.df <- data.frame(num.zero=unlist(number.of.zeros)) %>%
  rownames_to_column('SampleID') %>%
  select(SampleID, num.zero) %>%
  mutate(zero.ratio=num.zero/nrow(dat[!grepl('X|Y', rownames(dat)), ])) %>%
  mutate(too.many.zeros=ifelse(zero.ratio >= zero.col.thresh, TRUE, FALSE)) %>%
  left_join(sif, by='SampleID')
saveRDS(qc.df, file=here('02_TCGA_data_preparation/tmp/02_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'qc.df.rds'), compress=FALSE)

dat.sample.flt <- dat[, qc.df %>% filter(too.many.zeros==FALSE) %>% pull(SampleID)]

## 1-2. Probes (rows)
probes.with.bad.marker <- pos.used %>%
  mutate(zero.num=length(dat.sample.flt) - apply(dat.sample.flt, 1, Matrix::nnzero)) %>%
  mutate(zero.ratio=zero.num / length(dat.sample.flt)) %>%
  mutate(bad.marker=case_when(zero.ratio >= zero.row.thresh & Chr!='Y' ~ TRUE, TRUE ~ FALSE)) # Basically ChrY should be 0 in female samples, so bad markers on ChrY cannot be determined in this way.

dat.sample.flt.male.chrY <- dat.sample.flt[grepl('Y', rownames(dat.sample.flt)), qc.df %>% filter(SampleID %in% colnames(dat.sample.flt)) %>% filter(gender=='male') %>% pull(SampleID)]

probes.with.bad.marker.chrY <- pos.used %>%
  filter(Chr=='Y') %>%
  mutate(zero.num=length(dat.sample.flt.male.chrY) - apply(dat.sample.flt.male.chrY, 1, Matrix::nnzero)) %>%
  mutate(zero.ratio=zero.num / length(dat.sample.flt)) %>%
  mutate(bad.marker=case_when(zero.ratio >= zero.row.thresh ~ TRUE, TRUE ~ FALSE))

probe.num <- pos.used %>%
  group_by(Chr) %>%
  summarize(num.probes=n())

probes.with.bad.marker %>%
  filter(Chr!='Y') %>%
  bind_rows(probes.with.bad.marker.chrY) %>%
  filter(bad.marker==TRUE) %>%
  pull(Chr) %>%
  table.df() %>%
  as.data.frame.matrix() %>%
  rename(Chr='x', num.bad.probes='Freq') %>%
  mutate(Chr=factor(.$Chr, levels=.$Chr %>% unique() %>% gtools::mixedsort())) %>%
  arrange(Chr) %>%
  left_join(probe.num, by='Chr') %>%
  janitor::adorn_totals() %>%
  mutate(bad.rate=num.bad.probes/num.probes)
## There were no bad marker in chrY.

dat.sample.probe.flt <- dat.sample.flt[probes.with.bad.marker %>% filter(bad.marker==FALSE) %>% pull(locus), ]
saveRDS(dat.sample.probe.flt, file=here('02_TCGA_data_preparation/tmp/02_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'dat.sample.probe.flt.rds'), compress=FALSE)
