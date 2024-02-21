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


## 3. Remove outliers

hg38.interval <- read.delim(here('05_CCLE_data_preparation/output/01_lifgOver_Hg19_to_GRCh38/exomeAgilentCaptureKit_GRCh38_data.txt'))
dat.flt <- readRDS(file=here('05_CCLE_data_preparation/output/03_2_DOC_Preprocessing_removeSexMislabeledSamples', 'dat.flt.rds'))


## A point is considered an outlier if ALL of the following apply:
## 1) value != NA
## 2) abs (value - median_of_last_5_points) > 0.3 
## 3) abs (value - median_of_next_5_points) > 0.3 
## 4) abs (ln(value) - ln(median_of_last_5_points)) > 4
## 5) abs (ln(value) - ln(median_of_next_5_points)) > 4

pos.dat.flt <- hg38.interval %>%
  mutate(Start=Start+1) %>%
  mutate(locus=paste(Chr, paste(Start, Stop, sep='-'), sep=':')) %>%
  filter(locus %in% rownames(dat.flt))

trailingN <- 5
remove.outliers <- function(list) {
  list.new <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(pos.dat.flt %>% select(locus, Chr, Arm2)) %>%
    group_by(Chr, Arm2)

  for (i in 1:trailingN) {
    column.name.last <- paste0('last', i)
    column.name.next <- paste0('next', i)
    list.new <- list.new %>%
      mutate(!!column.name.last := lag(signal, n=i, default=first(signal))) %>%
      mutate(!!column.name.next := lead(signal, n=i, default=last(signal))) %>%
      select(locus, Chr, Arm2, signal, starts_with('last'), starts_with('next'))
  }

  list.outlier.detected <- list.new %>%
    ungroup() %>%
    mutate(last_N_median=select(., starts_with('last')) %>% as.matrix %>% matrixStats::rowMedians()) %>%
    mutate(next_N_median=select(., starts_with('next')) %>% as.matrix %>% matrixStats::rowMedians()) %>%
    mutate(outlier=case_when(signal!=0 &
                              abs(signal - last_N_median) > 0.3 &
                              abs(signal - next_N_median) > 0.3 &
                              abs(log(signal) - log(last_N_median)) > 4 &
                              abs(log(signal) - log(next_N_median)) > 4 ~ TRUE,
                              TRUE ~ FALSE))

  list.replaced <- list.outlier.detected %>%
    mutate(signal.new=case_when(outlier==TRUE ~ select(., c('last1', 'signal', 'next1')) %>% as.matrix %>% matrixStats::rowMedians(), TRUE ~ signal)) %>%
    mutate(index=1:n())

  list.outlier.detected.return <- list.outlier.detected %>%
    column_to_rownames('locus') %>%
    select(outlier)

  list.replaced.return <- list.replaced %>%
    column_to_rownames('locus') %>%
    select(signal.new)

  return(list(new.signal=list.replaced.return, outlier.detected=list.outlier.detected.return))
}

dat.outlier.removed.list <- parallel::mclapply(dat.flt, remove.outliers, mc.cores=parallel::detectCores()-2)

dat.outlier.detected.list <- lapply(dat.outlier.removed.list, '[[', 'outlier.detected')
dat.outlier.detected <- dat.outlier.detected.list %>%
  bind_cols() %>%
  setNames(colnames(dat.flt))
saveRDS(dat.outlier.detected, file=here('05_CCLE_data_preparation/output/03_3_DOC_Preprocessing_removeOutliers', 'dat.outlier.detected.rds'), compress=FALSE)

dat.new.signal.list <- lapply(dat.outlier.removed.list, '[[', 'new.signal')
dat.outlier.removed <- dat.new.signal.list %>%
  bind_cols() %>%
  setNames(colnames(dat.flt))
saveRDS(dat.outlier.removed, file=here('05_CCLE_data_preparation/output/03_3_DOC_Preprocessing_removeOutliers', 'dat.outlier.removed.rds'), compress=FALSE)
