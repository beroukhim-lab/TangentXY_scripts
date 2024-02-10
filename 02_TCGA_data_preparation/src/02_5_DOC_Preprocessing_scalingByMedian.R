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


## 5. Scaling by median of each sample
dat.log2 <- readRDS(file=here('02_TCGA_data_preparation/tmp/02_4_DOC_Preprocessing_log2Transformation', 'dat.log2.rds'))

medians <- sapply(dat.log2[!grepl('X|Y', rownames(dat.log2)), ], median)
dat.scaled <- t(t(dat.log2) - medians)

saveRDS(dat.scaled, file=here('02_TCGA_data_preparation/tmp/02_5_DOC_Preprocessing_scaleByMedian', 'dat.scaled.rds'), compress=FALSE)
