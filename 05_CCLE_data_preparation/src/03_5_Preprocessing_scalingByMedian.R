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


## 5. Scaling by mean/median of each sample
dat.log2 <- readRDS(file=here('05_CCLE_data_preparation/output/03_4_DOC_Preprocessing_log2Transformation', 'dat.log2.rds'))

means <- sapply(dat.log2[!grepl('chrX|chrY', rownames(dat.log2)), ], mean)
dat.scaled <- t(t(dat.log2) - means)

saveRDS(dat.scaled, file=here('05_CCLE_data_preparation/output/03_5_DOC_Preprocessing_scalingByMedian', 'dat.scaled.rds'), compress=FALSE)
