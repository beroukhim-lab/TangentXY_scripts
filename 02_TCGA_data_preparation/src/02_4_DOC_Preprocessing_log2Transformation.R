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


## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation

## Options
floor.val.frac <- 0.001 # Threshold value specified as double within [0,1].
#       CN values will be floored at floor.val.frac x (data mean)

dat.outlier.removed <- readRDS(file=here('02_TCGA_data_preparation/tmp/02_3_DOC_Preprocessing_removeOutliers', 'dat.outlier.removed.rds'))

col.means <- sapply(dat.outlier.removed[!grepl('X|Y', rownames(dat.outlier.removed)),], mean) %>%
  as.data.frame() %>%
  setNames('mean')

floor.vals <- col.means %>%
  mutate(floor.vals=mean * floor.val.frac)
dat.0rep <- dat.outlier.removed
for (i in 1:length(dat.outlier.removed)) {
  print(i)
  min.val <- floor.vals[i, 'floor.vals']
  vals <- dat.outlier.removed[,i] %>% unlist()
  dat.0rep[vals < min.val, i] <- min.val
}


## 4.2. Log2 transformation

dat.log2 <- log2(dat.0rep)
saveRDS(dat.log2, file=here('02_TCGA_data_preparation/tmp/02_4_DOC_Preprocessing_log2Transformation', 'dat.log2.rds'), compress=FALSE)
