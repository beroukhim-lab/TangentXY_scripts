library(tidyverse)
library(here)

doc.n <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_T_QCed_commonCNVremoved.rds'))

N <- doc.n %>%
  as.matrix()

T <- doc.t %>%
  as.matrix()

## Get the origin in the normal subspace
N.means <- apply(N, 1, mean)
N0 <- N - N.means
T0 <- T - N.means

Npi <- pracma::pinv(N0)
weights <- Npi %*% T0
proj <- N0 %*% weights
Tn <- T0 - proj

## Re-scaling after Tangent (by median)
Tn.auto.medians <- Tn[!grepl('chrX|chrY', rownames(Tn)),] %>%
  apply(., 2, median)

Tn.normalized <- t(t(Tn)- Tn.auto.medians)

saveRDS(Tn.normalized, file=here('06_CCLE_TangentXY/output/04_PrototypeTangent', 'Tn_PrototypeTangent.rds'), compress=FALSE)
