library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

doc.t <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_T_QCed_commonCNVremoved.rds'))

T <- doc.t %>%
  as.matrix()

## Tangent on autosome & chrX
T.autox <- T[!grepl('chrY', rownames(T)),]

dimension <- 5
N.lower <- readRDS(file=here('06_CCLE_TangentXY/output/02_SVD_on_normals', paste0('N.lower_23normals_', dimension, 'dimensions.rds')))

## Get the origin in the normal subspace
N.lower.means <- apply(N.lower, 1, mean)
N.lower0 <- N.lower - N.lower.means
T.autox0 <- T.autox - N.lower.means

Npi.autox <- pracma::pinv(N.lower0)
weights.autox <- Npi.autox %*% T.autox0
proj.autox <- N.lower0 %*% weights.autox
Tn.autox <- T.autox0 - proj.autox

## Re-scaling after Tangent (by median)
Tn.auto.medians <- Tn.autox[!grepl('chrX', rownames(Tn.autox)),] %>%
  apply(., 2, median)

Tn.autox.normalized <- t(t(Tn.autox)- Tn.auto.medians)

saveRDS(Tn.autox.normalized, file=here('06_CCLE_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_23normals_', dimension, 'dimensions.rds')), compress=FALSE)


## Tangent on chrY (sex-matched Tangent on male samples)
doc.n <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_N_QCed_commonCNVremoved.rds'))

male.normals <- sif %>%
  filter(ModelID %in% colnames(doc.n)) %>%
  filter(Sex=='Male') %>%
  pull(ModelID)

male.tumors <- sif %>%
  filter(ModelID %in% colnames(doc.t)) %>%
  filter(Sex=='Male') %>%
  pull(ModelID)

N.male <- doc.n[, male.normals] %>%
  as.matrix()

T.male <- doc.t[, male.tumors] %>%
  as.matrix()

## Get the origin in the normal subspace
N.male.means <- apply(N.male, 1, mean)
N.male0 <- N.male - N.male.means
T.male0 <- T.male - N.male.means

## Run Tangent
Npi.male <- pracma::pinv(N.male0)
weights <- Npi.male %*% T.male0
proj.male <- N.male0 %*% weights
Tn.male <- T.male0 - proj.male

## Re-scaling after Tangent (by median)
Tn.male.medians <- Tn.male[!grepl('chrX|chrY', rownames(Tn.male)),] %>%
  apply(., 2, median)

Tn.male.normalized <- t(t(Tn.male)- Tn.male.medians)

saveRDS(Tn.male.normalized, file=here('06_CCLE_TangentXY/output/03_TangentXY', 'Tn_sexMatchedTangentOnMale.rds'), compress=FALSE)

## Combine Tangent-normalized autosomes & chrX and chrY

Tn.male.y <- Tn.male.normalized[grepl('chrY', rownames(Tn.male.normalized)),]

## Adjust chrY to be relative to ploidy for each sample
Tn.male.y.adj <- Tn.male.y - 1

Tn.combined <- as.data.frame(Tn.autox.normalized) %>%
  bind_rows(as.data.frame(Tn.male.y.adj)) %>%
  replace(is.na(.), -10) %>%
  as.matrix()

saveRDS(Tn.combined, file=here('06_CCLE_TangentXY/output/03_TangentXY', 'Tn.rds'), compress=FALSE)
