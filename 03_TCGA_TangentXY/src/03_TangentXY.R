library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

T <- doc.t %>%
  as.matrix()

## Tangent on autosome & chrX
T.autox <- T[!grepl('Y', rownames(T)),]

dimensions <- c(10, 30, 50, 100, 200, 500, 5000, 10441)
for (i in 1:length(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  N.lower.i <- readRDS(file=here('03_TCGA_TangentXY/output/02_SVD_on_normals', paste0('N.lower_', dim.i, 'dimensions.rds')))

  ## Get the origin in the normal subspace
  N.means <- apply(N.lower.i, 1, mean)
  N.lower.i0 <- N.lower.i - N.means
  T.autox0 <- T.autox - N.means

  Npi.autox <- pracma::pinv(N.lower.i0)
  weights.autox <- Npi.autox %*% T.autox0
  proj.autox <- N.lower.i0 %*% weights.autox
  Tn.autox <- T.autox0 - proj.autox

  ## Re-scaling after Tangent (by median)
  Tn.auto.medians <- Tn.autox[!grepl('X', rownames(Tn.autox)),] %>%
    apply(., 2, median)

  Tn.autox.scaled <- t(t(Tn.autox)- Tn.auto.medians)

  saveRDS(Tn.autox.scaled, file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds')), compress=FALSE)
}


## Tangent on chrY (sex-matched Tangent on male samples)
doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))

male.normals <- sif %>%
  filter(SampleID %in% colnames(doc.n) & Gender=='Male') %>%
  pull(SampleID)

male.tumors <- sif %>%
  filter(SampleID %in% colnames(doc.t) & Gender=='Male') %>%
  pull(SampleID)

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
Tn.male.medians <- Tn.male[!grepl('X|Y', rownames(Tn.male)),] %>%
  apply(., 2, median)

Tn.male.normalized <- t(t(Tn.male)- Tn.male.medians)

saveRDS(Tn.male.normalized, file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn_sexMatchedTangentOnMale.rds'), compress=FALSE)


## Combine Tangent-normalized autosomes & chrX and chrY
## Load Tangent-normalized autosomes & chrX signal that is normalized with re-constructed N with a desired number of latent factors
dimensions.to.be.used <- 50
Tn.autox.normalized <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dimensions.to.be.used, 'dimensions.rds')))

Tn.male.y <- Tn.male.normalized[grepl('Y', rownames(Tn.male.normalized)),]

## Adjust chrY to be relative to ploidy for each sample
Tn.male.y.adj <- Tn.male.y - 1

Tn.combined <- as.data.frame(Tn.autox.normalized) %>%
  bind_rows(as.data.frame(Tn.male.y.adj)) %>%
  replace(is.na(.), -10) %>%
  as.matrix()

saveRDS(Tn.combined, file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn.rds'), compress=FALSE)


## Divide Tn into each tumor type and save as R object for CBS
tumor.types <- sif$project %>% unique()
sample.ids <- colnames(Tn.combined)
for (i in seq_along(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))

  sample.ids.i <- sample.ids[grepl(paste0('^', tumor.type.i), sample.ids)]
  Tn.i <- Tn.combined[, sample.ids.i]

  saveRDS(Tn.i, file=here('03_TCGA_TangentXY/output/03_TangentXY/byTumorType', paste0('Tn_', tumor.type.i, '_', length(sample.ids.i), '.rds')), compress=FALSE)
}
