library(tidyverse)
library(here)

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

load(file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

tumors.to.be.analyzed <- c(tumors.female, tumors.male)

T <- doc.t[!grepl('X|Y', rownames(doc.t)), tumors.to.be.analyzed] %>%
  as.matrix()

n.num <- c(10, 50, 100, 200, 500)
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  female.samples <- normals.female[1:(n.num.i/2)]
  male.samples <- normals.male[1:(n.num.i/2)]

  N <- doc.n[!grepl('X|Y', rownames(doc.n)), c(female.samples, male.samples)] %>%
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
  Tn.medians <- Tn %>%
    apply(., 2, median)

  Tn.normalized <- t(t(Tn)- Tn.medians)

  saveRDS(Tn.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/10_PrototypeTangent', paste0('Tn_PrototypeTangent_', n.num.i, '.rds')), compress=FALSE)
}
