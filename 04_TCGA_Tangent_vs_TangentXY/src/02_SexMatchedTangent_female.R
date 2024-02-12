library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

load(file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

N <- doc.n[!grepl('Y', rownames(doc.n)), ] %>%
  as.matrix()

T.f <- doc.t[!grepl('Y', rownames(doc.t)), ] %>%
  select(tumors.female) %>%
  as.matrix()

n.num <- c(10, 50, 100, 200, 500)
for (i in 1:length(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  N.f <- N[, normals.female[1:(n.num.i/2)]]

  ## Get the origin in the normal subspace
  N.f.means <- apply(N.f, 1, mean)
  N.f0 <- N.f - N.f.means
  T.f0 <- T.f - N.f.means

  Npi <- pracma::pinv(N.f0)
  weights <- Npi %*% T.f0
  proj <- N.f0 %*% weights
  Tn.f <- T.f0 - proj

  ## Re-scaling after Tangent (by median)
  Tn.f.medians <- Tn.f[!grepl('X', rownames(Tn.f)),] %>%
    apply(., 2, median)

  Tn.f.normalized <- t(t(Tn.f)- Tn.f.medians)

  saveRDS(Tn.f.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/02_SexMatchedTangent_female', paste0('Tn_SMTangent_female_', n.num.i/2, '.rds')), compress=FALSE)
}
