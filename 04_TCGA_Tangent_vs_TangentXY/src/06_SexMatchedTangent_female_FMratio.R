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

fm.table <- data.frame(case=1:8, female=c(5, 25, 50, 150, 350, 450, 475, 495), male=c(495, 475, 450, 350, 150, 50, 25, 5))
for (i in 1:nrow(fm.table)) {
  n.num.i <- fm.table[i, 'female']
  print(paste(i, n.num.i))

  N.f <- N[, normals.female[1:n.num.i]]

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

  saveRDS(Tn.f.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/06_SexMatchedTangent_female_FMratio', paste0('Tn_SMTangent_female_', n.num.i, '.rds')), compress=FALSE)
}
