library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

load(file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

N <- doc.n %>%
  as.matrix()

T.m <- doc.t %>%
  select(tumors.male) %>%
  as.matrix()

fm.table <- data.frame(case=1:8, female=c(5, 25, 50, 150, 350, 450, 475, 495), male=c(495, 475, 450, 350, 150, 50, 25, 5))
for (i in 1:nrow(fm.table)) {
  n.num.i <- fm.table[i, 'male']
  print(paste(i, n.num.i))

  N.m <- N[, normals.male[1:n.num.i]]

  ## Get the origin in the normal subspace
  N.m.means <- apply(N.m, 1, mean)
  N.m0 <- N.m - N.m.means
  T.m0 <- T.m - N.m.means

  Npi <- pracma::pinv(N.m0)
  weights <- Npi %*% T.m0
  proj <- N.m0 %*% weights
  Tn.m <- T.m0 - proj

  ## Re-scaling after Tangent (by median)
  Tn.m.medians <- Tn.m[!grepl('X|Y', rownames(Tn.m)),] %>%
    apply(., 2, median)

  Tn.m.rescaled <- t(t(Tn.m)- Tn.m.medians)

  ## Adjust chrY so that it is relative to CN=2
  Tn.m.rescaled.auto <- Tn.m.rescaled[!grepl('X|Y', rownames(Tn.m.rescaled)),]
  Tn.m.rescaled.xy <- Tn.m.rescaled[grepl('X|Y', rownames(Tn.m.rescaled)),]
  Tn.m.rescaled.xy.adjusted <- Tn.m.rescaled.xy - 1
  Tn.m.normalized <- Tn.m.rescaled.auto %>%
    as.data.frame() %>%
    bind_rows(Tn.m.rescaled.xy.adjusted %>% as.data.frame()) %>%
    as.matrix()

  saveRDS(Tn.m.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/07_SexMatchedTangent_male_FMratio', paste0('Tn_SMTangent_male_', n.num.i, '.rds')), compress=FALSE)
}
