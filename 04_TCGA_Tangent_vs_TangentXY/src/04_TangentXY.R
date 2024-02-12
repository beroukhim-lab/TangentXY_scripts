library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

load(file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

tumors.to.be.analyzed <- c(tumors.female, tumors.male)

T <- doc.t[!grepl('Y', rownames(doc.t)), tumors.to.be.analyzed] %>%
  as.matrix()

## Z-score transformation
n.num <- c(10, 50, 100, 200, 500)
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  female.samples <- normals.female[1:(n.num.i/2)]
  male.samples <- normals.male[1:(n.num.i/2)]

  female.x.mean <- doc.n[grepl('X', rownames(doc.n)),] %>%
    select(female.samples) %>%
    as.matrix() %>%
    mean()

  male.x.mean <- doc.n[grepl('X', rownames(doc.n)),] %>%
    select(male.samples) %>%
    as.matrix() %>%
    mean()
  
  female.x.sd <- doc.n[grepl('X', rownames(doc.n)),] %>%
    select(female.samples) %>%
    as.matrix() %>%
    sd()
  
  male.x.sd <- doc.n[grepl('X', rownames(doc.n)),] %>%
    select(male.samples) %>%
    as.matrix() %>%
    sd()

  doc.n.x.transformed <- doc.n[grepl('X', rownames(doc.n)), c(female.samples, male.samples)] %>%
    rownames_to_column('locus') %>%
    separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
    pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
    mutate(signal=case_when(SampleID %in% male.samples & Chr=='X' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                            TRUE ~ signal)) %>%
    pivot_wider(names_from='SampleID', values_from='signal') %>%
    unite(col=locus, c('Chr', 'pos'), sep=':') %>%
    column_to_rownames('locus')

  doc.n.transformed <- doc.n[!grepl('X|Y', rownames(doc.n)), c(female.samples, male.samples)] %>%
    bind_rows(doc.n.x.transformed)
  saveRDS(doc.n.transformed, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/LinearTransformationOnNormals', paste0('N_LinearTransformed_', n.num.i, 'normals.rds')), compress=FALSE)
}

## SVD
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  doc.n.transformed <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/LinearTransformationOnNormals', paste0('N_LinearTransformed_', n.num.i, 'normals.rds')))
  N.fm <- doc.n.transformed %>%
    as.matrix()

  ## SVD
  N.fm.svd <- svd(N.fm)
  saveRDS(N.fm.svd, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/SVD', paste0('N.svd_female', n.num.i/2, 'male', n.num.i/2, '.rds')), compress=FALSE)

  d.i <- N.fm.svd$d %>%
    as.data.frame() %>%
    setNames('r') %>%
    mutate(n=1:n()) %>%
    mutate(num.sample=n.num.i)

  if (i==1) {
    d <- d.i
  } else {
    d <- d %>% bind_rows(d.i)
  }
}
saveRDS(d, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/SVD', 'd.rds'), compress=FALSE)

g <- ggplot(d %>% filter(n <= 30), aes(x=n, y=r)) +
  geom_line() +
  facet_wrap(~ num.sample, scale='free', ncol=1) +
  geom_vline(xintercept=seq(1, 10, 1), linetype='dashed') +
  scale_x_continuous(breaks=c(1, 5, 10, 20, 30), labels=c(1, 5, 10, 20, 30)) +
  theme_classic(base_size=20)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/SVD', 'latentFactorImportance.png'), dpi=100, width=8, height=20)

for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  N.fm.svd <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/SVD', paste0('N.svd_female', n.num.i/2, 'male', n.num.i/2, '.rds')))
  if (n.num.i=='10') {
    n.lf <- c(1, 2, 3, 4, 5, 10)
  } else if (n.num.i=='50') {
    n.lf <- c(1, 5, 10, 15, 20, 50)
  } else if (n.num.i=='100') {
    n.lf <- c(1, 5, 10, 20, 30, 100)
  } else if (n.num.i=='200') {
    n.lf <- c(1, 5, 10, 20, 50, 200)
  } else if (n.num.i=='500') {
    n.lf <- c(1, 5, 10, 20, 100, 500)
  }

  for (j in 1:length(n.lf)) {
    n.lf.j <- n.lf[j]
    print(paste0('-- dim:', n.lf.j))

    if (n.lf.j==1) {
      N.lower.i <- as.matrix(N.fm.svd$u[, n.lf.j]) %*% N.fm.svd$d[n.lf.j] %*% t(N.fm.svd$v[, n.lf.j])
    } else {
      N.lower.i <- N.fm.svd$u[, 1:n.lf.j] %*% diag(N.fm.svd$d[1:n.lf.j]) %*% t(N.fm.svd$v[, 1:n.lf.j])
    }

    ## Get the origin in the normal subspace
    N.means <- apply(N.lower.i, 1, mean)
    N.lower.i0 <- N.lower.i - N.means
    T0 <- T - N.means

    Npi <- pracma::pinv(N.lower.i0)
    weights <- Npi %*% T0
    proj <- N.lower.i0 %*% weights
    Tn <- T0 - proj

    ## Re-scaling after Tangent (by median)
    Tn.medians <- Tn[!grepl('X', rownames(Tn)),] %>%
      apply(., 2, median)

    Tn.normalized <- t(t(Tn)- Tn.medians)

    Tn.chrX.signal.j <- Tn.normalized %>%
      as.data.frame() %>%
      rownames_to_column('locus') %>%
      filter(grepl('^X', locus)) %>%
      pivot_longer(names_to='SampleID', values_to='signal', cols=-'locus') %>%
      mutate(dimension=n.lf.j)


    if (j==1) {
      Tn.chrX.signal <- Tn.chrX.signal.j
    } else {
      Tn.chrX.signal <- bind_rows(Tn.chrX.signal, Tn.chrX.signal.j)
    }
  }

  chrX.median <- Tn.chrX.signal %>%
    left_join(sif, by='SampleID') %>%
    group_by(SampleID, dimension) %>%
    mutate(median.signal=median(signal)) %>%
    select(SampleID, Gender, median.signal, dimension) %>%
    distinct()

  g <- ggplot(chrX.median, aes(x=median.signal)) +
    geom_density(aes(fill=Gender), alpha=0.25) +
    geom_vline(xintercept=0, col='red', linetype='dashed') +
    geom_vline(xintercept=-1, col='blue', linetype='dashed') +
    coord_flip() +
    facet_wrap(~dimension, nrow=1) +
    labs(title='ChrX signal distribution') +
    theme_classic(base_size=30)
  ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/ChrX_signal', paste0('Normals', n.num.i, '_ChrX_signalDistribution_afterTangent.png')), dpi=100, width=24, height=5)
}

## Manually pick up a good number of latent factors to be used for each normal subspace based on the plot above
num.lf <- data.frame(n.num=n.num, num.lf=c(10, 25, 50, 50, 50))

## Run ZTangentXY
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  N.fm.svd <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/SVD', paste0('N.svd_female', n.num.i/2, 'male', n.num.i/2, '.rds')))

  dim.i <- num.lf %>%
    filter(n.num==n.num.i) %>%
    pull(num.lf)

  N.lower.i <- N.fm.svd$u[, 1:dim.i] %*% diag(N.fm.svd$d[1:dim.i]) %*% t(N.fm.svd$v[, 1:dim.i])

  ## Get the origin in the normal subspace
  N.means <- apply(N.lower.i, 1, mean)
  N.lower.i0 <- N.lower.i - N.means
  T0 <- T - N.means

  Npi <- pracma::pinv(N.lower.i0)
  weights <- Npi %*% T0
  proj <- N.lower.i0 %*% weights
  Tn <- T0 - proj

  ## Re-scaling after Tangent (by median)
  Tn.medians <- Tn[!grepl('X', rownames(Tn)),] %>%
    apply(., 2, median)

  Tn.normalized <- t(t(Tn)- Tn.medians)

  saveRDS(Tn.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/TangentXY', paste0('Tn_TangentXY_', n.num.i, '.rds')), compress=FALSE)
}
