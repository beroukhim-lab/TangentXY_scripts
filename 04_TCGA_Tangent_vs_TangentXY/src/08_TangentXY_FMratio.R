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
fm.table <- data.frame(case=1:8, female=c(5, 25, 50, 150, 350, 450, 475, 495), male=c(495, 475, 450, 350, 150, 50, 25, 5))
for (i in 1:nrow(fm.table)) {
  n.num.f <- fm.table[i, 'female']
  n.num.m <- fm.table[i, 'male']
  print(paste0(i, ' F:', n.num.f, ', M:', n.num.m))

  fm <- paste0('F', n.num.f, 'M', n.num.m)

  female.samples <- normals.female[1:n.num.f]
  male.samples <- normals.male[1:n.num.m]

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
  saveRDS(doc.n.transformed, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/LinearTransformation', paste0('N_LinearTransformed_', fm, '.rds')), compress=FALSE)
}

## SVD
for (i in 1:nrow(fm.table)) {
  n.num.f <- fm.table[i, 'female']
  n.num.m <- fm.table[i, 'male']
  print(paste0(i, ' F:', n.num.f, ', M:', n.num.m))

  fm <- paste0('F', n.num.f, 'M', n.num.m)

  doc.n.transformed <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/LinearTransformation', paste0('N_LinearTransformed_', fm, '.rds')))
  N.fm <- doc.n.transformed %>%
    as.matrix()

  ## SVD
  N.fm.svd <- svd(N.fm)
  saveRDS(N.fm.svd, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/SVD', paste0('N.svd_', fm, '.rds')), compress=FALSE)

  d.i <- N.fm.svd$d %>%
    as.data.frame() %>%
    setNames('r') %>%
    mutate(n=1:n()) %>%
    mutate(fm=fm)

  if (i==1) {
    d <- d.i
  } else {
    d <- d %>% bind_rows(d.i)
  }
}
saveRDS(d, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/SVD', 'd.rds'), compress=FALSE)

g <- ggplot(d %>% filter(n <= 30), aes(x=n, y=r)) +
  geom_line() +
  facet_wrap(~ fm, scale='free', ncol=1) +
  geom_vline(xintercept=seq(1, 10, 1), linetype='dashed') +
  scale_x_continuous(breaks=c(1, 5, 10, 20, 30), labels=c(1, 5, 10, 20, 30)) +
  theme_bw(base_size=20)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/SVD', 'latentFactorImportance_FMratio.png'), dpi=100, width=10, height=28)

n.lf <- c(1, 5, 10, 20, 50, 100, 500)
for (i in 1:nrow(fm.table)) {
  n.num.f <- fm.table[i, 'female']
  n.num.m <- fm.table[i, 'male']
  print(paste0(i, ' F:', n.num.f, ', M:', n.num.m))

  fm <- paste0('F', n.num.f, 'M', n.num.m)

  N.fm.svd <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/SVD', paste0('N.svd_', fm, '.rds')))

  for (j in seq_along(n.lf)) {
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
    ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/ChrX_signal', paste0(fm, '_ChrX_signalDistribution_afterTangent.png')), dpi=100, width=24, height=5)
  }
}

## Manually pick up a good number of latent factors to be used for each normal subspace based on the plot above
num.lf <- 50

## Run TangentXY
for (i in 1:nrow(fm.table)) {
  n.num.f <- fm.table[i, 'female']
  n.num.m <- fm.table[i, 'male']
  print(paste0(i, ' F:', n.num.f, ', M:', n.num.m))

  fm <- paste0('F', n.num.f, 'M', n.num.m)

  N.fm.svd <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/SVD', paste0('N.svd_', fm, '.rds')))

  dim.i <- num.lf

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

  saveRDS(Tn.normalized, file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/TangentXY', paste0('Tn_TangentXY_', fm, '.rds')), compress=FALSE)
}
