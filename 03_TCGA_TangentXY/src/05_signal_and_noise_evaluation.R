library(tidyverse)
library(here)

sif <- readRDS(file=here('02_TCGA_data_preparation/output/00_format_sif', 'sif.rds'))
probes <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds'))

gender.known.samples <- sif %>%
  filter(!is.na(Gender)) %>%
  pull(SampleID)

doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds')) %>%
  select(any_of(gender.known.samples))

female.tumors <- sif %>%
  filter(SampleID %in% colnames(doc.t)) %>%
  filter(Gender=='Female') %>%
  pull(SampleID)

doc.t[grepl('Y', rownames(doc.t)), female.tumors] <- NA

Tn.male.normalized <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn_sexMatchedTangentOnMale.rds'))
Tn.male.y <- Tn.male.normalized[grepl('Y', rownames(Tn.male.normalized)),] %>%
  as.data.frame()

## Signal and noise
## Make function for calculating signal and noise
calc.signal.noise <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(probes) %>%
    filter(!is.na(signal))

  y_nrow <- data %>%
    filter(chr=='Y') %>%
    nrow()
  
  signal.auto <- data %>%
    filter(!chr %in% c('X', 'Y')) %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  signal.x <- data %>%
    filter(chr=='X') %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  if (y_nrow!=0) {
    signal.y <- data %>%
      filter(chr=='Y') %>%
      group_by(chr, arm) %>%
      summarize(arm.median=median(signal)) %>%
      ungroup() %>%
      summarize(signal=sd(arm.median)) %>%
      pull(signal)
  } else {
    signal.y <- NA
  }

  signal.xy <- data %>%
    filter(chr %in% c('X', 'Y')) %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  signal.whole <- data %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sd(arm.median)) %>%
    pull(signal)

  noise.auto <- data %>% 
    filter(!chr %in% c('X', 'Y')) %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.x <- data %>% 
    filter(chr=='X') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  if (y_nrow!=0) {
    noise.y <- data %>% 
      filter(chr=='Y') %>%
      pull(signal) %>%
      diff() %>%
      abs() %>%
      median()
  } else {
    noise.y <- NA
  }
  
  noise.xy <- data %>% 
    filter(chr %in% c('X', 'Y')) %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  noise.whole <- data %>% 
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.auto=signal.auto,
                          signal.x=signal.x,
                          signal.y=signal.y,
                          signal.xy=signal.xy,
                          signal.whole=signal.whole,
                          noise.auto=noise.auto,
                          noise.x=noise.x,
                          noise.y=noise.y,
                          noise.xy=noise.xy,
                          noise.whole=noise.whole)

  return(result.df)
}


## Check chrX signal distribution after SVD (different number of dimensions)
dimensions <- c('Pre-norm', 10, 30, 50, 100, 200, 500, 5000, 10441)
for (i in seq_along(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  if (dim.i=='Pre-norm') {
    # Tn.i <- doc.t[!grepl('Y', rownames(doc.t)),]
    Tn.i <- doc.t
  } else {
    # Tn.i <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds'))) %>%
    #   as.data.frame()
    Tn.autox.normalized <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds'))) %>%
      as.data.frame()
    Tn.i <- Tn.autox.normalized %>%
      bind_rows(Tn.male.y)
  }

  # Tn.signal.noise.list.i <- parallel::mclapply(Tn.i, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.signal.noise.list.i <- lapply(Tn.i, calc.signal.noise)
  Tn.signal.noise.df.i <- Tn.signal.noise.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.signal.noise.list.i)) %>%
    mutate(dim=dim.i)

  if (i==1) {
    Tn.signal.noise.df <- Tn.signal.noise.df.i
  } else {
    Tn.signal.noise.df <- Tn.signal.noise.df %>% bind_rows(Tn.signal.noise.df.i)
  }

  rm(Tn.i)
  rm(Tn.signal.noise.list.i)
  rm(Tn.signal.noise.df.i)
  gc()
  gc()
}
saveRDS(Tn.signal.noise.df, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Tn.signal.noise.df.rds'), compress=FALSE)

signal.noise <- Tn.signal.noise.df %>%
  mutate(sn.auto=signal.auto/noise.auto) %>%
  mutate(sn.x=signal.x/noise.x) %>%
  mutate(sn.y=signal.y/noise.y) %>%
  mutate(sn.xy=signal.xy/noise.xy) %>%
  mutate(sn.whole=signal.whole/noise.whole) %>%
  mutate(dim=factor(.$dim, levels=.$dim %>% unique())) %>%
  left_join(sif, by='SampleID')

## Noise (violin plot)
g <- ggplot(signal.noise, aes(x=dim, y=noise.whole)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.1) +
  ylim(0, NA) +
  labs(x='# latent factors (k) in reference plane', y='Noise') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2b.png'), dpi=100, width=8, height=6)
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2b.pdf'), width=8, height=6)

## Signal-to-noise (violin plot)
g <- ggplot(signal.noise, aes(x=dim, y=sn.whole)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.1) +
  labs(x='# latent factors (k) in reference plane', y='Signal/Noise') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2c.png'), dpi=100, width=8, height=6)
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2c.pdf'), width=8, height=6)

signal.noise.prenorm <- signal.noise %>%
  filter(dim=='Pre-norm') %>%
  select(c(SampleID, matches('^signal|^noise'))) %>%
  column_to_rownames('SampleID') %>%
  setNames(paste0(colnames(.), '.prenorm')) %>%
  rownames_to_column('SampleID')

signal.noise2 <- signal.noise %>%
  filter(dim!='Pre-norm') %>%
  left_join(signal.noise.prenorm, by='SampleID')

g <- ggplot(signal.noise2, aes(x=signal.auto.prenorm, y=signal.auto)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, col='blue', linetype='dashed') +
  coord_fixed(ratio=1) +
  facet_wrap(~dim, nrow=1) +
  theme_classic(base_size=20)
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2d.png'), dpi=100, width=16, height=6)



chrx.diff <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(probes %>% filter(chr!='Y'))
  
  auto.median <- data %>%
    filter(chr!='X') %>%
    pull(signal) %>%
    median()
  
  chrx.p.median <- data %>%
    filter(chr=='X' & arm=='p') %>%
    pull(signal) %>%
    median()
  
  chrx.q.median <- data %>%
    filter(chr=='X' & arm=='q') %>%
    pull(signal) %>%
    median()
  
  delta.p <- chrx.p.median - auto.median
  delta.q <- chrx.q.median - auto.median

  result.df <- data.frame(delta.p=delta.p,
                          delta.q=delta.q)

  return(result.df)
}

for (i in seq_along(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  if (dim.i=='Pre-norm') {
    Tn.i <- doc.t[!grepl('Y', rownames(doc.t)),]
  } else {
    Tn.i <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds'))) %>%
      as.data.frame()
  }

  # chrx.diff.lilst.i <- parallel::mclapply(Tn.i, chrx.diff, mc.cores=parallel::detectCores()-2)
  chrx.diff.lilst.i <- lapply(Tn.i, chrx.diff)
  chrx.diff.df.i <- chrx.diff.lilst.i %>%
    bind_rows() %>%
    mutate(SampleID=names(chrx.diff.lilst.i)) %>%
    mutate(dim=dim.i)

  if (i==1) {
    chrx.diff.df <- chrx.diff.df.i
  } else {
    chrx.diff.df <- chrx.diff.df %>% bind_rows(chrx.diff.df.i)
  }

  rm(Tn.i)
  rm(chrx.diff.lilst.i)
  rm(chrx.diff.df.i)
  gc()
  gc()
}
saveRDS(chrx.diff.df, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'chrx.diff.df.rds'), compress=FALSE)

chrx.diff.df.l <- chrx.diff.df %>%
  pivot_longer(names_to='arm', values_to='diff', cols=c(delta.p, delta.q)) %>%
  left_join(sif, by='SampleID') %>%
  mutate(dim=factor(.$dim, levels=.$dim %>% unique())) %>%
  mutate(arm=sub('delta.', '', arm))

g <- ggplot(chrx.diff.df.l, aes(x=dim, y=diff)) +
  geom_hline(yintercept=0, col='red', linetype='dashed') +
  geom_hline(yintercept=-1, col='blue', linetype='dashed') +
  geom_boxplot(aes(fill=arm)) +
  facet_wrap(~Gender, nrow=2, scales='free_y') +
  labs(x='# latent factors (k) in reference plane', y='median(chrX arm) - median(autosomes)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2d.png'), dpi=100, width=12, height=8)

chrx.diff.prenorm <- chrx.diff.df.l %>%
  filter(dim=='Pre-norm') %>%
  select(SampleID, arm, diff) %>%
  rename(diff.prenorm=diff)

chrx.diff.df.l2 <- chrx.diff.df.l %>%
  filter(dim!='Pre-norm') %>%
  left_join(chrx.diff.prenorm, by=c('SampleID', 'arm'))

g <- ggplot(chrx.diff.df.l2, aes(x=diff.prenorm, y=diff)) +
  geom_point(aes(col=arm)) +
  geom_abline(intercept=0, slope=1, col='blue', linetype='dashed') +
  coord_fixed(ratio=1) +
  facet_wrap(~dim, nrow=1) +
  labs(title='median(chrX arm) - median(autosomes)', x='Pre-normalization', y='TamgemtXY') +
  theme_classic(base_size=20)
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2e.png'), dpi=100, width=16, height=6)
