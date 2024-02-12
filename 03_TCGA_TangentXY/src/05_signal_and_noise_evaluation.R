library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))
probes <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds'))

doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

## Signal and noise
## Make function for calculating signal and noise
calc.signal.noise <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(probes %>% filter(chr!='Y'))

  signal <- data %>%
    filter(!chr %in% c('X', 'Y')) %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sds(arm.median)) %>%
    pull(signal)

  noise <- data %>% 
    filter(!chr %in% c('X', 'Y')) %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal=signal, noise=noise)

  return(result.df)
}


## Check chrX signal distribution after SVD (different number of dimensions)
dimensions <- c('Pre-norm', 10, 30, 50, 100, 200, 500, 5000, 10441)
for (i in seq_along(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  if (dim.i=='Pre-norm') {
    Tn.i <- doc.t[!grepl('Y', rownames(doc.t)),]
  } else {
    Tn.i <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds'))) %>%
      as.data.frame()
  }

  Tn.signal.noise.list.i <- parallel::mclapply(Tn.i, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.signal.noise.df.i <- Tn.signal.noise.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.signal.noise.list.i)) %>%
    mutate(dim=dim.i)

  if (i==1) {
    Tn.signal.noise.df <- Tn.signal.noise.df.i
  } else {
    Tn.signal.noise.df <- Tn.signal.noise.df %>% bind_rows(Tn.signal.noise.df.i)
  }
}
saveRDS(Tn.signal.noise.df, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Tn.signal.noise.df.rds'), compress=FALSE)


signal.noise <- Tn.signal.noise.df %>%
  mutate(sn=signal/noise) %>%
  mutate(dim=factor(.$dim, levels=.$dim %>% unique())) %>%
  left_join(sif, by='SampleID')

## Noise (violin plot)
g <- ggplot(signal.noise %>% filter(Gender!='NA'), aes(x=dim, y=noise)) +
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
g <- ggplot(signal.noise %>% filter(Gender!='NA'), aes(x=dim, y=sn)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.1) +
  labs(x='# latent factors (k) in reference plane', y='Signal/Noise') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2c.png'), dpi=100, width=8, height=6)
ggsave(g, file=here('03_TCGA_TangentXY/output/05_signal_and_noise_evaluation', 'Fig2c.pdf'), width=8, height=6)
