library(tidyverse)
library(here)

sif <- readRDS(file=here('02_TCGA_data_preparation/output/00_format_sif', 'sif.rds'))

gender.known.samples <- sif %>%
  filter(!is.na(Gender)) %>%
  pull(SampleID)

## Check the distribution of ChrX median signal
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds')) %>%
  select(any_of(gender.known.samples))

T.chrX.signal.median <- doc.t %>%
  rownames_to_column('locus') %>%
  filter(grepl('^X', locus)) %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-'locus') %>%
  left_join(sif, by='SampleID') %>%
  group_by(SampleID, Gender) %>%
  summarize(median.signal=median(signal)) %>%
  mutate(dimension='Pre-norm')

dimensions <- c(10, 30, 50, 100, 200, 500, 5000, 10441)
for (i in seq_along(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  Tn.i <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', paste0('Tn_autox_svd_', dim.i, 'dimensions.rds'))) %>%
    as.data.frame()

  Tn.chrX.signal.median.i <- Tn.i %>%
    rownames_to_column('locus') %>%
    filter(grepl('^X', locus)) %>%
    pivot_longer(names_to='SampleID', values_to='signal', cols=-'locus') %>%
    left_join(sif, by='SampleID') %>%
    group_by(SampleID, Gender) %>%
    summarize(median.signal=median(signal)) %>%
    mutate(dimension=dim.i)

  if (i==1) {
    Tn.chrX.signal.median <- Tn.chrX.signal.median.i
  } else {
    Tn.chrX.signal.median <- bind_rows(Tn.chrX.signal.median, Tn.chrX.signal.median.i)
  }

  rm(Tn.i)
  rm(Tn.chrX.signal.median.i)
  gc()
  gc()
}

chrX.signal.median <- T.chrX.signal.median %>%
  bind_rows(Tn.chrX.signal.median %>% mutate(dimension=as.character(dimension))) %>%
  ungroup() %>%
  mutate(dimension=factor(.$dimension, levels=.$dimension %>% unique()))
saveRDS(chrX.signal.median, file=here('03_TCGA_TangentXY/output/04_normalized_chrX_signal', 'chrX.signal.median.rds'), compress=FALSE)

dimensions.labels <- c('Pre-norm'='Pre-norm', '10'='k=10', '30'='k=30', '50'='k=50', '100'='k=100', '200'='k=200', '500'='k=500', '5000'='k=5000', '10441'='k=10441\n(No SVD on N)')
g <- ggplot(chrX.signal.median %>% filter(Gender!='NA'), aes(x=median.signal, group=Gender)) +
  geom_density(aes(fill=Gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  scale_y_continuous(breaks=c(0, 10, 20)) +
  coord_flip() +
  lemon::facet_rep_wrap(~dimension, nrow=1, labeller=as_labeller(dimensions.labels), repeat.tick.labels=TRUE) +
  labs(x=expression(atop('Median of', paste({log[2]}, '[Relative copy-number]', sep=''))), y='Density') +
  # labs(x=expression(paste0('Median of ', '\n', {log[2]}, '[Relative copy-number]')), y='Density') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('03_TCGA_TangentXY/output/04_normalized_chrX_signal', 'Fig2a.png'), dpi=100, width=22, height=5)
ggsave(g, file=here('03_TCGA_TangentXY/output/04_normalized_chrX_signal', 'Fig2a.pdf'), width=22, height=5)
