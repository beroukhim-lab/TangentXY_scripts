library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

probes <- readRDS(file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'probes.rds')) %>%
  rename(Chr=chr) %>%
  rename(chr=chrom) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

## Signal and noise
## Make function for calculating signal and noise
calc.signal.noise <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(probes)

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


## Before Tangent
doc.t <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_T_QCed_commonCNVremoved.rds'))

T.signal.noise.list <- parallel::mclapply(doc.t, calc.signal.noise, mc.cores=parallel::detectCores()-2)
T.signal.noise.df <- T.signal.noise.list %>%
  bind_rows() %>%
  mutate(ModelID=names(T.signal.noise.list)) %>%
  mutate(norm='Pre-Tangent')
saveRDS(T.signal.noise.df, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'T.signal.noise.df.rds'), compress=FALSE)

## After Tangent
Tn.combined <- readRDS(file=here('06_CCLE_TangentXY/output/03_TangentXY', 'Tn.rds')) %>%
  as.data.frame()

Tn.signal.noise.list <- parallel::mclapply(Tn.combined, calc.signal.noise, mc.cores=parallel::detectCores()-2)
Tn.signal.noise.df <- Tn.signal.noise.list %>%
  bind_rows() %>%
  mutate(ModelID=names(Tn.signal.noise.list)) %>%
  mutate(norm='Post-Tangent')
saveRDS(Tn.signal.noise.df, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'Tn.signal.noise.df.rds'), compress=FALSE)


## Signal (scatter plot)
signal <- T.signal.noise.df %>%
  select(ModelID, signal) %>%
  left_join(Tn.signal.noise.df %>% select(ModelID, signal), by='ModelID') %>%
  left_join(sif %>% select(ModelID, Sex, OncotreeLineage), by='ModelID')

g <- ggplot(signal %>% filter(Sex!='Unknown'), aes(x=signal.x, y=signal.y)) +
  geom_point(aes(col=OncotreeLineage, shape=Sex), alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  labs(title='Signal', x='Pre-Tangennt', y='Post-Tangent') +
  coord_fixed(ratio = 1, xlim=c(0, NA), ylim=c(0, NA)) +
  guides(shape=guide_legend(order=1), color=guide_legend(order=2)) +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'Signal_scatter.png'), dpi=100, width=16, height=8)


signal.noise <- T.signal.noise.df %>%
  bind_rows(Tn.signal.noise.df) %>%
  mutate(sn=signal/noise) %>%
  mutate(norm=factor(.$norm, levels=c('Pre-Tangent', 'Post-Tangent'))) %>%
  left_join(sif %>% select(ModelID, Sex, OncotreeLineage), by='ModelID')

## Signal (violin plot)
g <- ggplot(signal.noise %>% filter(Sex!='Unknown'), aes(x=norm, y=signal)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.05) +
  ylim(0, NA) +
  labs(y='Signal') +
  theme_classic(base_size=20) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'Signal_violin.png'), dpi=100, width=10, height=6)

## Noise (violin plot)
g <- ggplot(signal.noise %>% filter(Sex!='Unknown'), aes(x=norm, y=noise)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.05) +
  ylim(0, NA) +
  theme_classic(base_size=20) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'Noise_violin.png'), dpi=100, width=10, height=6)

## Signal-to-noise (violin plot)
g <- ggplot(signal.noise %>% filter(Sex!='Unknown'), aes(x=norm, y=sn)) +
  geom_violin(fill='red') +
  geom_boxplot(outlier.shape=NA, width=0.05) +
  ggbeeswarm::geom_beeswarm(aes(col=OncotreeLineage, shape=Sex), show.legend=F) +
  geom_line(aes(group=ModelID)) +
  labs(x='Number of dimensions in reference plane', y='S/N') +
  theme_classic(base_size=20) +
  theme(axis.title.x=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/06_signal_and_noise_evaluation', 'SN_violin.png'), dpi=100, width=10, height=6)
