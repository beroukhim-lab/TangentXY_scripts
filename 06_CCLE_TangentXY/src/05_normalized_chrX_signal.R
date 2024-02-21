library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))

## Check the distribution of ChrX median signal
doc.t <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_T_QCed_commonCNVremoved.rds'))

T.chrX.signal <- doc.t %>%
  rownames_to_column('locus') %>%
  filter(grepl('^chrX', locus)) %>%
  pivot_longer(names_to='DepMapID', values_to='signal', cols=-'locus') %>%
  mutate(norm='Pre-Tangent')

Tn.normalized <- readRDS(file=here('06_CCLE_TangentXY/output/03_TangentXY', 'Tn.rds')) %>%
  as.data.frame()

Tn.chrX.signal <- Tn.normalized %>%
  rownames_to_column('locus') %>%
  filter(grepl('^chrX', locus)) %>%
  pivot_longer(names_to='DepMapID', values_to='signal', cols=-'locus') %>%
  mutate(norm='Post-Tangent')

chrX.signal <- T.chrX.signal %>%
  bind_rows(Tn.chrX.signal) %>%
  mutate(norm=factor(.$norm, levels=c('Pre-Tangent', 'Post-Tangent'))) %>%
  left_join(sif %>% select(DepMap_ID, sex, lineage), by=c('DepMapID'='DepMap_ID'))

chrX.signal.median <- chrX.signal %>%
  filter(!is.na(sex)) %>%
  group_by(DepMapID, norm, sex) %>%
  summarize(median.signal=median(signal))

g <- ggplot(chrX.signal.median, aes(x=median.signal, group=sex)) +
  geom_density(aes(fill=sex), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  geom_vline(xintercept=-0.9, col='black', linetype='dashed') +
  coord_flip() +
  facet_wrap(~norm, nrow=2) +
  labs(title='ChrX median signal distribution') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/05_normalized_chrX_signal', 'ChrX_medianDistribution_afterTangent.png'), dpi=100, width=12, height=10)
