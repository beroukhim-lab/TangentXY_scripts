library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled by biological gender (female <-> male)
## 3. Remove outliers (Replace outlier signal with merginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by median of each sample 
## 6. Removal of common germline CNVs


## 2. Check if there are samples that have been mislabeled (female <-> male)
sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))
qc.df <- readRDS(file=here('05_CCLE_data_preparation/output/03_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'qc.df.rds'))

dat.sample.probe.flt <- readRDS(file=here('05_CCLE_data_preparation/output/03_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'dat.sample.probe.flt.rds'))

normal.samples <- sif %>%
  filter(DepMap_ID %in% colnames(dat.sample.probe.flt)) %>%
  filter(primary_disease=='Non-Cancerous') %>%
  filter(!is.na(sex)) %>%
  pull(DepMap_ID)

doc.xy.g <- dat.sample.probe.flt[grepl('chrX|chrY', rownames(dat.sample.probe.flt)), ] %>%
  select(normal.samples) %>%
  rownames_to_column('probe') %>%
  separate(data=., col=probe, into=c('Chr', 'pos'), sep=':') %>%
  separate(data=., col=pos, into=c('Start', 'Stop'), sep='-') %>%
  pivot_longer(names_to='DepMap_ID', values_to='signal', -c('Chr', 'Start', 'Stop')) %>%
  left_join(sif %>% select(DepMap_ID, cell_line_name, stripped_cell_line_name, sex, primary_disease, lineage), by='DepMap_ID')

g <- ggplot(doc.xy.g, aes(x=Chr, y=signal, col=sex)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_grid(lineage~DepMap_ID, scales='free_y') +
  theme_bw(base_size=20) +
  theme(strip.text.x=element_text(angle=90)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('05_CCLE_data_preparation/output/03_2_DOC_Preprocessing_removeSexMislabeledSamples', 'XY_signal_distribution.png'), dpi=100, width=18, height=6)

xy <- doc.xy.g %>%
  group_by(DepMap_ID, lineage, Chr, sex) %>%
  summarize(median=median(signal)) %>%
  ungroup() %>%
  pivot_wider(names_from='Chr', values_from='median') %>%
  mutate(y.x.ratio=chrY/chrX)

g <- ggplot(xy, aes(x=chrX, y=chrY)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_abline(data=. %>% filter(sex=='Male'), aes(slope=2, intercept=0), linetype='dashed') +
  geom_abline(data=. %>% filter(sex=='Male'), aes(slope=0.5, intercept=0), linetype='dashed') +
  geom_point(aes(col=lineage), size=2) +
  lims(x=c(0, NA), y=c(0, NA)) +
  ggrepel::geom_text_repel(data=. %>% filter(sex=='Male' & chrY < 10), aes(label=DepMap_ID), col='black', box.padding=1, show.legend=FALSE) +
  labs(x='X (median of depth)', y='Y (median of depth)') +
  facet_wrap(~sex) +
  theme_bw(base_size=20)
ggsave(g, file=here('05_CCLE_data_preparation/output/03_2_DOC_Preprocessing_removeSexMislabeledSamples', 'XY_signal.png'), dpi=100, width=16, height=6)

## Manually pick up suspicious samples
suspicious.samples <- c('ACH.000224', 'ACH.000531')
qc.df2 <- qc.df %>%
  mutate(suspicious.gender=case_when(DepMap_ID %in% suspicious.samples ~ TRUE, TRUE ~ FALSE))
saveRDS(qc.df2, file=here('05_CCLE_data_preparation/output/03_2_DOC_Preprocessing_removeSexMislabeledSamples', 'qc.df2.RData'), compress=FALSE)

samples.to.be.analyzed <- qc.df2 %>%
  filter(suspicious.gender==FALSE) %>%
  pull(DepMap_ID)

dat.flt <- dat.sample.probe.flt %>% select(samples.to.be.analyzed)
saveRDS(dat.flt, file=here('05_CCLE_data_preparation/output/03_2_DOC_Preprocessing_removeSexMislabeledSamples', 'dat.flt.rds'), compress=FALSE)
