library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

doc.n <- readRDS(file=here('05_CCLE_data_preparation/output/03_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'CCLE_WES_hg38_N_QCed_commonCNVremoved.rds'))

hg38 <- rCGH::hg38 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))

probes <- doc.n %>%
  mutate(locus=rownames(.)) %>%
  select(locus) %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  separate(col=pos, into=c('start', 'end'), sep='-') %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  mutate(probe=1:n()) %>%
  mutate(chrom=sub('chr', '', chr)) %>%
  left_join(hg38 %>% select(chrom, centromerStart, centromerEnd), by='chrom') %>%
  mutate(arm=case_when(end < centromerStart + 1500000 ~ 'p',
                        start > centromerEnd - 1500000 ~ 'q'))
saveRDS(probes, file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'probes.rds'), compress=FALSE)

## Check the signal distribution of chrX and chrY in normal samples (Before Z-score conversion)
signal.x <- doc.n[grepl('chrX', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  left_join(probes, by='locus') %>%
  pivot_longer(names_to='ModelID', values_to='signal', cols=colnames(doc.n)) %>%
  left_join(sif %>% select(ModelID, CellLineName, StrippedCellLineName, Sex, OncotreePrimaryDisease, OncotreeLineage), by='ModelID')

g <- ggplot(signal.x, aes(x=signal, group=ModelID)) +
  geom_density(aes(fill=Sex), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='ChrX signal (Before transformation)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'NormalSamples_ChrX_signalDistribution.png'), dpi=100, width=8, height=8)


## Linear transformation on male chrX signals
## Only male chrX is linear transformed so that it has the same mean and SD as female chrX
female.samples <- sif %>%
  filter(ModelID %in% colnames(doc.n)) %>%
  filter(Sex=='Female') %>%
  pull(ModelID)

male.samples <- sif %>%
  filter(ModelID %in% colnames(doc.n)) %>%
  filter(Sex=='Male') %>%
  pull(ModelID)

female.x.mean <- doc.n[grepl('chrX', rownames(doc.n)),] %>%
  select(female.samples) %>%
  as.matrix() %>%
  mean()

male.x.mean <- doc.n[grepl('chrX', rownames(doc.n)),] %>%
  select(male.samples) %>%
  as.matrix() %>%
  mean()

female.x.sd <- doc.n[grepl('chrX', rownames(doc.n)),] %>%
  select(female.samples) %>%
  as.matrix() %>%
  sd()

male.x.sd <- doc.n[grepl('chrX', rownames(doc.n)),] %>%
  select(male.samples) %>%
  as.matrix() %>%
  sd()

doc.n.xy.transformed <- doc.n[grepl('chrX|chrY', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='ModelID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  mutate(signal=case_when(ModelID %in% male.samples & Chr=='chrX' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='ModelID', values_from='signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

doc.n.transformed <- doc.n[!grepl('chrX|chrY', rownames(doc.n)), ] %>%
  bind_rows(doc.n.xy.transformed)
saveRDS(doc.n.transformed, file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'CCLE_WES_hg38_N_Transformed.rds'), compress=FALSE)

## Check the signal distribution of chrX after Z-score conversion
signal.x.lt <- doc.n.transformed[grepl('chrX', rownames(doc.n.transformed)),] %>%
  pivot_longer(names_to='ModelID', values_to='signal', cols=everything()) %>%
  left_join(sif %>% select(ModelID, CellLineName, StrippedCellLineName, Sex, OncotreePrimaryDisease, OncotreeLineage), by='ModelID')

g <- ggplot(signal.x.lt, aes(x=signal, group=ModelID)) +
  geom_density(aes(fill=Sex), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='ChrX signal (After transformation)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'NormalSamples_ChrX_transformedSignalDistribution.png'), dpi=100, width=8, height=8)
