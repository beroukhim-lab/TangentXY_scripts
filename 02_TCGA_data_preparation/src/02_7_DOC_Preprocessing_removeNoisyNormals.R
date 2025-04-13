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
## 7. Removal of noisy normal samples


## 7. Removal of noisy normal samples
sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

gender.known.samples <- sif %>%
  filter(Gender!='NA') %>%
  pull(SampleID)

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds')) %>%
  select(any_of(gender.known.samples))

pos <- readRDS(file=here('02_TCGA_data_preparation/output/01_WES_probe_annotation', 'probes.hg19.annotated.rds')) %>%
  mutate(locus=paste(Chr, paste(Start, Stop, sep='-'), sep=':'))

arm.medians.normal <- doc.n %>%
  rownames_to_column('locus') %>%
  left_join(pos %>% select(locus, Chr, arm), by='locus') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('locus', 'Chr', 'arm')) %>%
  group_by(SampleID, Chr, arm) %>%
  summarize(mean=mean(signal), median=median(signal))

arm.diff.normal <- doc.n %>%
  ungroup() %>%
  filter(!Chr %in% c('X', 'Y')) %>%
  mutate(gain=case_when(median > log2(1.5) ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(loss=case_when(median < log2(0.5) ~ TRUE, TRUE ~ FALSE)) %>%
  group_by(SampleID) %>%
  summarize(arm.gain=length(gain[gain==TRUE]), arm.loss=length(loss[loss==TRUE])) %>%
  mutate(total=arm.gain + arm.loss)



sample.of.interest <- 'ACC.OR.A5J1.NB'
sample.of.interest <- 'KIRC.AK.3436.NB'
sample.of.interest <- 'BLCA.FD.A43U.NB'

gender <- sif %>%
  filter(SampleID==sample.of.interest) %>%
  pull(Gender)

points.pr <- doc.n %>%
  select(sample.of.interest) %>%
  setNames('log2rcn') %>%
  rownames_to_column('locus') %>%
  left_join(pos, by='locus') %>%
  mutate(loci=row_number()) %>%
  mutate(chr.class=case_when(Chr %in% c(seq(1, 21, 2), 'X') ~ 'odd', TRUE ~ 'even')) %>%
  mutate(chr.class.col=case_when(chr.class=='odd' ~ '#0CB702', TRUE ~ '#000000'))

segments.pr <- arm.medians.normal %>%
  filter(SampleID==sample.of.interest) %>%
  left_join(points.pr %>% group_by(Chr, arm) %>% summarize(start=min(loci), end=max(loci)) %>% select(Chr, arm, start, end), by=c('Chr', 'arm'))

g <- ggplot(points.pr, aes(x=loci, y=log2rcn)) +
  ggrastr::rasterize(geom_point(aes(col=chr.class), size=1, show.legend=FALSE), dpi=300, dev='ragg_png') +
  geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
  geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
  geom_hline(yintercept=log2(1.5), col='red', linetype='dashed', linewidth=1) +
  geom_segment(data=segments.pr, aes(x=start, xend=end, y=median, yend=median), col='purple', linewidth=3) +
  scale_color_manual(values=c('odd'='limegreen', 'even'='black')) +
  coord_cartesian(ylim=c(-2.5, 2.5)) +
  labs(title=sample.of.interest, x='Genomic position', y=expression(paste({log[2]}, '[Relative copy-number]', sep=''))) +
  theme_bw(base_size=20) +
  theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(plot.title=element_text(hjust=0.5))
ggsave(g, file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', paste0(sample.of.interest, '_', gender, '.png')), bg='white', dpi=100, width=10, height=6)


doc.n.autosomes <- doc.n[!grepl('^X|^Y', rownames(doc.n)),]
doc.n.autosomes.t <- t(doc.n.autosomes)

pca <- prcomp(doc.n.autosomes.t)
saveRDS(pca, here('02_TCGA_data_preparation/output/02_7_DOC_Preprocessing_removeNoisyNormals', 'pca.rds'), compress=FALSE)

rob.pca <- rrcov::PcaHubert(doc.n.autosomes.t)
saveRDS(rob.pca, here('02_TCGA_data_preparation/output/02_7_DOC_Preprocessing_removeNoisyNormals', 'rob.pca.rds'), compress=FALSE)

grid.pca <- rrcov::PcaGrid(doc.n.autosomes.t)
saveRDS(grid.pca, here('02_TCGA_data_preparation/output/02_7_DOC_Preprocessing_removeNoisyNormals', 'grid.pca.rds'), compress=FALSE)


normal.samples.signal.noise <- readRDS(here('../../project/Tangent/20220725_Tangent_sex_TCGA_WES/output/SignalNoise/N.signal.noise.df.RData'))

pca.x <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column('SampleID') %>%
  left_join(sif, by='SampleID') %>%
  left_join(normal.samples.signal.noise, by='SampleID')

g <- ggplot(pca.x, aes(x=PC1, y=PC2)) +
  # geom_point(aes(col=signal.auto, shape=Gender), size=3) +
  geom_point(aes(col=type, shape=Gender), size=3) +
  facet_wrap(~project) +
  theme_bw(base_size=20)
ggsave(g, file=here('02_TCGA_data_preparation/output/02_7_DOC_Preprocessing_removeNoisyNormals', 'PCA.png'), bg='white', dpi=100, width=20, height=20)

rob.pca.scores <- rob.pca@scores %>%
  as.data.frame() %>%
  rownames_to_column('SampleID') %>%
  left_join(sif, by='SampleID') %>%
  left_join(normal.samples.signal.noise, by='SampleID')

g <- ggplot(rob.pca.scores, aes(x=PC1, y=PC2)) +
  # geom_point(aes(col=signal.auto, shape=Gender), size=3) +
  geom_point(aes(col=type, shape=Gender), size=3) +
  facet_wrap(~project) +
  theme_bw(base_size=20)
ggsave(g, file=here('02_TCGA_data_preparation/output/02_7_DOC_Preprocessing_removeNoisyNormals', 'Robust_PCA.png'), bg='white', dpi=100, width=20, height=20)





