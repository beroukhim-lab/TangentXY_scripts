library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

gender.known.samples <- sif %>%
  filter(Gender!='NA') %>%
  pull(SampleID)

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds')) %>%
  select(any_of(gender.known.samples))

hg19 <- rCGH::hg19 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))

probes <- doc.n %>%
  as.data.frame() %>%
  mutate(locus=rownames(.)) %>%
  select(locus) %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  separate(col=pos, into=c('start', 'end'), sep='-') %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  mutate(probe=1:n()) %>%
  mutate(chrom=paste0('chr', chr)) %>%
  left_join(hg19 %>% select(chrom, centromerStart, centromerEnd), by=c('chr'='chrom')) %>%
  mutate(arm=case_when(end < centromerStart ~ 'p',
                        start > centromerEnd ~ 'q'))
saveRDS(probes, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds'), compress=FALSE)

signal.x <- doc.n[grepl('X', rownames(doc.n)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')
saveRDS(signal.x, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'signal.x.rds'), compress=FALSE)

g <- ggplot(signal.x, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(x=expression(paste({log[2]}, '[Relative copy number]', sep='')), title='ChrX signal (Before transformation)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1a.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1a.pdf'), width=8, height=8)


signal.y <- doc.n[grepl('Y', rownames(doc.n)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')
saveRDS(signal.y, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'signal.y.rds'), compress=FALSE)

g <- ggplot(signal.y, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(x=expression(paste({log[2]}, '[Relative copy number]', sep='')), title='ChrY signal') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c.pdf'), width=8, height=8)


## Linear transformation on male chrX signals
## Only male chrX is linear transformed so that it has the same mean and SD as female chrX
female.samples <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(Gender=='Female') %>%
  pull(SampleID)

male.samples <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(Gender=='Male') %>%
  pull(SampleID)

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

doc.n.xy.transformed <- doc.n[grepl('X|Y', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  mutate(signal=case_when(SampleID %in% male.samples & Chr=='X' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='SampleID', values_from='signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

doc.n.transformed <- doc.n[!grepl('X|Y', rownames(doc.n)), ] %>%
  bind_rows(doc.n.xy.transformed)
saveRDS(doc.n.transformed, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'TCGA_WES_hg19_N_Transformed.rds'), compress=FALSE)

## Check the signal distribution of chrX after linear transformation
signal.x.lt <- doc.n.transformed[grepl('X', rownames(doc.n.transformed)),] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.lt, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(x=expression(paste({log[2]}, '[Relative copy number]', sep='')), title='ChrX signal (After transformation)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1b.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1b.pdf'), width=8, height=8)
