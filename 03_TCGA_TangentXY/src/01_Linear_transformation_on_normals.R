library(tidyverse)
library(here)

sif <- readRDS(file=here('02_TCGA_data_preparation/output/00_format_sif', 'sif.rds'))

sif <- sif %>%
  mutate(age_bin=case_when(Age >= 10 & Age < 20 ~ '10-19',
                            Age >= 20 & Age < 30 ~ '20-29',
                            Age >= 30 & Age < 40 ~ '30-39',
                            Age >= 40 & Age < 50 ~ '40-49',
                            Age >= 50 & Age < 60 ~ '50-59',
                            Age >= 60 & Age < 70 ~ '60-69',
                            Age >= 70 & Age < 80 ~ '70-79',
                            Age >= 80 & Age < 90 ~ '80-89',
                            Age >= 90 ~ '90-'))

gender.known.samples <- sif %>%
  filter(!is.na(Gender)) %>%
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
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrX signal (Before shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1a.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1a.pdf'), width=8, height=8)


signal.y <- doc.n[grepl('Y', rownames(doc.n)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')
saveRDS(signal.y, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'signal.y.rds'), compress=FALSE)

g1 <- ggplot(signal.y, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrY signal (Before shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g1, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c_main.png'), width=8, height=8)
ggsave(g1, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c_main.pdf'), width=8, height=8)

g2 <- ggplot(signal.y, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25, show.legend=FALSE), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(-6, 6), ylim=c(0, 0.8)) +
  scale_x_continuous(breaks=c(-5, 0, 5)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8)) +
  theme_classic(base_size=20) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank())
ggsave(g2, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c_sub.png'), width=5, height=7)
ggsave(g2, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1c_sub.pdf'), width=5, height=7)


## Linear transformation on male chrX signals
## Only male chrX is linear transformed so that it has the same mean and SD as female chrX
female.normals <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(Gender=='Female') %>%
  pull(SampleID)

male.normals <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(Gender=='Male') %>%
  pull(SampleID)

female.x.mean <- doc.n[grepl('X', rownames(doc.n)),] %>%
  select(all_of(female.normals)) %>%
  as.matrix() %>%
  mean()

male.x.mean <- doc.n[grepl('X', rownames(doc.n)),] %>%
  select(all_of(male.normals)) %>%
  as.matrix() %>%
  mean()

doc.n.xy.shifted <- doc.n[grepl('X|Y', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  mutate(signal=case_when(SampleID %in% male.normals & Chr=='X' ~ signal - male.x.mean + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='SampleID', values_from='signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

doc.n.shifted <- doc.n[!grepl('X|Y', rownames(doc.n)), ] %>%
  bind_rows(doc.n.xy.shifted)
saveRDS(doc.n.shifted, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'TCGA_WES_hg19_N_Shifted.rds'), compress=FALSE)

## Linear transformation on male chrY signals
male.y.mean <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(all_of(male.normals)) %>%
  as.matrix() %>%
  mean()

male.y.mean.mode <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(all_of(male.normals)) %>%
  as.matrix() %>%
  apply(., 2, mean) %>%
  density() %>%
  {.$x[which.max(.$y)]}

doc.n.male.y.shifted <- doc.n[grepl('Y', rownames(doc.n)), male.normals] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  group_by(SampleID) %>%
  mutate(sample_mean=mean(signal), sample_sd=sd(signal)) %>%
  ungroup() %>%
  mutate(zs_signal=(signal - sample_mean + male.y.mean.mode)) %>%
  select(Chr, pos, SampleID, zs_signal) %>%
  pivot_wider(names_from='SampleID', values_from='zs_signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

male.y.shifted.mean <- doc.n.male.y.shifted %>%
  as.matrix() %>%
  mean()

doc.n.male.shifted <- doc.n[!grepl('Y', rownames(doc.n)), male.normals] %>%
  bind_rows(doc.n.male.y.shifted)
saveRDS(doc.n.male.shifted, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'TCGA_WES_hg19_N_Shifted_males.rds'), compress=FALSE)

## Check the signal distribution of chrX after linear transformation
signal.x.shifted <- doc.n.shifted[grepl('X', rownames(doc.n.shifted)),] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.shifted, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrX signal (After shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1b.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1b.pdf'), width=8, height=8)

## Check the signal distribution of chrY after linear transformation
signal.y.shifted <- doc.n.male.y.shifted %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.y.shifted, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  scale_fill_manual(values=c('Male'='#00BFC4')) +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrY signal (After shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1d.png'), width=8, height=8)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'FigS1d.pdf'), width=8, height=8)



## Shift on chrX of male tumors
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds')) %>%
  select(any_of(gender.known.samples))

signal.x.t <- doc.t[grepl('X', rownames(doc.t)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.t, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrX signal in tumor samples (Before shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'Tumor_chrX_signal_before_shift.png'), width=8, height=8)

male.tumors <- sif %>%
  filter(SampleID %in% colnames(doc.t)) %>%
  filter(Gender=='Male') %>%
  pull(SampleID)

doc.t.xy.shifted <- doc.t[grepl('X|Y', rownames(doc.t)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  mutate(signal=case_when(SampleID %in% male.tumors & Chr=='X' ~ signal - male.x.mean + female.x.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='SampleID', values_from='signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

doc.t.shifted <- doc.t[!grepl('X|Y', rownames(doc.t)), ] %>%
  bind_rows(doc.t.xy.shifted)
saveRDS(doc.t.shifted, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'TCGA_WES_hg19_T_Shifted.rds'), compress=FALSE)

signal.x.t.shifted <- doc.t.shifted[grepl('X', rownames(doc.t.shifted)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.t.shifted, aes(x=signal, group=SampleID)) +
  ggrastr::rasterize(geom_density(aes(fill=Gender), alpha=0.25), dpi=300, dev='ragg_png') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip(xlim=c(NA, 8.3)) +
  labs(x=expression(paste({log[2]}, '[Relative copy-number]', sep='')), y='Density', title='ChrX signal in tumor samples (After shifting)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'Tumor_chrX_signal_after_shift.png'), width=8, height=8)

male.x.mean.shifted <- doc.n.shifted[grepl('X', rownames(doc.n)),] %>%
  select(all_of(male.normals)) %>%
  as.matrix() %>%
  mean()

male.x.mean.t <- doc.t[grepl('X', rownames(doc.t)),] %>%
  select(all_of(male.tumors)) %>%
  as.matrix() %>%
  mean()

female.tumors <- sif %>%
  filter(SampleID %in% colnames(doc.t)) %>%
  filter(Gender=='Female') %>%
  pull(SampleID)

female.x.mean.t <- doc.t[grepl('X', rownames(doc.t)),] %>%
  select(all_of(female.tumors)) %>%
  as.matrix() %>%
  mean()

male.x.mean.t.shifted <- doc.t.shifted[grepl('X', rownames(doc.t.shifted)),] %>%
  select(all_of(male.tumors)) %>%
  as.matrix() %>%
  mean()


chr_medians.n <- doc.n %>%
  select(any_of(gender.known.samples)) %>%
  rownames_to_column('locus') %>%
  left_join(probes %>% select(locus, chr, arm), by='locus') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('locus', 'chr', 'arm')) %>%
  left_join(sif %>% select(SampleID, Gender), by='SampleID') %>%
  group_by(chr, arm, Gender) %>%
  summarize(median_signal=median(signal)) %>%
  ungroup() %>%
  mutate(type='Normal')

chr_medians.t <- doc.t %>%
  select(any_of(gender.known.samples)) %>%
  rownames_to_column('locus') %>%
  left_join(probes %>% select(locus, chr, arm), by='locus') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('locus', 'chr', 'arm')) %>%
  left_join(sif %>% select(SampleID, Gender), by='SampleID') %>%
  group_by(chr, arm, Gender) %>%
  summarize(median_signal=median(signal)) %>%
  ungroup() %>%
  mutate(type='Tumor')

chr_medians <- chr_medians.n %>%
  bind_rows(chr_medians.t) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique() %>% gtools::mixedsort()))


g <- ggplot(chr_medians %>% filter(!(Gender=='Female' & chr=='Y')), aes(x=chr, y=median_signal)) +
  geom_point(aes(col=Gender), size=3) +
  geom_hline(yintercept=0, col='red', linetype='dashed') +
  geom_hline(yintercept=-1, col='blue', linetype='dashed') +
  facet_grid(arm~type) +
  theme_bw(base_size=20)
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'Chr_signal_medians.png'), width=16, height=8)




## Distribution of shifted chrY means
shifted.male.y.medians <- doc.n.male.y.shifted %>%
  sapply(median) %>%
  as.data.frame() %>%
  setNames('chrY_median') %>%
  rownames_to_column('SampleID') %>%
  left_join(sif, by='SampleID')

g <- ggplot(shifted.male.y.medians, aes(x=chrY_median)) +
  geom_histogram(binwidth=0.01) +
  coord_flip() +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  theme_bw(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'Shifted_ChrY_medians.png'), width=12, height=8)

Tn.male.normalized <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn_sexMatchedTangentOnMale.rds'))
normalized.T.male.y.medians <- Tn.male.normalized %>%
  .[grepl('^Y', rownames(.)),] %>%
  sapply(median) %>%
  as.data.frame() %>%
  setNames('chrY_median') %>%
  rownames_to_column('SampleID') %>%
  left_join(sif, by='SampleID')

g <- ggplot(normalized.T.male.y.medians, aes(x=chrY_median)) +
  geom_histogram(binwidth=0.01) +
  coord_flip() +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  theme_bw(base_size=20) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'Normalized_Tumor_ChrY_medians.png'), width=12, height=8)


