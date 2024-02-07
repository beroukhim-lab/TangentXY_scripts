library(tidyverse)
library(here)

sif <- read.delim(file=here('../20220725_Tangent_sex_TCGA_WES/output', 'sif.txt'))
pos <- readRDS(file=here('../20220725_Tangent_sex_TCGA_WES/output/DOC/merged', 'probes.hg19.annotated.RData'))

gender.known.samples <- sif %>%
  filter(gender!='na') %>%
  pull(SampleID)

doc.n <- readRDS(file=here('output/07_7_preProcessing_7_removeNTMislabeledSamples', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.RData')) %>%
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
saveRDS(probes, file=here('output/08_ZscoreTransformation', 'probes.RData'), compress=FALSE)
# probes <- readRDS(file=here('output/08_ZscoreTransformation', 'probes.RData'))

## Check the signal distribution of chrX and chrY in normal samples (Before Z-score conversion)
# signal.all <- doc.n %>%
#   rownames_to_column('locus') %>%
#   pivot_longer(names_to='SampleID', values_to='signal', cols=-'locus') %>%
#   separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
#   left_join(sif, by='SampleID')

# g <- ggplot(signal.all, aes(x=signal, group=SampleID)) +
#   geom_density(aes(fill=gender), alpha=0.25) +
#   geom_vline(xintercept=0, col='red', linetype='dashed') +
#   geom_vline(xintercept=-1, col='blue', linetype='dashed') +
#   facet_wrap(~Chr, nrow=4) +
#   #lims(x=c(-10, 6.3)) +
#   coord_flip(ylim=c(NA, 0.8)) +
#   # coord_flip() +
#   labs(title='Normal samples signal distribution\n(scaled & log2 transformed)') +
#   theme_bw(base_size=20)
# show(g)
# ggsave(g, file=here('output/ZscoreTransformation', 'NormalSamples_SignalDistribution.png'), dpi=100, width=16, height=12)

# sd.all <- signal.all %>%
#   group_by(SampleID, Chr) %>%
#   mutate(sd=sd(signal)) %>%
#   ungroup() %>%
#   select(-signal) %>%
#   distinct()

# g <- ggplot(sd.all, aes(x=Chr, y=sd)) +
#   geom_boxplot(aes(fill=gender)) +
#   theme_bw(base_size=20)
# show(g)

signal.x <- doc.n[grepl('X', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  left_join(probes, by='locus') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=colnames(doc.n)) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='Normal samples ChrX signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_signalDistribution2.png'), dpi=100, width=8, height=8)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_signalDistribution2.pdf'), width=8, height=8, compression='lzw')

g <- ggplot(signal.x, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_wrap(~gender, nrow=1) +
  labs(title='Normal samples ChrX signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_signalDistribution_genderSeparated2.png'), dpi=100, width=12, height=8)

g <- ggplot(signal.x, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_grid(arm~gender) +
  labs(title='Normal samples ChrX signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_signalDistribution_genderSeparated_arm.png'), dpi=100, width=12, height=12)

signal.x.mean.median <- signal.x %>%
  group_by(SampleID) %>%
  summarize(mean=mean(signal), median=median(signal)) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.mean.median, aes(x=mean, y=median)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_point(aes(col=gender)) +
  coord_fixed(ratio = 1) +
  labs(title='Normal samples ChrX signal\n(scaled & log2 transformed)') +
  theme_bw(base_size=30)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_signal_meanVsMedian.png'), dpi=100, width=12, height=12)

## Define a function to find a peak value for chrX signal distribution
detect.peak <- function(df) {
  signal.density <- density(df$signal)
  max.dens <- which.max(signal.density$y)
  x <- signal.density$x[max.dens]
  return(x)
}

# Check if the peak of chrX distribution
peaks <- signal.x %>%
  group_by(SampleID) %>%
  nest() %>%
  mutate(peak=map_dbl(data, detect.peak)) %>%
  left_join(sif, by='SampleID') %>%
  select(-data) %>%
  arrange(abs(peak)) %>%
  ungroup()

g <- ggplot(peaks %>% mutate(SampleID=factor(.$SampleID, levels=.$SampleID)), aes(x=SampleID, y=peak)) +
  geom_point(aes(col=gender, shape=type), size=3) +
  scale_x_discrete(expand=expansion(add=200)) +
  geom_hline(yintercept=0.5, col='red', linetype='dashed') +
  geom_hline(yintercept=-0.5, col='red', linetype='dashed') +
  labs(x='SampleID (in order of |peak|)', col='Gender', shape='Type') +
  theme_bw(base_size=30) +
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank())
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_peaks.png'), dpi=100, width=12, height=8)


signal.y <- doc.n[grepl('Y', rownames(doc.n)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.y, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  # lims(x=c(-5, 6.3)) +
  labs(title='Normal samples ChrY signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrY_signalDistribution2.png'), dpi=100, width=8, height=8)
# ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrY_signalDistribution_zoomin2.png'), dpi=100, width=8, height=8)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrY_signalDistribution2.pdf'), width=8, height=8)

g <- ggplot(signal.y %>% filter(gender=='male'), aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='Normal samples ChrY signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrY_signalDistribution_Male.png'), dpi=100, width=8, height=8)

## Also check the signal distribution in tumor samples just in case
doc.t <- readRDS(file=here('output/07_7_preProcessing_7_removeNTMislabeledSamples', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.RData'))

signal.x.t <- doc.t[grepl('X', rownames(doc.t)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.x.t, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  scale_fill_manual(values=c('#F8766D', '#00BFC4', '#00BA38')) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_wrap(~gender) +
  labs(title='Tumor samples ChrX signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'TumorSamples_ChrX_signalDistribution2.png'), dpi=100, width=16, height=8)

signal.y.t <- doc.t[grepl('Y', rownames(doc.t)), ] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(signal.y.t, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  scale_fill_manual(values=c('#F8766D', '#00BFC4', '#00BA38')) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  # lims(x=c(-5, 6.3)) +
  labs(title='Tumor samples ChrY signal distribution\n(scaled & log2 transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'TumorSamples_ChrY_signalDistribution2.png'), dpi=100, width=8, height=8)
# ggsave(g, file=here('output/08_ZscoreTransformation', 'TumorSamples_ChrY_signalDistribution_zoomin2.png'), dpi=100, width=8, height=8)



## Make biased Z-score
## Firstly, compare mean and median of chrX and chrY
xy.summary <- doc.n[grepl('X|Y', rownames(doc.n)),] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  group_by(SampleID, Chr) %>%
  summarize(mean=mean(signal), median=median(signal)) %>%
  left_join(sif, by='SampleID')

g <- ggplot(xy.summary, aes(x=mean, y=median)) +
  geom_point(aes(col=gender)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  facet_wrap(~Chr, scale='free') +
  theme_bw(base_size=20)
show(g)



## Rameen's method
## Only male chrX is Z-score transformed so that it has the same mean and SD as female chrX
## Female chrY is also transformed just in case
female.samples <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(gender=='female') %>%
  pull(SampleID)
male.samples <- sif %>%
  filter(SampleID %in% colnames(doc.n)) %>%
  filter(gender=='male') %>%
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

female.y.mean <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(female.samples) %>%
  as.matrix() %>%
  mean()
male.y.mean <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(male.samples) %>%
  as.matrix() %>%
  mean()
female.y.sd <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(female.samples) %>%
  as.matrix() %>%
  sd()
male.y.sd <- doc.n[grepl('Y', rownames(doc.n)),] %>%
  select(male.samples) %>%
  as.matrix() %>%
  sd()

doc.n.xy.ztransformed <- doc.n[grepl('X|Y', rownames(doc.n)), ] %>%
  rownames_to_column('locus') %>%
  separate(col=locus, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  mutate(signal=case_when(SampleID %in% male.samples & Chr=='X' ~ ((signal - male.x.mean)/male.x.sd) * female.x.sd + female.x.mean,
                          # SampleID %in% female.samples & Chr=='Y' ~ ((signal - female.y.mean)/female.y.sd) * male.y.sd + male.y.mean,
                          TRUE ~ signal)) %>%
  pivot_wider(names_from='SampleID', values_from='signal') %>%
  unite(col=locus, c('Chr', 'pos'), sep=':') %>%
  column_to_rownames('locus')

doc.n.ztransformed <- doc.n[!grepl('X|Y', rownames(doc.n)), ] %>%
  bind_rows(doc.n.xy.ztransformed)
# write.table(doc.n.ztransformed, file=here('output/08_ZscoreTransformation', 'TCGA_WES_hg19_N_Ztransformed.txt'), sep='\t', row.names=FALSE)
saveRDS(doc.n.ztransformed, file=here('output/08_ZscoreTransformation', 'TCGA_WES_hg19_N_Ztransformed.RData'), compress=FALSE)
# doc.n.ztransformed <- readRDS(file=here('output/08_ZscoreTransformation', 'TCGA_WES_hg19_N_Ztransformed.RData'))

## Check the signal distribution of chrX after Z-score conversion
z.signal.x <- doc.n.ztransformed[grepl('X', rownames(doc.n.ztransformed)),] %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
  left_join(sif, by='SampleID')

g <- ggplot(z.signal.x, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='Normal samples ChrX signal distribution\n(Z-score transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_ZscoreDistribution2.png'), dpi=100, width=8, height=8)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_ZscoreDistribution2.pdf'), width=8, height=8)

g <- ggplot(z.signal.x, aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  facet_wrap(~gender, nrow=1) +
  labs(title='Normal samples ChrX signal distribution\n(Z-score transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_ZscoreDistribution_genderSeparated2.png'), dpi=100, width=12, height=8)

z.peaks <- z.signal.x %>%
  group_by(SampleID) %>%
  nest() %>%
  mutate(peak=map_dbl(data, detect.peak)) %>%
  left_join(sif, by='SampleID') %>%
  select(-data) %>%
  arrange(abs(peak)) %>%
  ungroup()

g <- ggplot(z.peaks %>% mutate(SampleID=factor(.$SampleID, levels=.$SampleID)), aes(x=SampleID, y=peak)) +
  geom_point(aes(col=gender, shape=type), size=3) +
  scale_x_discrete(expand=expansion(add=200)) +
  geom_hline(yintercept=0.5, col='red', linetype='dashed') +
  geom_hline(yintercept=-0.5, col='red', linetype='dashed') +
  labs(x='SampleID (in order of |peak|)', col='Gender', shape='Type') +
  theme_bw(base_size=30) +
  theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank())
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_Zscore_peaks.png'), dpi=100, width=12, height=8)

# Create the function.
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

z.signal.x.mean.median <- z.signal.x %>%
  group_by(SampleID) %>%
  summarize(mean=mean(signal), median=median(signal), mode=getmode(signal)) %>%
  left_join(sif, by='SampleID')

g <- ggplot(z.signal.x.mean.median, aes(x=mean, y=median)) +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_point(aes(col=gender)) +
  coord_fixed(ratio = 1) +
  labs(title='Normal samples ChrX signal\n(Z-score transformed)') +
  theme_bw(base_size=30)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_ZscoreTransformedSignal_meanVsMedian.png'), dpi=100, width=12, height=12)

gpub <- ggpubr::ggscatterhist(
  z.signal.x.mean.median, x='mean', y='median',
  color='gender',
  margin.plot='histogram',
  margin.params=list(fill='gender'),
  title='Z-score transformed chrX signal: Mean vs Median',
  ggtheme=theme_bw(base_size=30)
)
gpub$sp <- gpub$sp +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  geom_vline(xintercept=0, linetype='dashed') +
  geom_hline(yintercept=0, linetype='dashed')
ggpubr::ggexport(gpub, filename=here('output/08_ZscoreTransformation', 'NormalSamples_ChrX_ZscoreTransformedSignal_meanVsMedian2.png'), width=1000, height=1000)

g <- ggplot(z.signal.x %>% filter(SampleID=='KIRC.DV.5576.NB'), aes(x=signal, group=SampleID)) +
  geom_density(aes(fill=gender), alpha=0.25) +
  # geom_vline(data=z.signal.x %>% filter(SampleID=='KIRC.DV.5576.NB'), mapping=aes(xintercept=median)) +
  geom_vline(xintercept=0.3443613, col='blue', linetype='dashed') +
  geom_vline(xintercept=0, col='red', linetype='dashed') +
  geom_vline(xintercept=-1, col='blue', linetype='dashed') +
  coord_flip() +
  labs(title='Normal samples ChrX signal distribution\n(Z-score transformed)') +
  theme_bw(base_size=20)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'test.png'), dpi=100, width=8, height=8)


# z.signal.y <- doc.n.ztransformed[grepl('Y', rownames(doc.n.ztransformed)),] %>%
#   pivot_longer(names_to='SampleID', values_to='signal', cols=everything()) %>%
#   left_join(sif, by='SampleID')

# g <- ggplot(z.signal.y, aes(x=signal, group=SampleID)) +
#   geom_density(aes(fill=gender), alpha=0.25) +
#   geom_vline(xintercept=0, col='red', linetype='dashed') +
#   geom_vline(xintercept=-1, col='blue', linetype='dashed') +
#   # lims(x=c(-10, 6.3), y=c(NA, 0.8)) +
#   coord_flip() +
#   labs(title='Normal samples ChrY signal distribution\n(Z-score transformed)') +
#   theme_bw(base_size=20)
# # show(g)
# ggsave(g, file=here('output', 'NormalSamples_ChrY_ZscoreDistribution.png'), dpi=100, width=8, height=8)



## Compare the entire signal before and after Z score transformation
doc.n.fm <- doc.n %>%
  select(1:2) %>%
  rownames_to_column('probe') %>%
  separate(col=probe, into=c('Chr', 'pos'), sep=':') %>%
  mutate(index=1:n()) %>%
  select(-'pos') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'index')) %>%
  mutate(signal.type='Before Z-score transformation')

doc.n.z.fm <- doc.n.ztransformed %>%
  select(1:2) %>%
  rownames_to_column('probe') %>%
  separate(col=probe, into=c('Chr', 'pos'), sep=':') %>%
  mutate(index=1:n()) %>%
  select(-'pos') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'index')) %>%
  mutate(signal.type='After Z-score transformation')

doc.n.ba <- doc.n.fm %>%
  bind_rows(doc.n.z.fm) %>%
  # filter(Chr %in% c('X', 'Y')) %>%
  left_join(sif, by='SampleID') %>%
  mutate(signal.type=factor(.$signal.type, levels=c('Before Z-score transformation', 'After Z-score transformation'))) %>%
  mutate(Chr.class=case_when(Chr %in% c('1', '3', '5', '7', '9', '11', '13', '15', '17', '19', '21', 'X') ~ 'odd',
                            Chr %in% c('2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '22', 'Y') ~ 'even'))

doc.n.ba.medians <- doc.n.ba %>%
  group_by(SampleID, Chr, signal.type) %>%
  mutate(median=median(signal)) %>%
  mutate(i.min=min(index)) %>%
  mutate(i.max=max(index)) %>%
  mutate(n=1:n()) %>%
  filter(n==1)

g <- ggplot(doc.n.ba, aes(x=index, y=signal)) +
  geom_point(aes(col=Chr.class), show.legend=FALSE) +
  geom_hline(yintercept=0, col='darkgray', linetype='dashed') +
  geom_hline(yintercept=-1, col='blue', linetype='dashed') +
  geom_segment(data=doc.n.ba.medians, aes(x=i.min, xend=i.max, y=median, yend=median), col='red') +
  scale_color_manual(values=c('odd'='green', 'even'='black')) +
  # ylim(-5, 2) +
  # labs(title=array) +
  facet_grid(gender~signal.type) +
  theme_bw(base_size=15)
# show(g)
ggsave(g, file=here('output/08_ZscoreTransformation', 'NormalSamples_signalDistribution_beforeAfterZscoreTransformation2.png'), dpi=100, width=14, height=8)

