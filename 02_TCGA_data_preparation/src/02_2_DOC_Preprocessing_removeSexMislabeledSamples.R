library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled (female <-> male)
## 3. Remove outliers (Replace outlier signal with merginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by mean/median of each sample 
## 6. Removal of common germline CNVs
## 7. Remove mislabeled samples (NB <-> TP)

## Options
floor.val.frac <- 0.001 # Threshold value specified as double within [0,1].
#       CN values will be floored at floor.val.frac x (data mean)
zero.row.thresh <- 0.4 # Threshold value specifed as double within [0,1].  
#       Warning is issued for markers for which the corresponding row in data frame contains
#       more than zero.row.thresh x (number of samples) zeros. Default value for exome: 0.4
zero.col.thresh <- 0.05 # Threshold value specified as double within [0,1].
#       Samples whose corresponding column contains more than zero.col.thresh x (number of markers) zeros are removed.
#       Default value for exome: 0.05.
scale.method <- 'median' # Function used to scale each column. Can be mean or median.


## 2. Check if there are samples that have been mislabeled (female <-> male)
sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

dat.sample.probe.flt <- readRDS(file=here('02_TCGA_data_preparation/tmp/02_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'dat.sample.probe.flt.rds'))

doc.xy.g <- dat.sample.probe.flt[grepl('X|Y', rownames(dat.sample.probe.flt)), ] %>%
  select(ends_with('NB') | ends_with('NBC') | ends_with('NT')) %>%
  rownames_to_column('probe') %>%
  separate(data=., col=probe, into=c('Chr', 'pos'), sep=':') %>%
  separate(data=., col=pos, into=c('Start', 'Stop'), sep='-') %>%
  gather(key='SampleID', value='signal', -c('Chr', 'Start', 'Stop')) %>%
  left_join(sif, by='SampleID')

g <- ggplot(doc.xy.g %>% filter(gender=='male'), aes(x=Chr, y=signal, col=gender)) +
  # geom_boxplot() +
  # ggbeeswarm::geom_beeswarm() +
  geom_violin() +
  facet_grid(type~ID, scales='free_y')
show(g)

# xy.summary <- doc.xy.g %>%
#   group_by(ID, TYPE, Chr, gender, SITE) %>%
#   summarize(max=max(signal), min=min(signal), median=median(signal)) %>%
#   ungroup() %>%
#   gather(key='value', value='signal', -c('ID', 'TYPE', 'Chr', 'gender', 'SITE'))

# g <- ggplot(xy.summary %>% filter(gender=='female') %>% filter(ID %in% c('CESC.C5.A0TN', 'CESC.C5.A1BF')), aes(x=Chr, y=signal)) +
#   geom_point(aes(col=Chr, shape=value)) +
#   facet_grid(TYPE~ID, scales='free_y', space='free_y')
# show(g)

xy <- doc.xy.g %>%
  group_by(SampleID, ID, type, Chr, gender, project) %>%
  summarize(median=median(signal)) %>%
  ungroup() %>%
  pivot_wider(names_from='Chr', values_from='median') %>%
  mutate(y.x.ratio=Y/X) %>%
  left_join(qc.df %>% select(SampleID, TCGA.ID), by='SampleID')

## Supplementary data from Sadagopan 2022 Cell Systems
sadagopan.supp1.file <- here('data/XIST_Sadagopan2022CellSystems', '1-s2.0-S2405471222004033-mmc2.xlsx')
sadagopan.supp1 <- readxl::read_xlsx(sadagopan.supp1.file, skip=2) %>%
  as.data.frame() %>%
  setNames(gsub(' ', '.', colnames(.))) %>%
  mutate(Barcode=gsub('-', '.', Barcode))
sadagopan.male.suspicious.pts <- sadagopan.supp1 %>%
  filter(TCGA.Annotated.Sex=='MALE' & Determined.Sex=='FEMALE') %>%
  pull(Barcode) %>%
  sub('\\.[0-9]{2}$', '', .) %>%
  unique()
  
Klinefelter.ID <- sadagopan.supp1 %>%
  filter(Determined.Sex=='KLINEFELTER') %>%
  pull(Barcode) %>%
  sub('\\.[0-9]{2}$', '', .) %>%
  unique()

## Supplementary data from Meifang's paper
meifang.supp.file <- here('../20221114_ZSTangent_sex_TCGA_WES/data/Qi2022', 'media-2.xlsx')
meifang.mY <- readxl::read_xlsx(meifang.supp.file, sheet='Supp_Table_1') %>%
  select(-12)
meifang.fX <- readxl::read_xlsx(meifang.supp.file, sheet='Supp_Table_2')

meifang.mY.suspicious.pts <- meifang.mY %>%
  filter(Flag=='sex_discordant_or_abmigious_copy_number_sex_chromsome') %>%
  pull(case_id) %>%
  gsub('-', '.', .)
meifang.fX.suspicious.pts <- meifang.fX %>%
  filter(Flag=='sex_discordant_or_abmigious_copy_number_sex_chromsome') %>%
  pull(case_id) %>%
  gsub('-', '.', .)

kei.female.suspicious.pts <- xy %>%
  filter(gender=='female' & Y > 50) %>%
  pull(TCGA.ID) %>%
  unique()

kei.male.suspicious.pts <- xy %>%
  filter(gender=='male' & Y < 5) %>%
  pull(TCGA.ID) %>%
  unique()

g <- ggplot(xy, aes(x=X, y=Y)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=2, intercept=0), linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=0.5, intercept=0), linetype='dashed') +
  geom_point(aes(col=project), size=2) +
  lims(x=c(0, NA), y=c(0, NA)) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% kei.female.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% kei.male.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  # ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% sadagopan.male.suspicious.pts), aes(label=ID), col='blue', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% Klinefelter.ID), aes(label=ID), col='darkgreen', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% meifang.mY.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% meifang.fX.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  labs(x='X (median of depth)', y='Y (median of depth)') +
  facet_grid(type~gender) +
  theme_bw(base_size=20)
show(g)
ggsave(g, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'XY_signal.png'), dpi=100, width=16, height=12)

g <- ggplot(xy %>% filter(gender=='male' & type %in% c('NB', 'NT')), aes(x=X, y=Y)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=2, intercept=0), linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=0.5, intercept=0), linetype='dashed') +
  geom_point(aes(col=project), size=1, show.legend=FALSE) +
  lims(x=c(0, NA), y=c(0, NA)) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% kei.female.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% kei.male.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  # ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% sadagopan.male.suspicious.pts), aes(label=ID), col='blue', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% Klinefelter.ID), aes(label=ID), col='darkgreen', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% meifang.mY.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% meifang.fX.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  labs(x='X (median of depth)', y='Y (median of depth)') +
  facet_grid(type~project) +
  theme_bw(base_size=20)
show(g)
ggsave(g, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'XY_signal_byTumorType.png'), dpi=100, width=40, height=5)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

VennDiagram::venn.diagram(list(Kei=kei.male.suspicious.pts, Meifang=meifang.mY.suspicious.pts, Sadagopan=sadagopan.male.suspicious.pts),
  filename=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'Suspicious_male.tiff'),
  height=1000, width=1000,
  fill=c('blue', 'blue', 'blue'),
  alpha=c(0.5, 0.5, 0.5),
  lwd=0,
  # cat.col=c('green3', 'blue'),
  cat.dist=c(0.1, 0.1, 0.1),
  # cat.pos=c(45, 315), # for flipping
  # ext.pos=c(180, 180), # for flipping
  # rotation.degree=180, # for flipping
  margin=0.2)

VennDiagram::venn.diagram(list(Kei=kei.female.suspicious.pts, Meifang=meifang.fX.suspicious.pts),
  filename=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'Suspicious_female.tiff'),
  height=1000, width=1000,
  fill=c('red', 'red'),
  alpha=c(0.5, 0.5),
  lwd=0,
  # cat.col=c('green3', 'blue'),
  cat.dist=c(0.1, 0.1),
  # cat.pos=c(45, 315), # for flipping
  # ext.pos=c(180, 180), # for flipping
  # rotation.degree=180, # for flipping
  margin=0.2)

g <- ggplot(doc.xy.g %>% filter(SampleID %in% suspicious.samples), aes(x=Chr, y=signal, col=gender)) +
  # geom_boxplot() +
  # ggbeeswarm::geom_beeswarm() +
  geom_violin() +
  facet_grid(type~ID, scales='free_y')
show(g)

# g <- ggplot(doc.xy.g %>% filter(SampleID %in% suspicious.samples), aes(x=type, y=signal)) +
#   geom_boxplot(aes(fill=Chr)) +
#   # ggbeeswarm::geom_beeswarm() +
#   # geom_violin(aes(fill=Chr)) +
#   ggrepel::geom_text_repel(data=. %>% filter(SampleID %in% suspicious.samples) %>% group_by(ID, type, Chr) %>% summarize(signal=median(signal)), aes(x=type, label=signal, col=Chr), nudge_x=0.5, box.padding=1) +
#   labs(y='Depth') +
#   facet_grid(gender~ID, scales='free_y') +
#   # facet_wrap(~ID, nrow=2, scales='free') +
#   theme_bw(base_size=20) +
#   theme(axis.title.x = element_blank())
# show(g)
# ggsave(g, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'XY_signal_suspicious_pts.png'), dpi=100, width=24, height=12)

# g <- ggplot(doc.xy.g %>% filter(SampleID %in% suspicious.samples) %>% filter(gender=='female'), aes(x=type, y=signal)) +
#   geom_boxplot(aes(fill=Chr)) +
#   ggrepel::geom_text_repel(data=. %>% filter(SampleID %in% suspicious.samples) %>% group_by(ID, type, Chr) %>% summarize(signal=median(signal)), aes(x=type, label=signal, col=Chr), nudge_x=0.5, box.padding=1) +
#   labs(y='Depth') +
#   facet_wrap(~ID, nrow=1, scales='free') +
#   theme_bw(base_size=20) +
#   theme(axis.title.x = element_blank())
# show(g)
# ggsave(g, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'XY_signal_suspicious_pts_female.png'), dpi=100, width=12, height=6)

# g <- ggplot(doc.xy.g %>% filter(SampleID %in% suspicious.samples) %>% filter(gender=='male'), aes(x=type, y=signal)) +
#   geom_boxplot(aes(fill=Chr)) +
#   ggrepel::geom_text_repel(data=. %>% filter(SampleID %in% suspicious.samples) %>% group_by(ID, type, Chr) %>% summarize(signal=median(signal)), aes(x=type, label=signal, col=Chr), nudge_x=0.5, box.padding=1) +
#   labs(y='Depth') +
#   facet_wrap(~ID, nrow=1, scales='free') +
#   theme_bw(base_size=20) +
#   theme(axis.title.x = element_blank())
# show(g)
# ggsave(g, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'XY_signal_suspicious_pts_male.png'), dpi=100, width=24, height=12)

qc.df <- readRDS(file=here('output/07_1_preProcessing_1_removeBadProbesAndBadSamples', 'qc.df.RData'))
qc.df2 <- qc.df %>%
  mutate(suspicious.gender.kei=ifelse(TCGA.ID %in% c(kei.female.suspicious.pts, kei.male.suspicious.pts), TRUE, FALSE)) %>%
  mutate(suspicious.gender.meifang=ifelse(TCGA.ID %in% c(meifang.fX.suspicious.pts, meifang.mY.suspicious.pts), TRUE, FALSE)) %>%
  mutate(suspicious.gender.sadagopan=ifelse(TCGA.ID %in% sadagopan.male.suspicious.pts, TRUE, FALSE)) %>%
  mutate(suspicious.klinefelter=ifelse(TCGA.ID %in% Klinefelter.ID, TRUE, FALSE))
saveRDS(qc.df2, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'qc.df2.RData'), compress=FALSE)
# qc.df2 <- readRDS(file=here('output/02_preProcessing_2_removeSexMislabeledSamples', 'qc.df2.RData'))


## Select samples to be analyzed
samples.to.be.analyzed <- qc.df2 %>%
  filter(too.many.zeros==FALSE) %>%
  filter(suspicious.gender.kei==FALSE) %>%
  filter(suspicious.gender.meifang==FALSE) %>%
  filter(suspicious.gender.sadagopan==FALSE) %>%
  filter(suspicious.klinefelter==FALSE) %>%
  pull(SampleID)

types <- sif$type %>% unique()
for (i in 1:length(types)) {
  summary.i <- sif %>%
    filter(SampleID %in% samples.to.be.analyzed) %>%
    filter(type==types[i]) %>%
    select(ID, project, gender) %>%
    distinct() %>%
    select(project, gender) %>%
    table() %>%
    as.data.frame.matrix() %>%
    mutate(total=rowSums(.)) %>%
    rbind(total=colSums(.))
  print(types[i])
  print(summary.i)
}

dat.flt <- dat.sample.probe.flt %>% select(samples.to.be.analyzed)
saveRDS(dat.flt, file=here('output/07_2_preProcessing_2_removeSexMislabeledSamples', 'dat.flt.RData'), compress=FALSE)
