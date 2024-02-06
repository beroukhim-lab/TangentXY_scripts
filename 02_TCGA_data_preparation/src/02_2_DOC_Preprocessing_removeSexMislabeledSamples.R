library(tidyverse)
library(here)

## Preprocessing for Tangent
## 1. Remove samples (columns) with too many 0s and probes (rows) with too may 0s
## 2. Check if there are samples that have been mislabeled (female <-> male)
## 3. Remove outliers (Replace outlier signal with merginal median)
## 4.1. Replace 0 and small values with floor values to avoid -Inf in log2 transformation
## 4.2. Log2 transformation
## 5. Scaling by median of each sample 
## 6. Removal of common germline CNVs


## 2. Check if there are samples that have been mislabeled (female <-> male)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))
qc.df <- readRDS(here('02_TCGA_data_preparation/tmp/02_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'qc.df.rds'))

dat.sample.probe.flt <- readRDS(file=here('02_TCGA_data_preparation/tmp/02_1_DOC_Preprocessing_removeBadSamplesAndBadProbes', 'dat.sample.probe.flt.rds'))

xy.summary <- doc.xy.g <- dat.sample.probe.flt[grepl('X|Y', rownames(dat.sample.probe.flt)), ] %>%
  select(ends_with('NB') | ends_with('NBC') | ends_with('NT')) %>%
  rownames_to_column('probe') %>%
  separate(data=., col=probe, into=c('Chr', 'pos'), sep=':') %>%
  pivot_longer(names_to='SampleID', values_to='signal', cols=-c('Chr', 'pos')) %>%
  group_by(SampleID, Chr) %>%
  summarize(median=median(signal)) %>%
  ungroup() %>%
  pivot_wider(names_from='Chr', values_from='median') %>%
  mutate(y.x.ratio=Y/X) %>%
  left_join(sif, by='SampleID')

## Pick up likely mislabeled samples from literatures
## Supplementary data from Sadagopan et al., 2022 Cell Systems
sadagopan.supp.file <- here('02_TCGA_data_preparation/data/XIST_Sadagopan2022CellSystems', 'mmc2.xlsx') # Downloaded from Sadagopan et al., 2022 Cell Systems
sadagopan.supp <- readxl::read_xlsx(sadagopan.supp1.file, skip=2) %>%
  as.data.frame() %>%
  setNames(gsub(' ', '.', colnames(.))) %>%
  mutate(Barcode=gsub('-', '.', Barcode))
sadagopan.male.suspicious.pts <- sadagopan.supp %>%
  filter(TCGA.Annotated.Sex=='MALE' & Determined.Sex=='FEMALE') %>%
  pull(Barcode) %>%
  sub('\\.[0-9]{2}$', '', .) %>%
  unique()
  
Klinefelter.ID <- sadagopan.supp %>%
  filter(Determined.Sex=='KLINEFELTER') %>%
  pull(Barcode) %>%
  sub('\\.[0-9]{2}$', '', .) %>%
  unique()

## Supplementary data from Qi et al., 2023 Cell
qi.supp.file <- here('02_TCGA_data_preparation/data/LOY_Qi2023Cell', '1-s2.0-S0092867423006463-mmc1.xlsx') # Downloaded from Qi et al., 2023 Cell
qi.mY <- readxl::read_xlsx(qi.supp.file, sheet='Table S1A') %>%
  select(-13)
qi.fX <- readxl::read_xlsx(qi.supp.file, sheet='Table S1B')


qi.mY.suspicious.pts <- qi.mY %>%
  filter(Flag=='sex_discordant_or_ambiguous_copy_number_sex_chromsome') %>%
  pull(case_id) %>%
  gsub('-', '.', .)
qi.fX.suspicious.pts <- qi.fX %>%
  filter(Flag=='sex_discordant_or_ambiguous_copy_number_sex_chromsome') %>%
  pull(case_id) %>%
  gsub('-', '.', .)

## Pick up likely mislabeled samples based on my own criteria
enomoto.female.suspicious.pts <- xy.summary %>%
  filter(gender=='female' & Y > 50) %>%
  pull(TCGA.ID) %>%
  unique()

enomoto.male.suspicious.pts <- xy.summary %>%
  filter(gender=='male' & Y < 5) %>%
  pull(TCGA.ID) %>%
  unique()

# futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

## See how many suspicious samples are overlapped among studies
VennDiagram::venn.diagram(list(Enomoto=own.male.suspicious.pts, Qi=qi.mY.suspicious.pts, Sadagopan=sadagopan.male.suspicious.pts),
  filename=here('02_TCGA_data_preparation/output/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'Suspicious_male.png'),
  height=1000, width=1000,
  fill=c('blue', 'blue', 'blue'),
  alpha=c(0.5, 0.5, 0.5),
  lwd=0,
  cat.dist=c(0.1, 0.1, 0.1),
  margin=0.2)

VennDiagram::venn.diagram(list(Enomoto=own.female.suspicious.pts, Qi=qi.fX.suspicious.pts),
  filename=here('02_TCGA_data_preparation/output/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'Suspicious_female.png'),
  height=1000, width=1000,
  fill=c('red', 'red'),
  alpha=c(0.5, 0.5),
  lwd=0,
  cat.dist=c(0.1, 0.1),
  margin=0.2)

g <- ggplot(xy.summary, aes(x=X, y=Y)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=2, intercept=0), linetype='dashed') +
  geom_abline(data=. %>% filter(gender=='male'), aes(slope=0.5, intercept=0), linetype='dashed') +
  geom_point(aes(col=project), size=2) +
  lims(x=c(0, NA), y=c(0, NA)) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% enomoto.female.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% enomoto.male.suspicious.pts), aes(label=ID), col='black', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  # ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% sadagopan.male.suspicious.pts), aes(label=ID), col='blue', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% Klinefelter.ID), aes(label=ID), col='darkgreen', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% qi.mY.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  ggrepel::geom_text_repel(data=. %>% filter(TCGA.ID %in% qi.fX.suspicious.pts), aes(label=ID), col='pink', box.padding=1, position=position_jitter(), xlim=c(150, NA), max.overlaps=Inf, show.legend=FALSE) +
  labs(x='X (median of depth)', y='Y (median of depth)') +
  facet_grid(type~gender) +
  theme_bw(base_size=20)
ggsave(g, file=here('02_TCGA_data_preparation/output/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'XY_signal.png'), dpi=100, width=16, height=12)

qc.df2 <- qc.df %>%
  mutate(suspicious.gender.enomoto=ifelse(TCGA.ID %in% c(enomoto.female.suspicious.pts, enomoto.male.suspicious.pts), TRUE, FALSE)) %>%
  mutate(suspicious.gender.qi=ifelse(TCGA.ID %in% c(qi.fX.suspicious.pts, meifang.mY.suspicious.pts), TRUE, FALSE)) %>%
  mutate(suspicious.gender.sadagopan=ifelse(TCGA.ID %in% sadagopan.male.suspicious.pts, TRUE, FALSE)) %>%
  mutate(suspicious.klinefelter=ifelse(TCGA.ID %in% Klinefelter.ID, TRUE, FALSE)) %>%
  mutate(used.for.analysis=case_when(too.many.zeros==FALSE & suspicious.gender.enomoto==FALSE & suspicious.gender.qi==FALSE & suspicious.gender.sadagopan==FALSE & suspicious.klinefelter==FALSE ~ TRUE, TRUE ~ FALSE))
saveRDS(qc.df2, file=here('02_TCGA_data_preparation/tmp/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'qc.df2.rds'), compress=FALSE)


## Select samples to be analyzed
samples.to.be.analyzed <- qc.df2 %>%
  filter(used.for.analysis==TRUE) %>%
  pull(SampleID)

dat.flt <- dat.sample.probe.flt %>% select(samples.to.be.analyzed)
saveRDS(dat.flt, file=here('02_TCGA_data_preparation/tmp/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'dat.flt.rds'), compress=FALSE)
