library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))
qc <- readRDS(file=here('02_TCGA_data_preparation/output/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'qc.df2.rds'))
probes <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds')) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

absolute.file <- here('07_SCNAs_in_chrX_and_chrY/data', 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
absolute <- read.delim(absolute.file) %>%
  rename(barcode=sample) %>%
  separate(col=array, into=c('project', 'tss', 'participant', 'sample'), sep='-') %>%
  unite(col=TCGA.ID, c('project', 'tss', 'participant'), sep='.') %>%
  mutate(type=case_when(sample=='01' ~ 'TP',
                        sample=='02' ~ 'TR',
                        sample=='03' ~ 'TB',
                        sample=='05' ~ 'TAP',
                        sample=='06' ~ 'TM')) %>%
  left_join(sif, by=c('TCGA.ID', 'type'))

hg19 <- rCGH::hg19 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom))) %>%
  mutate(length=case_when(chrom=='Y' ~ 28800000, TRUE ~ length)) # Exclude heterochromatin region of chrY

arm.length <- hg19 %>%
  mutate(p=case_when(chrom=='Y' ~ centromerStart - 2600000, TRUE ~ centromerStart - 1)) %>% # Exclude PAR1 fo chrY
  mutate(q=length - centromerEnd) %>%
  rename(chr=chrom) %>%
  select(chr, p, q) %>%
  pivot_longer(names_to='arm', values_to='length', cols=c('p', 'q'))

arm.coverage <- probes %>%
  group_by(chr) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  group_by(chr, arm) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  mutate(start.pos=case_when(arm=='p' ~ start - 1, arm=='q' ~ as.numeric(centromerEnd))) %>%
  mutate(end.pos=case_when(arm=='p' ~ as.numeric(centromerStart), arm=='q' ~ end)) %>%
  mutate(covered.length=end.pos - start.pos) %>%
  left_join(arm.length, by=c('chr', 'arm')) %>%
  mutate(coverage=covered.length/length) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique())) %>%
  select(chr, arm, start.pos, end.pos, covered.length, length, coverage)

segment.smoothed.CNA.obj <- readRDS(file=here('03_TCGA_TangentXY/output/06_CBS', 'CBS_all.rds'))

segment <- segment.smoothed.CNA.obj %>%
  rename(SampleID=ID, chr=chrom) %>%
  left_join(probes %>% select(c('chr', 'start', 'end', 'centromerStart', 'centromerEnd')), by=c('chr', 'loc.end'='start')) %>%
  mutate(arm=case_when(end < centromerStart ~ 'p',
                      loc.start > centromerEnd ~ 'q',
                      loc.start < centromerStart & end > centromerStart & end < centromerEnd ~ 'p',
                      loc.start > centromerStart & loc.start < centromerEnd & end > centromerEnd ~ 'q',
                      loc.start < centromerStart & end > centromerEnd ~ 'overlap')) %>%
  mutate(centromer.overlap=case_when(arm=='overlap' ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

segment.centromer.overlap <- segment %>%
  filter(centromer.overlap==TRUE) %>%
  mutate(arm='q') %>%
  mutate(loc.start=centromerEnd + 1)

segment.arm <- segment %>%
  mutate(arm=case_when(centromer.overlap == TRUE ~ 'p', TRUE ~ arm)) %>%
  mutate(end=case_when(centromer.overlap == TRUE ~ centromerStart - 1, TRUE ~ end)) %>%
  bind_rows(segment.centromer.overlap) %>%
  left_join(qc, by='SampleID') %>%
  filter(Gender!='NA') %>%
  filter(used.for.analysis==TRUE) %>%
  filter(!(Gender=='Female' & chr=='Y')) %>%
  mutate(SampleID=factor(.$SampleID, levels=.$SampleID %>% unique())) %>%
  arrange(SampleID, chr, loc.start)
saveRDS(segment.arm, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'segment.arm.rds'), compress=FALSE)

## Assign Amp/Del to each segment
amp.del.offset <- 0.2
segment.amp.del <- segment.arm %>%
  mutate(alt.class=case_when(!(Gender=='Male' & chr %in% c('X', 'Y')) & seg.mean > amp.del.offset ~ 'Amp',
                              !(Gender=='Male' & chr %in% c('X', 'Y')) & seg.mean < (-1 * amp.del.offset) ~ 'Del',
                              Gender=='Male' & chr=='X' & seg.mean > -0.9 + amp.del.offset ~ 'Amp',
                              Gender=='Male' & chr=='X' & seg.mean < -0.9 - amp.del.offset ~ 'Del',
                              chr=='Y' & seg.mean > -1 + amp.del.offset ~ 'Amp',
                              chr=='Y' & seg.mean < -1 - amp.del.offset ~ 'Del',
                              TRUE ~ 'no.alt')) %>%
  mutate(seg.length = end - loc.start + 1) %>%
  left_join(arm.coverage %>% select(chr, arm, covered.length, length), by=c('chr', 'arm')) %>%
  mutate(seg.ratio = seg.length/covered.length)

alt.summary <- segment.amp.del %>%
  group_by(SampleID, chr, arm, alt.class) %>%
  summarize(alt.class.ratio=sum(seg.ratio))

arm.alt.thresh <- 0.5

arm.classifier <- function(df) {
  alt.summary <- df %>%
    group_by(alt.class) %>%
    summarize(alt.class.ratio=sum(seg.ratio))

  amp.ratio <- alt.summary %>% filter(alt.class=='Amp') %>% pull(alt.class.ratio)
  del.ratio <- alt.summary %>% filter(alt.class=='Del') %>% pull(alt.class.ratio)

  if (length(amp.ratio)!=1) {
    amp.ratio <- 0
  }
  if (length(del.ratio)!=1) {
    del.ratio <- 0
  }

  if (amp.ratio > arm.alt.thresh) {
    arm.class <- 'Amp'
  } else if (del.ratio > arm.alt.thresh) {
    arm.class <- 'Del'
  } else {
    arm.class <- 'No.Alt'
  }

  return(arm.class)
}

arm.amp.del <- segment.amp.del %>%
  group_by(SampleID, chr, arm) %>%
  nest() %>%
  mutate(arm.class=map_chr(data, arm.classifier)) %>%
  select(-data) %>%
  ungroup() %>%
  left_join(sif, by='SampleID')
saveRDS(arm.amp.del, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'arm.amp.del.rds'), compress=FALSE)

karyo.classifier <- function(df) {
  alt.class.unique <- df$arm.class %>% unique()
  alt.class.unified <- alt.class.unique %>% sort() %>% paste(collapse='_')
  if (length(alt.class.unique)==1) {
    karyo.detail <- NA
    if (alt.class.unified=='No.Alt') {
      karyo <- 'No.Alt'
    } else if (alt.class.unified=='Amp') {
      karyo <- 'Whole.Amp'
    } else if (alt.class.unified=='Del') {
      karyo <- 'Whole.Del'
    }
  } else {
    karyo.detail <- df %>%
      filter(arm.class!='No.Alt') %>%
      mutate(arm.karyo=paste(arm, arm.class, sep='_')) %>%
      pull(arm.karyo) %>%
      paste(collapse='&')
    if (alt.class.unified=='Amp_No.Alt') {
      karyo <- 'Arm.Amp'
    } else if (alt.class.unified=='Del_No.Alt') {
      karyo <- 'Arm.Del'
    } else if (alt.class.unified=='Amp_Del') {
      karyo <- 'Amp.Del'
    }
  }
  karyo.df <- data.frame(karyo.class=karyo, karyo.detail=karyo.detail)
  return(karyo.df)
}

sample.amp.del <- arm.amp.del %>%
  group_by(SampleID, chr) %>%
  nest() %>%
  mutate(karyo=map_df(data, karyo.classifier)) %>%
  select(-data) %>%
  unnest(cols=c(karyo)) %>%
  ungroup() %>%
  left_join(sif, by='SampleID') %>%
  mutate(chr.type=case_when(!chr %in% c('X', 'Y') ~ 'autosome',
                            chr=='X' ~ 'chrX',
                            chr=='Y' ~ 'chrY')) %>%
  left_join(absolute %>% select(SampleID, purity, ploidy, Genome.doublings), by='SampleID') %>%
  as.data.frame()
saveRDS(sample.amp.del, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'sample.amp.del.rds'), compress=FALSE)

sample.amp.del.count <- sample.amp.del %>%
  group_by(project, chr, Gender, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(project, chr, Gender) %>%
  mutate(total=sum(n)) %>%
  ungroup() %>%
  mutate(fraction=n/total) %>%
  mutate(Gender=factor(.$Gender, levels=c('Female', 'Male'))) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Alt')))

polyploidy.threshold <- 2.5

chrx.project.order <-  sample.amp.del.count %>%
  filter(chr=='X') %>%
  filter(Gender=='Female') %>%
  filter(karyo.class=='No.Alt') %>%
  arrange(fraction) %>%
  pull(project) %>%
  append(c('PRAD', 'TGCT'))

g <- ggplot(sample.amp.del.count %>%
      filter(chr=='X') %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chrx.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('No.Alt'='gray', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Amp.Del'='purple', 'Arm.Del'='cyan', 'Whole.Del'='blue')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrX', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3a.png'), dpi=100, width=32, height=8)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3a.pdf'), width=32, height=8)

sample.amp.del.chrx <- sample.amp.del %>%
  filter(chr=='X') %>%
  mutate(ploidy.class=case_when(ploidy < polyploidy.threshold & Genome.doublings==0 ~ 'Diploid', ploidy >= polyploidy.threshold | Genome.doublings > 0 ~ 'Polyploid', TRUE ~ 'NA')) %>%
  group_by(Gender, ploidy.class, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, ploidy.class) %>%
  mutate(total=sum(n)) %>%
  mutate(ploidy.class.total=paste0(ploidy.class, ' (n=', total, ')')) %>%
  ungroup() %>%
  mutate(fraction=n/total) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Alt'))) %>%
  mutate(Gender=factor(.$Gender, levels=c('Female', 'Male'))) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid', 'NA'))) %>%
  arrange(Gender, ploidy.class) %>%
  mutate(ploidy.class.total=factor(.$ploidy.class.total, levels=unique(.$ploidy.class.total)))

g <- ggplot(sample.amp.del.chrx, aes(x=ploidy.class.total, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity') +
  scale_fill_manual(values=c('No.Alt'='gray', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Amp.Del'='purple', 'Arm.Del'='cyan', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  lemon::facet_rep_wrap(~Gender, nrow=1, scales='free_x', repeat.tick.labels=TRUE) +
  labs(title='ChrX', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'FigS4a.png'), dpi=100, width=9, height=6)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'FigS4a.pdf'), width=9, height=6)

## ChrY
chry.project.order <-  sample.amp.del.count %>%
  filter(chr=='Y') %>%
  filter(Gender=='Male') %>%
  filter(karyo.class=='No.Alt') %>%
  arrange(fraction) %>%
  pull(project)

g <- ggplot(sample.amp.del.count %>%
      filter(chr=='Y') %>%
      mutate(project.n=paste0(project, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chry.project.order)),
    aes(x=project.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('No.Alt'='gray', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Amp.Del'='purple', 'Arm.Del'='cyan', 'Whole.Del'='blue')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrY', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3b.png'), dpi=100, width=32, height=8)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3b.pdf'), width=32, height=8)

sample.amp.del.chry <- sample.amp.del %>%
  filter(chr=='Y') %>%
  mutate(ploidy.class=case_when(ploidy < polyploidy.threshold & Genome.doublings==0 ~ 'Diploid', ploidy >= polyploidy.threshold | Genome.doublings > 0 ~ 'Polyploid', TRUE ~ 'NA')) %>%
  group_by(Gender, ploidy.class, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, ploidy.class) %>%
  mutate(total=sum(n)) %>%
  mutate(ploidy.class.total=paste0(ploidy.class, ' (n=', total, ')')) %>%
  ungroup() %>%
  mutate(fraction=n/total) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Alt'))) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid', 'NA'))) %>%
  arrange(ploidy.class) %>%
  mutate(ploidy.class.total=factor(.$ploidy.class.total, levels=unique(.$ploidy.class.total)))

g <- ggplot(sample.amp.del.chry, aes(x=ploidy.class.total, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity') +
  scale_fill_manual(values=c('No.Alt'='gray', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Amp.Del'='purple', 'Arm.Del'='cyan', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  lemon::facet_rep_wrap(~Gender, nrow=1, scales='free_x', repeat.tick.labels=TRUE) +
  labs(title='ChrY', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'FigS4b.png'), dpi=100, width=6, height=6)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'FigS4b.pdf'), width=6, height=6)
