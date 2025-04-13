library(tidyverse)
library(here)

sif <- readRDS(file=here('02_TCGA_data_preparation/output/00_format_sif', 'sif.rds'))
qc <- readRDS(file=here('02_TCGA_data_preparation/output/02_2_DOC_Preprocessing_removeSexMislabeledSamples', 'qc.df2.rds'))
probes <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds')) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

polyploidy.threshold <- 2.5

absolute.file <- here('07_SCNAs_in_chrX_and_chrY/data', 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
absolute <- read.delim(absolute.file) %>%
  rename(barcode=sample) %>%
  mutate(D=(ploidy*purity) + 2*(1-purity)) %>%
  mutate(ploidy.class=case_when(ploidy < polyploidy.threshold & Genome.doublings==0 ~ 'Diploid', ploidy >= polyploidy.threshold | Genome.doublings > 0 ~ 'Polyploid', TRUE ~ 'NA')) %>%
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
  mutate(p=case_when(chrom=='Y' ~ centromerStart - 2600000, TRUE ~ centromerStart - 1)) %>% # Exclude PAR1 of chrY
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

  if (amp.ratio >= arm.alt.thresh) {
    arm.class <- 'Amp'
  } else if (del.ratio >= arm.alt.thresh) {
    arm.class <- 'Del'
  } else {
    arm.class <- 'No.Arm-level.Alt'
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
    if (alt.class.unified=='No.Arm-level.Alt') {
      karyo <- 'No.Arm-level.Alt'
    } else if (alt.class.unified=='Amp') {
      karyo <- 'Whole.Amp'
    } else if (alt.class.unified=='Del') {
      karyo <- 'Whole.Del'
    }
  } else {
    karyo.detail <- df %>%
      filter(arm.class!='No.Arm-level.Alt') %>%
      mutate(arm.karyo=paste(arm, arm.class, sep='_')) %>%
      pull(arm.karyo) %>%
      paste(collapse='&')
    if (alt.class.unified=='Amp_No.Arm-level.Alt') {
      karyo <- 'Arm.Amp'
    } else if (alt.class.unified=='Del_No.Arm-level.Alt') {
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
  left_join(absolute %>% select(SampleID, purity, ploidy, Genome.doublings, ploidy.class), by='SampleID') %>%
  mutate(ploidy.class=case_when(is.na(ploidy.class) ~ 'NA', TRUE ~ ploidy.class)) %>%
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
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Arm-level.Alt')))

polyploidy.threshold <- 2.5

chrx.project.order <-  sample.amp.del.count %>%
  filter(chr=='X') %>%
  filter(Gender=='Female') %>%
  filter(karyo.class=='No.Arm-level.Alt') %>%
  arrange(fraction) %>%
  pull(project) %>%
  append(c('PRAD', 'TGCT'))

g <- ggplot(sample.amp.del.count %>%
      filter(chr=='X') %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chrx.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('No.Arm-level.Alt'='gray', 'Whole.Amp'='#D7191C', 'Arm.Amp'='#FDAE61', 'Amp.Del'='#FFFFBF', 'Arm.Del'='#ABD9E9', 'Whole.Del'='#2C7BB6')) +
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
  group_by(Gender, ploidy.class, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, ploidy.class) %>%
  mutate(total=sum(n)) %>%
  mutate(ploidy.class.total=paste0(ploidy.class, ' (n=', total, ')')) %>%
  ungroup() %>%
  mutate(fraction=n/total) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Arm-level.Alt'))) %>%
  mutate(Gender=factor(.$Gender, levels=c('Female', 'Male'))) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid', 'NA'))) %>%
  arrange(Gender, ploidy.class) %>%
  mutate(ploidy.class.total=factor(.$ploidy.class.total, levels=unique(.$ploidy.class.total)))

g <- ggplot(sample.amp.del.chrx, aes(x=ploidy.class.total, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity') +
  scale_fill_manual(values=c('No.Arm-level.Alt'='gray', 'Whole.Amp'='#D7191C', 'Arm.Amp'='#FDAE61', 'Amp.Del'='#FFFFBF', 'Arm.Del'='#ABD9E9', 'Whole.Del'='#2C7BB6')) +
  scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  lemon::facet_rep_wrap(~Gender, nrow=1, scales='free_x', repeat.tick.labels=TRUE) +
  labs(title='ChrX', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3c.png'), dpi=100, width=9, height=6)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3c.pdf'), width=9, height=6)

## ChrY
chry.project.order <-  sample.amp.del.count %>%
  filter(chr=='Y') %>%
  filter(Gender=='Male') %>%
  filter(karyo.class!='No.Arm-level.Alt') %>%
  group_by(project) %>%
  summarize(fraction=sum(fraction)) %>%
  arrange(desc(fraction)) %>%
  pull(project)

g <- ggplot(sample.amp.del.count %>%
      filter(chr=='Y') %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chry.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('No.Arm-level.Alt'='gray', 'Whole.Amp'='#D7191C', 'Arm.Amp'='#FDAE61', 'Amp.Del'='#FFFFBF', 'Arm.Del'='#ABD9E9', 'Whole.Del'='#2C7BB6')) +
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
  group_by(Gender, ploidy.class, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, ploidy.class) %>%
  mutate(total=sum(n)) %>%
  mutate(ploidy.class.total=paste0(ploidy.class, ' (n=', total, ')')) %>%
  ungroup() %>%
  mutate(fraction=n/total) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Whole.Amp', 'Arm.Amp', 'Amp.Del', 'Arm.Del', 'Whole.Del', 'No.Arm-level.Alt'))) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid', 'NA'))) %>%
  arrange(ploidy.class) %>%
  mutate(ploidy.class.total=factor(.$ploidy.class.total, levels=unique(.$ploidy.class.total)))

g <- ggplot(sample.amp.del.chry, aes(x=ploidy.class.total, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity') +
  scale_fill_manual(values=c('No.Arm-level.Alt'='gray', 'Whole.Amp'='#D7191C', 'Arm.Amp'='#FDAE61', 'Amp.Del'='#FFFFBF', 'Arm.Del'='#ABD9E9', 'Whole.Del'='#2C7BB6')) +
  scale_y_continuous(limits=c(0, 1.0), breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  lemon::facet_rep_wrap(~Gender, nrow=1, scales='free_x', repeat.tick.labels=TRUE) +
  labs(title='ChrY', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3d.png'), dpi=100, width=6, height=6)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3d.pdf'), width=6, height=6)








## Shahab's method (See "male-bias-mutation-chrY-loss-correlation.R")
# Reading in TCGA tangent-normalized chrY copy number data
input_cn <- readRDS(file=here('../../../shahab/tangent-Y/Y_shifted_sexMatchedTangentOnMale.RData'))
Tn.male.normalized <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn_sexMatchedTangentOnMale.rds'))

# segment.smoothed.CNA.obj <- readRDS(file=here('03_TCGA_TangentXY/output/06_CBS', 'CBS_all.rds'))

# Formatting matrix cn data into a long data frame
# y_cn <- input_cn %>%
y_cn <- Tn.male.normalized %>%
  as.data.frame() %>%
  bind_cols(probes %>% select(chr,start,end)) %>%
  filter(chr=='Y') %>%
  pivot_longer(names_to = 'SampleID', values_to = 'signal', 
               cols = colnames(Tn.male.normalized))

# Using median signal as the chrY CN for each tumor
median_y_cn <- y_cn %>%
  group_by(SampleID) %>%
  summarize(median=median(signal))

# Retaining samples with known male gender
# Calculating chrY integer CN using absolute purity/ploidy estimates
integer_y_cn <- median_y_cn %>% 
  left_join(absolute, by='SampleID') %>%
  mutate(CN=(D*2^(median)-1+purity)/purity) %>%
  mutate(CN_status=ifelse(CN<0.5, "CN<0.5", "CN>=0.5")) %>%
  filter(!is.na(CN_status))

g <- ggplot(integer_y_cn, aes(x=CN)) +
  geom_histogram(binwidth=0.02) +
  geom_vline(xintercept=c(0, 1, 2, 3, 4), linetype='dashed') +
  facet_wrap(~Gender) +
  labs(title='ChrY integer CN') +
  theme_bw(base_size=20)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'ChrY_IntegerCN_distribution.png'), dpi=100, width=10, height=6)

integer_y_cn_class <- integer_y_cn %>%
  mutate(karyo.class=case_when(ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                ploidy.class=='Diploid' & CN >= 0.2 & CN < 0.8 ~ 'Subclonal Loss',
                                ploidy.class=='Diploid' & CN >= 0.8 & CN < 1.2 ~ 'Neutral',
                                ploidy.class=='Diploid' & CN >= 1.2 ~ 'Amp',
                                ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 2.2 ~ 'Amp',
                                ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 4.2 ~ 'Amp',
                                ))

integer_y_cn_class_summary <- integer_y_cn_class %>%
  filter(!is.na(karyo.class)) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Amp', 'Subclonal Loss', 'Clonal Loss', 'Neutral'))) %>%
  group_by(Gender, project, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, project) %>%
  mutate(total=sum(n)) %>%
  mutate(fraction=n/total) %>%
  ungroup()

chry.project.order <- integer_y_cn_class_summary %>%
  select(project, karyo.class, fraction) %>%
  distinct() %>%
  pivot_wider(names_from=karyo.class, values_from=fraction) %>%
  as.data.frame() %>%
  replace(is.na(.), 0) %>%
  arrange(Neutral, desc(`Clonal Loss`)) %>%
  pull(project)

g <- ggplot(integer_y_cn_class_summary %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chry.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('Neutral'='gray', 'Amp'='#D7191C', 'Subclonal Loss'='#ABD9E9', 'Clonal Loss'='#2C7BB6')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrY', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3b_new.png'), dpi=100, width=32, height=8)

## Arm-level
median_y_arm_cn <- y_cn %>%
  left_join(probes, by=c('chr', 'start', 'end')) %>%
  group_by(SampleID, arm) %>%
  summarize(median=median(signal))

integer_y_arm_cn <- median_y_arm_cn %>% 
  left_join(absolute, by='SampleID') %>%
  mutate(CN=(D*2^(median)-1+purity)/purity) %>%
  filter(!is.na(CN))

g <- ggplot(integer_y_arm_cn, aes(x=CN)) +
  geom_histogram(binwidth=0.02) +
  geom_vline(xintercept=c(0, 1, 2, 3, 4), linetype='dashed') +
  facet_grid(arm~Gender) +
  labs(title='ChrY integer CN of each arm') +
  theme_bw(base_size=20)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'ChrY_IntegerCN_distribution_arm.png'), dpi=100, width=10, height=6)

integer_y_arm_cn_class <- integer_y_arm_cn %>%
  mutate(karyo.class.arm=case_when(ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                  ploidy.class=='Diploid' & CN >= 0.2 & CN < 0.8 ~ 'Subclonal Loss',
                                  ploidy.class=='Diploid' & CN >= 0.8 & CN < 1.2 ~ 'Neutral',
                                  ploidy.class=='Diploid' & CN >= 1.2 ~ 'Amp',
                                  ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                  ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                  ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                  ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 2.2 ~ 'Amp',
                                  ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                  ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                  ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                  ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 4.2 ~ 'Amp',
                                  )) %>%
  select(SampleID, Gender, project, arm, karyo.class.arm) %>%
  pivot_wider(names_from='arm', values_from='karyo.class.arm') %>%
  as.data.frame() %>%
  mutate(p=factor(.$p, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp'))) %>%
  mutate(q=factor(.$q, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp')))

integer_y_arm_cn_class <- integer_y_arm_cn %>%
  mutate(buffer=ploidy * 0.5 * 0.2) %>%
  mutate(karyo.class.arm=case_when(CN < buffer ~ 'Clonal Loss',
                                  CN >= buffer & CN < ploidy * 0.5 - buffer ~ 'Subclonal Loss',
                                  CN >= ploidy * 0.5 - buffer & CN < ploidy * 0.5 + buffer ~ 'Neutral',
                                  CN >= ploidy * 0.5 + buffer ~ 'Amp'
                                )) %>%
  select(SampleID, Gender, project, arm, karyo.class.arm) %>%
  pivot_wider(names_from='arm', values_from='karyo.class.arm') %>%
  as.data.frame() %>%
  mutate(p=factor(.$p, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp'))) %>%
  mutate(q=factor(.$q, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp')))

integer_y_arm_cn_class_comb <- integer_y_arm_cn_class %>%
  filter(!is.na(p) & !is.na(q)) %>%
  mutate(karyo.class=case_when(p=='Clonal Loss' & q=='Clonal Loss' ~ 'Clonal Loss',
                                p=='Subclonal Loss' & q=='Subclonal Loss' ~ 'Subclonal Loss',
                                p=='Neutral' & q=='Neutral' ~ 'Neutral',
                                p=='Amp' & q=='Amp' ~ 'Amp',
                                (p=='Clonal Loss' & q=='Amp') | (p=='Amp' & q=='Clonal Loss') ~ 'Amp & Loss',
                                TRUE ~ 'Other')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Amp', 'Amp & Loss', 'Subclonal Loss', 'Clonal Loss', 'Other', 'Neutral')))

integer_y_arm_cn_class_comb_summary <- integer_y_arm_cn_class_comb %>%
  group_by(Gender, project, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, project) %>%
  mutate(total=sum(n)) %>%
  mutate(fraction=n/total) %>%
  ungroup()

chry.project.order <- integer_y_arm_cn_class_summary %>%
  select(project, karyo.class, fraction) %>%
  distinct() %>%
  pivot_wider(names_from=karyo.class, values_from=fraction) %>%
  as.data.frame() %>%
  replace(is.na(.), 0) %>%
  arrange(Neutral, desc(`Clonal Loss`)) %>%
  pull(project)

g <- ggplot(integer_y_arm_cn_class_summary %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chry.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('Amp'='#D7191C', 'Amp & Loss'='#FFFFBF', 'Subclonal Loss'='#ABD9E9', 'Clonal Loss'='#2C7BB6', 'Other'='purple', 'Neutral'='gray')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrY', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3b_new_arm.png'), dpi=100, width=32, height=8)




## ChrX
Tn <- readRDS(file=here('03_TCGA_TangentXY/output/03_TangentXY', 'Tn.rds'))

# Formatting matrix cn data into a long data frame
x_cn <- Tn %>%
  as.data.frame() %>%
  bind_cols(probes %>% select(chr,start,end)) %>%
  filter(chr=='X') %>%
  pivot_longer(names_to = 'SampleID', values_to = 'signal', 
               cols = colnames(Tn))

# Using median signal as the chrX CN for each tumor
median_x_cn <- x_cn %>%
  group_by(SampleID) %>%
  summarize(median=median(signal))

# Calculating chrX integer CN using absolute purity/ploidy estimates
integer_x_cn <- median_x_cn %>% 
  left_join(absolute, by='SampleID') %>%
  mutate(CN=case_when(Gender=='Female' ~ (D*2^(median)-2*(1-purity))/purity, Gender=='Male' ~ (D*2^(median)-1+purity)/purity)) %>%
  filter(!is.na(CN)) %>%
  filter(!is.na(Gender))

g <- ggplot(integer_x_cn, aes(x=CN)) +
  geom_histogram(binwidth=0.02) +
  geom_vline(xintercept=c(0, 1, 2, 3, 4), linetype='dashed') +
  facet_wrap(~Gender, nrow=2) +
  labs(title='ChrX integer CN') +
  theme_bw(base_size=20)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'ChrX_IntegerCN_distribution.png'), dpi=100, width=10, height=6)

integer_x_cn_class <- integer_x_cn %>%
  mutate(karyo.class=case_when(
                                Gender=='Female' & ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Female' & ploidy.class=='Diploid' & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                Gender=='Female' & ploidy.class=='Diploid' & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                Gender=='Female' & ploidy.class=='Diploid' & CN >= 2.2 ~ 'Amp',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 4.2 ~ 'Amp',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 7.8 ~ 'Subclonal Loss',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 7.8 & CN < 8.2 ~ 'Neutral',
                                Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 8.2 ~ 'Amp',
                                Gender=='Male' & ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Male' & ploidy.class=='Diploid' & CN >= 0.2 & CN < 0.8 ~ 'Subclonal Loss',
                                Gender=='Male' & ploidy.class=='Diploid' & CN >= 0.8 & CN < 1.2 ~ 'Neutral',
                                Gender=='Male' & ploidy.class=='Diploid' & CN >= 1.2 ~ 'Amp',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 2.2 ~ 'Amp',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 4.2 ~ 'Amp'
                                ))

integer_x_cn_class_summary <- integer_x_cn_class %>%
  filter(!is.na(karyo.class)) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Amp', 'Subclonal Loss', 'Clonal Loss', 'Neutral'))) %>%
  group_by(Gender, project, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, project) %>%
  mutate(total=sum(n)) %>%
  mutate(fraction=n/total) %>%
  ungroup()

chrx.project.order <- integer_x_cn_class_summary %>%
  filter(Gender=='Female') %>%
  select(project, karyo.class, fraction) %>%
  distinct() %>%
  pivot_wider(names_from=karyo.class, values_from=fraction) %>%
  as.data.frame() %>%
  replace(is.na(.), 0) %>%
  arrange(Neutral, desc(`Clonal Loss`)) %>%
  pull(project) %>%
  append(c('PRAD', 'TGCT'))

g <- ggplot(integer_x_cn_class_summary %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chrx.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('Neutral'='gray', 'Amp'='#D7191C', 'Subclonal Loss'='#ABD9E9', 'Clonal Loss'='#2C7BB6')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrX', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3a_new.png'), dpi=100, width=32, height=8)


## Arm-level
median_x_arm_cn <- x_cn %>%
  left_join(probes, by=c('chr', 'start', 'end')) %>%
  group_by(SampleID, arm) %>%
  summarize(median=median(signal))

integer_x_arm_cn <- median_x_arm_cn %>% 
  left_join(absolute, by='SampleID') %>%
  mutate(CN=case_when(Gender=='Female' ~ (D*2^(median)-2*(1-purity))/purity, Gender=='Male' ~ (D*2^(median)-1+purity)/purity)) %>%
  filter(!is.na(CN)) %>%
  filter(!is.na(Gender))

g <- ggplot(integer_x_arm_cn, aes(x=CN)) +
  geom_histogram(binwidth=0.02) +
  geom_vline(xintercept=c(0, 1, 2, 3, 4), linetype='dashed') +
  facet_grid(arm~Gender) +
  labs(title='ChrX integer CN of each arm') +
  theme_bw(base_size=20)
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'ChrX_IntegerCN_distribution_arm.png'), dpi=100, width=10, height=6)

integer_x_arm_cn_class <- integer_x_arm_cn %>%
  mutate(karyo.class.arm=case_when(
                                  Gender=='Female' & ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Female' & ploidy.class=='Diploid' & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                  Gender=='Female' & ploidy.class=='Diploid' & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                  Gender=='Female' & ploidy.class=='Diploid' & CN >= 2.2 ~ 'Amp',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 4.2 ~ 'Amp',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 7.8 ~ 'Subclonal Loss',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 7.8 & CN < 8.2 ~ 'Neutral',
                                  Gender=='Female' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 8.2 ~ 'Amp',
                                  Gender=='Male' & ploidy.class=='Diploid' & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Male' & ploidy.class=='Diploid' & CN >= 0.2 & CN < 0.8 ~ 'Subclonal Loss',
                                  Gender=='Male' & ploidy.class=='Diploid' & CN >= 0.8 & CN < 1.2 ~ 'Neutral',
                                  Gender=='Male' & ploidy.class=='Diploid' & CN >= 1.2 ~ 'Amp',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 0.2 & CN < 1.8 ~ 'Subclonal Loss',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 1.8 & CN < 2.2 ~ 'Neutral',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==1 & CN >= 2.2 ~ 'Amp',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN < 0.2 ~ 'Clonal Loss',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 0.2 & CN < 3.8 ~ 'Subclonal Loss',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 3.8 & CN < 4.2 ~ 'Neutral',
                                  Gender=='Male' & ploidy.class=='Polyploid' & Genome.doublings==2 & CN >= 4.2 ~ 'Amp'
                                )) %>%
  select(SampleID, Gender, project, arm, karyo.class.arm) %>%
  pivot_wider(names_from='arm', values_from='karyo.class.arm') %>%
  as.data.frame() %>%
  mutate(p=factor(.$p, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp'))) %>%
  mutate(q=factor(.$q, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp')))

integer_x_arm_cn_class <- integer_x_arm_cn %>%
  mutate(buffer=ploidy * 0.5 * 0.2) %>%
  mutate(karyo.class.arm=case_when(
                                  Gender=='Female' & CN < buffer ~ 'Clonal Loss',
                                  Gender=='Female' & CN >= buffer & CN < ploidy - buffer ~ 'Subclonal Loss',
                                  Gender=='Female' & CN >= ploidy - buffer & CN < ploidy + buffer ~ 'Neutral',
                                  Gender=='Female' & CN >= ploidy + buffer ~ 'Amp',
                                  Gender=='Male' & CN < buffer ~ 'Clonal Loss',
                                  Gender=='Male' & CN >= buffer & CN < ploidy * 0.5 - buffer ~ 'Subclonal Loss',
                                  Gender=='Male' & CN >= ploidy * 0.5 - buffer & CN < ploidy * 0.5 + buffer ~ 'Neutral',
                                  Gender=='Male' & CN >= ploidy * 0.5 + buffer ~ 'Amp'
                                )) %>%
  select(SampleID, Gender, project, arm, karyo.class.arm) %>%
  pivot_wider(names_from='arm', values_from='karyo.class.arm') %>%
  as.data.frame() %>%
  mutate(p=factor(.$p, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp'))) %>%
  mutate(q=factor(.$q, levels=c('Clonal Loss', 'Subclonal Loss', 'Neutral', 'Amp')))

integer_x_arm_cn_class_comb <- integer_x_arm_cn_class %>%
  filter(!is.na(p) & !is.na(q)) %>%
  mutate(karyo.class=case_when(p=='Clonal Loss' & q=='Clonal Loss' ~ 'Clonal Loss',
                                p=='Subclonal Loss' & q=='Subclonal Loss' ~ 'Subclonal Loss',
                                p=='Neutral' & q=='Neutral' ~ 'Neutral',
                                p=='Amp' & q=='Amp' ~ 'Amp',
                                (p=='Clonal Loss' & q=='Amp') | (p=='Amp' & q=='Clonal Loss') ~ 'Amp & Loss',
                                TRUE ~ 'Other')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('Amp', 'Amp & Loss', 'Subclonal Loss', 'Clonal Loss', 'Other', 'Neutral')))

integer_x_arm_cn_class_comb_summary <- integer_x_arm_cn_class_comb %>%
  group_by(Gender, project, karyo.class) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(Gender, project) %>%
  mutate(total=sum(n)) %>%
  mutate(fraction=n/total) %>%
  ungroup()

chrx.project.order <- integer_x_arm_cn_class_comb_summary %>%
  filter(Gender=='Female') %>%
  select(project, karyo.class, fraction) %>%
  distinct() %>%
  pivot_wider(names_from=karyo.class, values_from=fraction) %>%
  as.data.frame() %>%
  replace(is.na(.), 0) %>%
  arrange(Neutral, desc(`Clonal Loss`)) %>%
  pull(project) %>%
  append(c('PRAD', 'TGCT'))

g <- ggplot(integer_x_arm_cn_class_comb_summary %>%
      mutate(Gender.n=paste0(Gender, ' (n=', total, ')')) %>%
      mutate(project=factor(.$project, levels=chrx.project.order)),
    aes(x=Gender.n, y=fraction)) +
  geom_bar(aes(fill=karyo.class), stat='identity', position='fill') +
  scale_fill_manual(values=c('Amp'='#D7191C', 'Amp & Loss'='#FFFFBF', 'Subclonal Loss'='#ABD9E9', 'Clonal Loss'='#2C7BB6', 'Other'='purple', 'Neutral'='gray')) +
  scale_y_continuous(breaks=seq(0, 1.0, by=0.2), expand=c(0, 0)) +
  facet_wrap(~project, nrow=1, scales='free_x', strip.position='bottom') +
  labs(title='ChrX', y='Fraction of patients', fill='Alteration type') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'Fig3a_new_arm.png'), dpi=100, width=32, height=8)
