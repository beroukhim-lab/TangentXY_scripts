library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))

probes <- readRDS(file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'probes.rds')) %>%
  rename(Chr=chr) %>%
  rename(chr=chrom) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

hg38 <- rCGH::hg38 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom))) %>%
  mutate(length=case_when(chrom=='Y' ~ 26600000, TRUE ~ length)) # Exclude heterochromatin region of chrY (https://en.wikipedia.org/wiki/Y_chromosome)

arm.length <- hg38 %>%
  mutate(p=case_when(chrom=='Y' ~ centromerStart - 2781479 + 1500000, TRUE ~ centromerStart + 1500000)) %>% # Exclude PAR1 fo chrY (https://en.wikipedia.org/wiki/Pseudoautosomal_region)
  mutate(q=length - centromerEnd + 1500000) %>%
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
  mutate(start.pos=case_when(arm=='p' ~ start - 1, arm=='q' ~ as.numeric(centromerEnd) - 1500000)) %>%
  mutate(end.pos=case_when(arm=='p' ~ as.numeric(centromerStart) + 1500000, arm=='q' ~ end)) %>%
  mutate(covered.length=end.pos - start.pos) %>%
  left_join(arm.length, by=c('chr', 'arm')) %>%
  mutate(coverage=covered.length/length) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique())) %>%
  select(chr, arm, start.pos, end.pos, covered.length, length, coverage)

Mahmoud.supp.file <- here('05_CCLE_data_preparation/data/CCLE_Mahmoud2019Nature', '41586_2019_1186_MOESM4_ESM.xlsx')
annotations <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Cell Line Annotations') %>%
  mutate(DepMapID=sub('-', '.', depMapID))

datasets <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Datasets')

cbio.file <- here('10_chrX_SCNA_vs_gene_expression/data', 'ccle_broad_2019_clinical_data.tsv') # Downloaded from cBioPortal (https://www.cbioportal.org/study/clinicalData?id=ccle_broad_2019)
cbio <- read.delim(cbio.file) %>%
  setNames(gsub('\\.', '', colnames(.))) %>%
  mutate(DepMapID=gsub('-', '.', DepMapID)) %>%
  left_join(annotations %>% select(DepMapID, tcga_code), by='DepMapID')

segment.smoothed.CNA.obj <- readRDS(file=here('06_CCLE_TangentXY/output/07_CBS', 'CBS_TangentXY.rds'))

segment <- segment.smoothed.CNA.obj$output %>%
  rename(DepMapID=ID, Chr=chrom) %>%
  left_join(probes %>% select(c('Chr', 'chr', 'start', 'end', 'centromerStart', 'centromerEnd')), by=c('Chr', 'loc.end'='start')) %>%
  mutate(arm=case_when(end < centromerStart + 1500000 ~ 'p',
                      loc.start > centromerEnd - 1500000 ~ 'q',
                      loc.start < centromerStart + 1500000 & end > centromerEnd - 1500000 ~ 'overlap')) %>%
  mutate(centromer.overlap=case_when(arm=='overlap' ~ TRUE, TRUE ~ FALSE)) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

segment.centromer.overlap <- segment %>%
  filter(centromer.overlap==TRUE) %>%
  mutate(arm='q') %>%
  mutate(loc.start=centromerEnd - 1500000 + 1)

segment.arm <- segment %>%
  mutate(arm=case_when(centromer.overlap == TRUE ~ 'p', TRUE ~ arm)) %>%
  mutate(end=case_when(centromer.overlap == TRUE ~ centromerStart + 1500000, TRUE ~ end)) %>%
  bind_rows(segment.centromer.overlap) %>%
  left_join(sif %>% select(DepMap_ID, sex, lineage), by=c('DepMapID'='DepMap_ID')) %>%
  filter(!is.na(sex)) %>%
  filter(!(sex=='Female' & chr=='Y')) %>%
  mutate(DepMapID=factor(.$DepMapID, levels=.$DepMapID %>% unique())) %>%
  arrange(DepMapID, chr, loc.start)
saveRDS(segment.arm, file=here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'segment.arm.rds'), compress=FALSE)

## Assign Amp/Del to each segment
amp.del.offset <- 0.2
segment.amp.del <- segment.arm %>%
  mutate(alt.class=case_when(!(sex=='Male' & chr %in% c('X', 'Y')) & seg.mean > amp.del.offset ~ 'Amp',
                              !(sex=='Male' & chr %in% c('X', 'Y')) & seg.mean < (-1 * amp.del.offset) ~ 'Del',
                              sex=='Male' & chr=='X' & seg.mean > -0.9 + amp.del.offset ~ 'Amp',
                              sex=='Male' & chr=='X' & seg.mean < -0.9 - amp.del.offset ~ 'Del',
                              chr=='Y' & seg.mean > -1 + amp.del.offset ~ 'Amp',
                              chr=='Y' & seg.mean < -1 - amp.del.offset ~ 'Del',
                              TRUE ~ 'no.alt')) %>%
  mutate(seg.length = end - loc.start + 1) %>%
  left_join(arm.coverage %>% select(chr, arm, covered.length, length), by=c('chr', 'arm')) %>%
  mutate(seg.ratio = seg.length/covered.length)

alt.summary <- segment.amp.del %>%
  group_by(DepMapID, chr, arm, alt.class) %>%
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
  group_by(DepMapID, chr, arm) %>%
  nest() %>%
  mutate(arm.class=map_chr(data, arm.classifier)) %>%
  select(-data) %>%
  ungroup() %>%
  left_join(sif %>% select(DepMap_ID, lineage, sex), by=c('DepMapID'='DepMap_ID'))
saveRDS(arm.amp.del, file=here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'arm.amp.del.rds'), compress=FALSE)

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
  group_by(DepMapID, chr) %>%
  nest() %>%
  mutate(karyo=map_df(data, karyo.classifier)) %>%
  select(-data) %>%
  unnest(cols=c(karyo)) %>%
  ungroup() %>%
  left_join(sif %>% select(DepMap_ID, lineage, sex), by=c('DepMapID'='DepMap_ID')) %>%
  mutate(chr.type=case_when(!chr %in% c('X', 'Y') ~ 'autosome',
                            chr=='X' ~ 'chrX',
                            chr=='Y' ~ 'chrY')) %>%
  left_join(cbio %>% select(DepMapID, Purity, Ploidy, GenomeDoublings), by='DepMapID') %>%
  left_join(annotations %>% select(DepMapID, tcga_code)) %>%
  as.data.frame()
saveRDS(sample.amp.del, here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'sample.amp.del.rds'), compress=FALSE)
