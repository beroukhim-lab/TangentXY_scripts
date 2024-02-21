library(tidyverse)
library(here)

## Whole exome capture kit (developed by Agilent and customized by Broad Institute) can be downloaded from URL below
## https://bitbucket.org/cghub/cghub-capture-kit-info/src/b33c1cb33593/BI/vendor/Agilent/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed?at=master
exome.capture.kit.file <- here('02_TCGA_data_preparation/data', 'whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed')
exome.capture.kit.bed <- read.delim(exome.capture.kit.file, header=F) %>%
  setNames(c('Chr', 'Start', 'Stop', 'Target', 'Strand')) %>%
  mutate(Chr=paste0('chr', Chr)) %>%
  mutate(Chr=factor(.$Chr, levels=.$Chr %>% unique() %>% gtools::mixedsort())) %>%
  mutate(length=Stop - Start + 1)

## Liftover hg19 to GRCh38 using 'rtracklayer' package
library(rtracklayer)

grObject <- GRanges(seqnames=exome.capture.kit.bed$Chr, ranges=IRanges(start=exome.capture.kit.bed$Start, end=exome.capture.kit.bed$Stop))

chainObject <- import.chain(here('05_CCLE_data_preparation/data', 'hg19ToHg38.over.chain')) # Downloaded from UCSC Genome Browser

results <- as.data.frame(liftOver(grObject, chainObject)) %>%
  dplyr::rename(chr='seqnames')

results2 <- results %>%
  arrange(chr, start) %>%
  filter(start!=end)

chromosomes <- results$chr %>% levels()
for (i in 1:length(chromosomes)) {
  data.i <- results2 %>%
    filter(chr==chromosomes[i]) %>%
    select(group, chr, start, end) %>%
    mutate(row.num=1:n())
    j <- 2
  while (j <= nrow(data.i)) {
    data.j_1 <- data.i %>% filter(row.num==j-1)
    data.j <- data.i %>% filter(row.num==j)

    if (j >= 3) {
      data.top <- data.i %>% filter(row.num <= j - 2) %>%
        mutate(group=as.character(group))
    } else {
      data.top <- data.frame()
    }    
    data.bottom <- data.i %>% filter(row.num > j) %>%
        mutate(group=as.character(group))

    start.j <- data.j$start
    end.j_1 <- data.j_1$end

    if (start.j <= end.j_1) {
      new.data.j <- data.j_1 %>%
        bind_rows(data.j) %>%
        group_by(chr) %>%
        summarize(group=group %>% unique() %>% paste(collapse=':'),
          chr=chr %>% unique(),
          start=min(start),
          end=max(end),
          merge=row.num %>% unique() %>% paste(collapse=':'),
          row.num=max(row.num))
    } else {
       new.data.j <- data.j_1 %>%
        bind_rows(data.j) %>%
        mutate(group=as.character(group))
    }
    data.i <- data.top %>%
      bind_rows(new.data.j) %>%
      bind_rows(data.bottom)

    j <- j + 1
    print(paste0('chr', i, ':', j))
  }

  if (i==1) {
    results3 <- data.i
  } else {
    results3 <- results3 %>% bind_rows(data.i)
  }
}

results4 <- results3 %>%
  dplyr::rename(Chr='chr', Start='start', Stop='end')

## Add annotation about common CNV
common.germ.cnv.file <- here('05_CCLE_data_preparation/data', 'hg38_GDC_SNP6_CNV_list.161107.txt') # Got personally from Andrew Cherniack. Attached to the email from Andrew on Mar 8, 2022
common.germ.cnv <- read_delim(common.germ.cnv.file) %>%
  mutate(chr=paste0('chr', chr)) %>%
  mutate(length=end-start)
probes.with.flag <- results4 %>%
  left_join(common.germ.cnv[, c('chr', 'start', 'end')], by=c('Chr'='chr')) %>%
  mutate(overlap=case_when(Stop < start | Start > end ~ FALSE, TRUE ~ TRUE)) %>%
  mutate(overlap=case_when(Chr=='chrY' ~ FALSE, TRUE ~ overlap))
probes.to.be.excluded <- probes.with.flag %>%
  filter(overlap==TRUE) %>%
  select(Chr, Start, Stop, overlap) %>%
  distinct()
results5 <- results4 %>%
  left_join(probes.to.be.excluded, by=c('Chr', 'Start', 'Stop')) %>%
  dplyr::rename(freqcnv='overlap') %>%
  mutate(freqcnv=case_when(is.na(freqcnv) ~ FALSE, TRUE ~ freqcnv))

## Add annotation about centromere
hg38 <- rCGH::hg38 %>%
  select(chrom, length, centromerStart, centromerEnd) %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom))) %>%
  mutate(chrom=paste0('chr', chrom)) %>%
  mutate(centromerLength=centromerEnd - centromerStart) %>%
  mutate(centromerCenter=centromerStart + centromerLength/2)

results6 <- results5 %>%
  left_join(hg38, by=c('Chr'='chrom')) %>%
  mutate(Arm=case_when(Stop < centromerStart ~ 'p',
                        Start > centromerEnd ~ 'q',
                        Start > centromerStart & Stop < centromerEnd ~ 'c')) %>%
  mutate(Arm2=case_when(Stop < centromerCenter ~ 'p',
                        Start > centromerCenter ~ 'q'))

## Add annotation about PAR
par.file <- here('02_TCGA_data_preparation/data', 'PAR.txt')
par <- read_delim(par.file)
par.hg38 <- par %>%
  filter(Version=='GRCh38')

results7 <- results6 %>%
  mutate(Par=case_when(Chr=='chrX' & Stop < par.hg38 %>% filter(Chr=='X' & Name=='par1') %>% pull(Stop) ~ 'par1',
                        Chr=='chrX' & Start > par.hg38 %>% filter(Chr=='X' & Name=='par2') %>% pull(Start) ~ 'par2',
                        Chr=='chrY' & Stop < par.hg38 %>% filter(Chr=='Y' & Name=='par1') %>% pull(Stop) ~ 'par1',
                        Chr=='chrY' & Start > par.hg38 %>% filter(Chr=='Y' & Name=='par2') %>% pull(Start) ~ 'par2')) %>%
  mutate(Masked=case_when(Chr=='chrY' & !is.na(Par) ~ TRUE, TRUE ~ FALSE))

exome.capture.kit.hg38.data <- results7 %>%
  select(Chr, Start, Stop, Par, Masked, freqcnv, length, centromerStart, centromerEnd, centromerLength, centromerCenter, Arm, Arm2)
write.table(exome.capture.kit.hg38.data, file=here('05_CCLE_data_preparation/output/01_lifgOver_Hg19_to_GRCh38', 'exomeAgilentCaptureKit_GRCh38_data.txt'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
