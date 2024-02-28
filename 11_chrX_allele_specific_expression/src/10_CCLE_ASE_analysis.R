library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

gene.type.lo <- readRDS(file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'gene.type.lo.rds'))

exon.positions.df <- readRDS(file=here('11_chrX_allele_specific_expression/output/07_CCLE_annotate_SNPs', 'exon.positions.df.rds')) %>%
  mutate(contig=paste0('chr', CHROM)) %>%
  mutate(exon.intron='Exon')

sample.amp.del <- readRDS(here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'sample.amp.del.rds'))

## Ploidy
polyploidy.threshold <- 2.5

Mahmoud.supp.file <- here('05_CCLE_data_preparation/data/CCLE_Mahmoud2019Nature', '41586_2019_1186_MOESM4_ESM.xlsx')
annotations <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Cell Line Annotations') %>%
  mutate(ModelID=sub('-', '.', depMapID))
datasets <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Datasets')

sample.with.wes <- datasets %>%
  filter(WES_CCLE=='TRUE')

sample.with.wes.depmapid <- annotations %>%
  filter(CCLE_ID %in% sample.with.wes$CCLE_ID) %>%
  pull(ModelID)

# cbio.file <- here('data', 'ccle_broad_2019_clinical_data.tsv') # Downloaded from cBioPortal (https://www.cbioportal.org/study/clinicalData?id=ccle_broad_2019)
# cbio <- read.delim(cbio.file) %>%
#   setNames(gsub('\\.', '', colnames(.))) %>%
#   mutate(DepMapID=gsub('-', '.', DepMapID)) %>%
#   filter(!is.na(DepMapID)) %>%
#   filter(AnnotationSource=='CCLE') %>%
#   left_join(annotations %>% select(DepMapID, tcga_code), by='DepMapID') %>%
#   mutate(ploidy.class=case_when(is.na(Ploidy) & is.na(GenomeDoublings) ~ NA,
#                                 Ploidy < polyploidy.threshold & GenomeDoublings==0 ~ 'Diploid',
#                                 TRUE ~ 'Multiploid'))

# ploidy <- cbio %>%
#   select(DepMapID, Purity, Ploidy, GenomeDoublings, ploidy.class)


## Load Terra output
ase.dirs <- Sys.glob(file.path(here('11_chrX_allele_specific_expression/data/CCLE_Terra_output/ASEReadCounterWorkflow'), '*', 'call-ASEReadCounter'))
ase.files <- list.files(ase.dirs , pattern='_ASE.txt$', recursive=TRUE, full.names=TRUE)
ase.ids <- sub('_ASE.txt$', '', basename(ase.files)) %>%
  gsub('-', '.', .)
ase.list <- lapply(ase.files, read_delim) %>%
  setNames(ase.ids)

ase.df <- ase.list %>%
  Map(bind_cols, ., ModelID=names(.)) %>%
  bind_rows() %>%
  as.data.frame()
saveRDS(ase.df, file=here('11_chrX_allele_specific_expression/output/10_CCLE_ASE_analysis', 'ase.df.rds'), compress=FALSE)


## ChrX
ase.chrx.annot <- ase.df %>%
  filter(contig=='chrX') %>%
  left_join(gene.type.lo %>% mutate(contig=paste0('chr', Chr)), by='contig') %>%
  mutate(hit=case_when(start.hg38 <= position & position <= end.hg38 ~ 'hit')) %>%
  filter(!is.na(hit)) %>%
  left_join(exon.positions.df %>% filter(CHROM=='X'), by=c('contig', 'position'='POS')) %>%
  mutate(exon.intron=case_when(is.na(exon.intron) ~ 'Intron', TRUE ~ exon.intron)) %>%
  rename(gene.type=Combined_XCI_status) %>%
  filter(!is.na(gene.type)) %>%
  mutate(gene.type=str_to_title(gene.type)) %>%
  mutate(gene.type=factor(.$gene.type, levels=c('Inactive', 'Escape', 'Variable'))) %>%
  group_by(DepMapID, position) %>%
  mutate(gene.type.unique=gene.type %>% unique() %>% sort() %>% paste(collapse='_')) %>%
  mutate(gene.type.number=gene.type %>% unique() %>% length()) %>%
  ungroup() %>%
  filter(gene.type.number==1) %>%
  distinct(DepMapID, position, .keep_all=TRUE) %>%
  mutate(arm=case_when(position < 58632012 ~ 'p', position > 61632012 ~ 'q')) %>% #Check centromere position with rCGH::hg38
  mutate(region=case_when(position < 2781479 ~ 'PAR1', position > 156030895 ~ 'PAR2', TRUE ~ 'non-PAR')) %>%
  left_join(sample.amp.del %>% mutate(contig=paste0('chr', chr)), by=c('DepMapID', 'contig')) %>%
  mutate(Released.to.SRA=case_when(DepMapID %in% sample.with.wes.depmapid ~ TRUE, TRUE ~ FALSE)) %>%
  as.data.frame()
saveRDS(ase.chrx.annot, file=here('output/21_ASEReadCounter', 'ase.chrx.annot.RData'), compress=FALSE)
# ase.chrx.annot <- readRDS(file=here('output/21_ASEReadCounter', 'ase.chrx.annot.RData'))

ase.chrx.annot.major.rate <- ase.chrx.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  mutate(major.rate=case_when(ref > alt ~ ref, ref < alt ~ alt, ref==alt ~ ref))

ase.chrx.annot.major.rate.summary <- ase.chrx.annot.major.rate %>%
  filter(exon.intron=='Exon') %>%
  mutate(ploidy.class=case_when(Ploidy < 2.5 & GenomeDoublings==0 ~ 'Diploid', Ploidy >= 2.5 | GenomeDoublings > 0 ~ 'Polyploid')) %>%
  group_by(DepMapID, gene.type) %>%
  mutate(major.rate.median=median(major.rate)) %>%
  mutate(sd=sd(major.rate)) %>%
  mutate(points.no=n()) %>%
  ungroup()

lowBaseQ.threshold <- 1
min.points.per.sample <- 20

ase.chrx.annot.major.rate.summary2 <- ase.chrx.annot.major.rate.summary %>%
  filter(sex=='Female') %>%
  filter(!is.na(ploidy.class)) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid'))) %>%
  mutate(tcga_code=factor(.$tcga_code, levels=.$tcga_code %>% unique() %>% sort())) %>%
  mutate(lineage=factor(.$lineage, levels=.$lineage %>% unique() %>% sort())) %>%
  filter(lowBaseQ.rate < lowBaseQ.threshold) %>%
  filter(gene.type=='Inactive') %>%
  filter(karyo.class %in% c('No.Alt', 'Whole.Amp')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('No.Alt', 'Whole.Amp'))) %>%
  group_by(DepMapID, sex, karyo.class, ploidy.class, tcga_code, Released.to.SRA) %>%
  summarize(major.rate.median=median(major.rate), points.per.sample=n(), sd=unique(sd)) %>%
  ungroup() %>%
  filter(points.per.sample >= min.points.per.sample) %>%
  mutate(XCI=case_when(karyo.class=='No.Alt' & major.rate.median < 0.9 ~ 'Dysregulated', TRUE ~ 'Normal')) %>%
  mutate(XCI=factor(.$XCI, levels=c('Normal', 'Dysregulated')))

g <- ggplot(ase.chrx.annot.major.rate.summary2, aes(x=karyo.class, y=major.rate.median)) +
  geom_violin(aes(col=ploidy.class), position='dodge', show.legend=FALSE) +
  # geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge2(0.75, preserve='single')) +
  geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge(0.9), width=0.2, fatten=3, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class, fill=sd, size=points.per.sample), shape=21, position=position_jitterdodge(), alpha=0.75) +
  scale_fill_gradient2(low='limegreen', mid='white', high='red', midpoint=0.1) +
  ylim(0.5, 1) +
  facet_wrap(~tcga_code, nrow=5, drop=F) +
  labs(title=paste0('Female, Inactive genes, LowBaseQualityCutoff=', lowBaseQ.threshold, ', At least ', min.points.per.sample, ' points/sample')) +
  guides(size=guide_legend(order=1), color=guide_legend(order=2), fill=guide_colorbar(order=3)) +
  theme_bw(base_size=20) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/21_ASEReadCounter', 'ASE_chrX_TumorTypesSeparated.png'), dpi=100, width=28, height=20)

g <- ggplot(ase.chrx.annot.major.rate.summary2 %>% filter(Released.to.SRA==TRUE), aes(x=karyo.class, y=major.rate.median)) +
  geom_violin(aes(col=ploidy.class), position='dodge', show.legend=FALSE) +
  # geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge2(0.75, preserve='single')) +
  geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge(0.9), width=0.2, fatten=3, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class, fill=tcga_code, size=points.per.sample), shape=21, position=position_jitterdodge(), alpha=0.75) +
  # scale_fill_gradient2(low='limegreen', mid='white', high='red', midpoint=0.1) +
  ylim(0.5, 1) +
  labs(title=paste0('Female, Inactive genes, LowBaseQualityCutoff=', lowBaseQ.threshold, ', At least ', min.points.per.sample, ' points/sample')) +
  guides(size=guide_legend(order=1), color=guide_legend(order=2), fill=guide_legend(order=3, ncol=2)) +
  theme_bw(base_size=20) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/21_ASEReadCounter', 'ASE_chrX_TumorTypesMerged.png'), dpi=100, width=12, height=7)


ase.chrx.annot.l <- ase.chrx.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  pivot_longer(names_to='ref.alt', values_to='rate', cols=c(ref, alt)) %>%
  as.data.frame()

sample.of.interest <- 'ACH.000173'
sample.of.interest <- 'ACH.000020'
sample.of.interest <- 'ACH.000129'
sample.of.interest <- 'ACH.000205'
sample.of.interest <- 'ACH.000281'
sample.of.interest <- 'ACH.000630'
sample.of.interest <- 'ACH.000640'
sample.of.interest <- 'ACH.000808'
sample.of.interest <- 'ACH.000876'

g <- ggplot(ase.chrx.annot.l %>% filter(DepMapID==sample.of.interest), aes(x=position, y=rate)) +
  geom_point(aes(size=totalCount, fill=lowBaseQ.rate, shape=exon.intron), col='white') +
  geom_hline(yintercept=0.5, col='darkgray', linetype='dashed') +
  ggrepel::geom_label_repel(data=. %>% filter(Gene_name=='XIST'), aes(label=Gene_name), size=5, box.padding=1, max.overlaps=Inf, show.legend=FALSE) +
  facet_wrap(~gene.type, ncol=1) +
  scale_shape_manual(values=list('Exon'=21, 'Intron'=24)) +
  scale_x_continuous(breaks=seq(0, 155270560, 50*1e6), labels=paste0(seq(0, 150, 50), 'Mb')) +
  labs(title=sample.of.interest, x='Position on chrX', y='Rate', size='#Total read') +
  guides(shape=guide_legend(override.aes=list(shape=1, col='black', size=3), order=1), size=guide_legend(override.aes=list(shape=1, col='black'), order=2), fill=guide_colorbar(order=3)) +
  theme_bw(base_size=20)
ggsave(g, file=here('output/21_ASEReadCounter', paste0('ASE_chrX_alongWithLoci_', sample.of.interest, '.png')), dpi=100, width=10, height=6)


## Autosomes
ase.auto.annot <- ase.df %>%
  filter(!contig %in% c('chrX', 'chrY')) %>%
  left_join(exon.positions.df, by=c('contig', 'position'='POS')) %>%
  left_join(sample.amp.del %>% mutate(contig=paste0('chr', chr)), by=c('DepMapID', 'contig')) %>%
  mutate(Released.to.SRA=case_when(DepMapID %in% sample.with.wes.depmapid ~ TRUE, TRUE ~ FALSE)) %>%
  as.data.frame()
saveRDS(ase.auto.annot, file=here('output/21_ASEReadCounter', 'ase.auto.annot.RData'), compress=FALSE)
# ase.auto.annot <- readRDS(file=here('output/21_ASEReadCounter', 'ase.auto.annot.RData'))

ase.auto.annot.major.rate <- ase.auto.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  mutate(major.rate=case_when(ref > alt ~ ref, ref < alt ~ alt, ref==alt ~ ref))

ase.auto.annot.major.rate.summary <- ase.auto.annot.major.rate %>%
  filter(exon.intron=='Exon') %>%
  mutate(ploidy.class=case_when(Ploidy < 2.5 & GenomeDoublings==0 ~ 'Diploid', Ploidy >= 2.5 | GenomeDoublings > 0 ~ 'Polyploid')) %>%
  group_by(DepMapID) %>%
  mutate(major.rate.median=median(major.rate)) %>%
  mutate(sd=sd(major.rate)) %>%
  mutate(points.no=n()) %>%
  ungroup()

ase.auto.annot.major.rate.summary2 <- ase.auto.annot.major.rate.summary %>%
  mutate(contig=factor(.$contig, levels=.$contig %>% unique() %>% gtools::mixedsort())) %>%
  filter(!is.na(ploidy.class)) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid'))) %>%
  filter(lowBaseQ.rate < lowBaseQ.threshold) %>%
  filter(karyo.class %in% c('No.Alt', 'Whole.Amp')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('No.Alt', 'Whole.Amp'))) %>%
  group_by(DepMapID, sex, contig, karyo.class, ploidy.class, tcga_code, Released.to.SRA) %>%
  summarize(major.rate.median=median(major.rate), points.per.sample=n()) %>%
  ungroup() %>%
  filter(points.per.sample >= min.points.per.sample)

tcga.code <- 'GBM'
g <- ggplot(ase.auto.annot.major.rate.summary2 %>% filter(Released.to.SRA==TRUE & tcga_code==tcga.code), aes(x=karyo.class, y=major.rate.median)) +
  geom_violin(aes(col=ploidy.class), position='dodge', show.legend=FALSE) +
  # geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge2(0.75, preserve='single')) +
  geom_boxplot(aes(col=ploidy.class), outlier.shape=NA, position=position_dodge(0.9), width=0.2, fatten=3, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class, shape=sex), position=position_jitterdodge(), alpha=0.75) +
  # scale_fill_gradient2(low='limegreen', mid='white', high='red', midpoint=0.1) +
  facet_wrap(~ contig, nrow=4, drop=FALSE) +
  ylim(0.5, 1) +
  labs(title=tcga.code) +
  guides(size=guide_legend(order=1), color=guide_legend(order=2), fill=guide_legend(order=3, ncol=2)) +
  theme_bw(base_size=20) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output/21_ASEReadCounter', paste0('ASE_auto_', tcga.code, '.png')), dpi=100, width=12, height=7)
