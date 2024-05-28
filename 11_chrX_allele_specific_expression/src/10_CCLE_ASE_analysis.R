library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

gene.type.lo <- readRDS(file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'gene.type.lo.rds'))

exon.positions.df <- readRDS(file=here('11_chrX_allele_specific_expression/output/07_CCLE_annotate_SNPs', 'exon.positions.df.rds')) %>%
  mutate(contig=paste0('chr', CHROM)) %>%
  mutate(exon.intron='Exon')

sample.amp.del <- readRDS(here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'sample.amp.del.rds'))

Mahmoud.supp.file <- here('05_CCLE_data_preparation/data/CCLE_Mahmoud2019Nature', '41586_2019_1186_MOESM4_ESM.xlsx')
annotations <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Cell Line Annotations') %>%
  mutate(ModelID=sub('-', '.', depMapID))
datasets <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Datasets')

sample.with.wes <- datasets %>%
  filter(WES_CCLE=='TRUE')

sample.with.wes.depmapid <- annotations %>%
  filter(CCLE_ID %in% sample.with.wes$CCLE_ID) %>%
  pull(ModelID)


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
  group_by(ModelID, position) %>%
  mutate(gene.type.unique=gene.type %>% unique() %>% sort() %>% paste(collapse='_')) %>%
  mutate(gene.type.number=gene.type %>% unique() %>% length()) %>%
  ungroup() %>%
  filter(gene.type.number==1) %>%
  distinct(ModelID, position, .keep_all=TRUE) %>%
  mutate(arm=case_when(position < 58632012 ~ 'p', position > 61632012 ~ 'q')) %>% #Check centromere position with rCGH::hg38
  mutate(region=case_when(position < 2781479 ~ 'PAR1', position > 156030895 ~ 'PAR2', TRUE ~ 'non-PAR')) %>%
  left_join(sample.amp.del %>% mutate(contig=paste0('chr', chr)), by=c('ModelID', 'contig')) %>%
  mutate(Released.to.SRA=case_when(ModelID %in% sample.with.wes.depmapid ~ TRUE, TRUE ~ FALSE)) %>%
  as.data.frame()
saveRDS(ase.chrx.annot, file=here('11_chrX_allele_specific_expression/output/10_CCLE_ASE_analysis', 'ase.chrx.annot.rds'), compress=FALSE)

ase.chrx.annot.major.rate <- ase.chrx.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  mutate(major.rate=case_when(ref > alt ~ ref, ref < alt ~ alt, ref==alt ~ ref))

min.points.per.sample <- 20

ase.chrx.annot.major.rate.female <- ase.chrx.annot.major.rate %>%
  filter(Sex=='Female') %>%
  filter(exon.intron=='Exon') %>%
  filter(gene.type=='Inactive') %>%
  group_by(ModelID) %>%
  mutate(points.no=n()) %>%
  filter(points.no >= min.points.per.sample) %>%
  mutate(major.rate.median=median(major.rate)) %>%
  mutate(major.rate.mean=mean(major.rate)) %>%
  mutate(sd=sd(major.rate)) %>%
  ungroup() %>%
  filter(karyo.class %in% c('No.Arm-level.Alt', 'Whole.Amp')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('No.Arm-level.Alt', 'Whole.Amp'))) %>%
  filter(!is.na(ploidy.class)) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid'))) %>%
  mutate(tcga_code=factor(.$tcga_code, levels=.$tcga_code %>% unique() %>% sort())) %>%
  mutate(OncotreeLineage=factor(.$OncotreeLineage, levels=.$OncotreeLineage %>% unique() %>% sort())) %>%
  group_by(ModelID, karyo.class, ploidy.class, tcga_code, OncotreeLineage) %>%
  summarize(major.rate.median=unique(major.rate.median), major.rate.mean=unique(major.rate.mean), points.per.sample=n(), sd=unique(sd)) %>%
  ungroup()

## No CNAs/Whole amp vs ASE value (tumor types merged)
stat.test.uni <- ase.chrx.annot.major.rate.female %>%
  rstatix::wilcox_test(major.rate.median ~ karyo.class)

stat.test.sep <- ase.chrx.annot.major.rate.female %>%
  group_by(ploidy.class) %>%
  rstatix::wilcox_test(major.rate.median ~ karyo.class)

stat.test.ploidy <- ase.chrx.annot.major.rate.female %>%
  group_by(karyo.class) %>%
  rstatix::wilcox_test(major.rate.median ~ ploidy.class)

g <- ggplot(ase.chrx.annot.major.rate.female, aes(x=karyo.class, y=major.rate.median)) +
  geom_boxplot(aes(col=ploidy.class, group=interaction(ploidy.class, karyo.class, drop=FALSE)), position='dodge', outlier.shape=NA, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class), size=3, shape=21, position=position_jitterdodge(0.5), alpha=0.5) +
  scale_y_continuous(breaks=seq(0.5, 1, by=0.1)) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  coord_cartesian(ylim=c(0.5, 1), clip='off') +
  ggpubr::stat_pvalue_manual(stat.test.uni %>% rstatix::add_xy_position(x='karyo.class'), y.position=1.18, tip.length=0.1, label='P = {p}', col='black', size=5, show.legend=FALSE) +
  ggpubr::stat_pvalue_manual(stat.test.sep %>% rstatix::add_xy_position(x='karyo.class', group='ploidy.class'), y.position=1.05, step.increase=0.3, tip.length=0.1, label='P = {p}', col='ploidy.class', size=5, show.legend=FALSE) +
  labs(y='ASE value', col='Ploidy') +
  guides(col=guide_legend(override.aes=list(alpha=1))) +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(plot.margin=margin(5,1,2,1, unit='cm'))
ggsave(g, file=here('11_chrX_allele_specific_expression/output/10_CCLE_ASE_analysis', 'Fig4f.png'), dpi=100, width=9, height=7)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/10_CCLE_ASE_analysis', 'Fig4f.pdf'), width=9, height=7, useDingbats=TRUE)


## Confirm that all CCLE samples are publicly available
sra.file <- here('11_chrX_allele_specific_expression/data', 'SraRunTable.txt')
sra <- read.delim(sra.file, sep=',')

sra.wxs <- sra %>%
  filter(Assay.Type=='WXS')

sra.wxs.modelid <- sif %>%
  filter(CCLEName %in% sra.wxs$Sample.Name) %>%
  pull(ModelID)

intersect(ase.chrx.annot.major.rate.female$ModelID, sra.wxs.modelid) %>% length() #48
setdiff(ase.chrx.annot.major.rate.female$ModelID, sra.wxs.modelid) %>% length() #0
