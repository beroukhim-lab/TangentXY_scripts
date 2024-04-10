library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

gene.type.file <- here('10_chrX_SCNA_vs_gene_expression/data/XCI_Tukiainen2017Nature', 'Suppl.Table.1.xlsx')
gene.type <- readxl::read_xlsx(gene.type.file, skip=1) %>%
  as.data.frame() %>%
  setNames(gsub(' ', '_', colnames(.))) %>%
  mutate(ensembl_gene_id=sub('\\..*$', '', Gene_ID)) %>%
  mutate(Gene_name=case_when(Gene_name=='42619' ~ 'SEPTIN6', TRUE ~ Gene_name))

gene.type.list <- split(gene.type, seq(nrow(gene.type)))
expand.position <- function(df) {
  start <- df$Start_position
  end <- df$End_position
  diff <- end - start + 1
  df.rep <- df %>%
    slice(rep(1, each=diff)) %>%
    mutate(position=start:end)
  return(df.rep)
}
gene.type.expanded.list <- lapply(gene.type.list, expand.position)
gene.type.expanded <- gene.type.expanded.list %>%
  bind_rows()

snp.annot.df <- readRDS(file=here('11_chrX_allele_specific_expression/output/03_TCGA_annotate_SNPs', 'snp.annot.df.rds'))

sample.amp.del <- readRDS(file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'sample.amp.del.rds')) %>%
  mutate(TCGA.SampleID=paste(TCGA.ID, type, sep='.'))

terra.table <- read.delim(file=here('11_chrX_allele_specific_expression/output/02_TCGA_make_input_table_for_Terra', 'sample_table_for_Terra.tsv'))


## Analysis on chrX
ase.chrx.dirs <- Sys.glob(file.path(here('11_chrX_allele_specific_expression/data/TCGA_Terra_output'), 'chrX', '*', 'ASEReadCounterWorkflow', '*', 'call-ASEReadCounter'))
ase.chrx.files <- list.files(ase.chrx.dirs , pattern='_ASE_chrX.txt$', recursive=TRUE, full.names=TRUE)
ase.chrx.ids <- sub('_ASE_chrX.txt$', '', basename(ase.chrx.files)) %>%
  gsub('-', '.', .)
ase.chrx.list <- lapply(ase.chrx.files, read_delim) %>%
  setNames(ase.chrx.ids)

ase.chrx.df <- ase.chrx.list %>%
  Map(bind_cols, ., SampleID=names(.)) %>%
  bind_rows() %>%
  as.data.frame()
saveRDS(ase.chrx.df, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'ase.chrx.df.rds'), compress=FALSE)

ase.chrx.annot <- ase.chrx.df %>%
  mutate(Chr=sub('chr', '', contig)) %>%
  left_join(gene.type.expanded, by=c('Chr', 'position')) %>%
  rename(gene.type=Combined_XCI_status) %>%
  filter(!is.na(gene.type)) %>%
  mutate(gene.type=factor(.$gene.type, levels=c('inactive', 'escape', 'variable'))) %>%
  group_by(SampleID, position) %>%
  mutate(gene.type.unique=gene.type %>% unique() %>% sort() %>% paste(collapse='_')) %>%
  mutate(gene.type.number=gene.type %>% unique() %>% length()) %>%
  ungroup() %>%
  filter(gene.type.number==1) %>%
  distinct(SampleID, position, .keep_all=TRUE) %>%
  mutate(arm=case_when(position < 58632012 ~ 'p', position > 61632012 ~ 'q')) %>% #Check centromere position with rCGH::hg19
  mutate(region=case_when(position < 2699520 ~ 'PAR1', position > 155260560 ~ 'PAR2', TRUE ~ 'non-PAR')) %>%
  left_join(snp.annot.df %>% rename(exon.intron=region), by=c('Chr'='CHROM', 'position'='POS')) %>%
  left_join(sample.amp.del, by=c('SampleID', 'Chr'='chr')) %>%
  as.data.frame()
saveRDS(ase.chrx.annot, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'ase.chrx.annot.rds'), compress=FALSE)

ase.chrx.annot.major.rate <- ase.chrx.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  mutate(major.rate=case_when(ref > alt ~ ref, ref < alt ~ alt, ref==alt ~ ref))

min.points.per.sample <- 20

ase.chrx.annot.major.rate.female <- ase.chrx.annot.major.rate %>%
  filter(Gender=='Female') %>%
  filter(exon.intron=='exon') %>%
  filter(gene.type=='inactive') %>%
  group_by(SampleID) %>%
  mutate(points.no=n()) %>%
  filter(points.no >= min.points.per.sample) %>%
  mutate(major.rate.median=median(major.rate)) %>%
  mutate(major.rate.mean=mean(major.rate)) %>%
  mutate(sd=sd(major.rate)) %>%
  ungroup() %>%
  filter(karyo.class %in% c('No.Alt', 'Whole.Amp')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('No.Alt', 'Whole.Amp'))) %>%
  mutate(ploidy.class=case_when(ploidy < 2.5 & Genome.doublings==0 ~ 'Diploid', ploidy >= 2.5 | Genome.doublings > 0 ~ 'Polyploid')) %>%
  filter(!is.na(ploidy.class)) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid'))) %>%
  mutate(project=factor(.$project, levels=.$project %>% unique() %>% sort())) %>%
  group_by(SampleID, karyo.class, ploidy.class, project) %>%
  summarize(major.rate.median=unique(major.rate.median), major.rate.mean=unique(major.rate.mean), points.per.sample=n(), sd=unique(sd)) %>%
  ungroup()
saveRDS(ase.chrx.annot.major.rate.female, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'ase.chrx.annot.major.rate.female.rds'), compress=FALSE)

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
  geom_point(aes(col=ploidy.class), size=3, shape=21, alpha=0.5, position=position_jitterdodge(0.5)) +
  scale_x_discrete(labels=c('No CNAs', 'Amp')) +
  scale_y_continuous(breaks=seq(0.5, 1, by=0.1)) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  coord_cartesian(ylim=c(0.5, 1), clip='off') +
  ggpubr::stat_pvalue_manual(stat.test.uni %>% rstatix::add_xy_position(x='karyo.class'), y.position=1.18, label='P = {p}', col='black', size=5, show.legend=FALSE) +
  ggpubr::stat_pvalue_manual(stat.test.sep %>% rstatix::add_xy_position(x='karyo.class', group='ploidy.class'), y.position=1.05, step.increase=0.1, label='P = {p}', col='ploidy.class', size=5, show.legend=FALSE) +
  labs(y='ASE value', col='Ploidy') +
  guides(col=guide_legend(override.aes=list(alpha=1))) +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank()) +
  theme(plot.margin=margin(5,1,2,1, unit='cm'))
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'Fig4e.png'), dpi=100, width=9, height=7)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'Fig4e.pdf'), width=9, height=7, useDingbats=TRUE)


## No CNAs/Whole amp vs ASE value (tumor types separated)
dummy.data <- ase.chrx.annot.major.rate.female %>%
  group_by(project, karyo.class, ploidy.class, .drop=FALSE) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  filter(n==0) %>%
  mutate(major.rate.median=0, major.rate.mean=0) %>%
  slice(rep(1:n(), each=10))

object.for.plot <- ase.chrx.annot.major.rate.female %>%
  bind_rows(dummy.data)

g <- ggplot(object.for.plot, aes(x=karyo.class, y=major.rate.median)) +
  geom_boxplot(aes(col=ploidy.class, group=interaction(ploidy.class, karyo.class, drop=FALSE)), position='dodge', outlier.shape=NA, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class), size=3, shape=21, position=position_jitterdodge(), alpha=0.75) +
  scale_x_discrete(labels=c('No CNAs', 'Amp')) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  coord_cartesian(ylim=c(0.5, 1)) +
  lemon::facet_rep_wrap(~project, nrow=5, drop=F, repeat.tick.labels=TRUE) +
  labs(y='ASE value', col='Ploidy', title='ChrX') +
  guides(col=guide_legend(override.aes=list(alpha=1))) +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'FigS6a.png'), dpi=100, width=20, height=14)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'FigS6a.pdf'), width=20, height=14, useDingbats=TRUE)




## Analysis on autosomes
tumor.types <- c('COAD', 'GBM', 'LUAD')
for (i in seq_along(tumor.types)) {
  tumor_i <- tumor.types[i]
  print(tumor_i)

  ase.auto.dirs <- Sys.glob(file.path(here('11_chrX_allele_specific_expression/data/TCGA_Terra_output'), 'autosomes', tumor_i, '*', 'ASEReadCounterWorkflow', '*', 'call-ASEReadCounter'))
  ase.auto.files <- list.files(ase.auto.dirs, , pattern='_ASE_autosomes.txt$', recursive=TRUE, full.names=TRUE)
  ase.auto.ids <- sub('_ASE_autosomes.txt$', '', basename(ase.auto.files)) %>%
    gsub('-', '.', .)
  ase.auto.list <- lapply(ase.auto.files, read_delim) %>%
    setNames(ase.auto.ids)

  ase.auto.df <- ase.auto.list %>%
    Map(bind_cols, ., SampleID=names(.)) %>%
    bind_rows()
  saveRDS(ase.auto.df, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', paste0('ase.auto.df_', tumor_i, '.rds')), compress=FALSE)

  ase.auto.annot <- ase.auto.df %>%
    mutate(Chr=sub('chr', '', contig)) %>%
    left_join(snp.annot.df %>% rename(exon.intron=region), by=c('Chr'='CHROM', 'position'='POS')) %>%
    left_join(sample.amp.del, by=c('SampleID', 'Chr'='chr')) %>%
    as.data.frame()
  saveRDS(ase.auto.annot, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', paste0('ase.auto.annot_', tumor_i, '.rds')), compress=FALSE)
}


tumor.type <- 'GBM'
ase.auto.annot <- readRDS(file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', paste0('ase.auto.annot_', tumor.type, '.rds')))

ase.auto.annot.major.rate <- ase.auto.annot %>%
  mutate(lowMapQ.rate=lowMAPQDepth/rawDepth) %>%
  mutate(lowBaseQ.rate=lowBaseQDepth/rawDepth) %>%
  mutate(ref=refCount/totalCount) %>%
  mutate(alt=altCount/totalCount) %>%
  mutate(major.rate=case_when(ref > alt ~ ref, ref < alt ~ alt, ref==alt ~ ref))

ase.auto.annot.major.rate.summary <- ase.auto.annot.major.rate %>%
  filter(exon.intron=='exon') %>%
  group_by(SampleID, contig) %>%
  mutate(points.no=n()) %>%
  filter(points.no >= min.points.per.sample) %>%
  mutate(major.rate.median=median(major.rate)) %>%
  mutate(major.rate.mean=mean(major.rate)) %>%
  mutate(sd=sd(major.rate)) %>%
  ungroup() %>%
  mutate(contig=factor(.$contig, levels=.$contig %>% unique() %>% gtools::mixedsort())) %>%
  filter(karyo.class %in% c('No.Alt', 'Whole.Amp')) %>%
  mutate(karyo.class=factor(.$karyo.class, levels=c('No.Alt', 'Whole.Amp'))) %>%
  mutate(ploidy.class=case_when(ploidy < 2.5 & Genome.doublings==0 ~ 'Diploid', ploidy >= 2.5 | Genome.doublings > 0 ~ 'Polyploid')) %>%
  filter(!is.na(ploidy.class)) %>%
  mutate(ploidy.class=factor(.$ploidy.class, levels=c('Diploid', 'Polyploid'))) %>%
  group_by(SampleID, Gender, contig, karyo.class, ploidy.class, project, sd) %>%
  summarize(major.rate.median=unique(major.rate.median), major.rate.mean=unique(major.rate.mean)) %>%
  ungroup()

dummy.data.auto <- ase.auto.annot.major.rate.summary %>%
  group_by(contig, karyo.class, ploidy.class, .drop=FALSE) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  filter(n==0) %>%
  mutate(major.rate.median=0) %>%
  slice(rep(1:n(), each=10))

object.for.plot.auto <- ase.auto.annot.major.rate.summary %>%
  bind_rows(dummy.data.auto)

g <- ggplot(object.for.plot.auto, aes(x=karyo.class, y=major.rate.median)) +
  geom_boxplot(aes(col=ploidy.class, group=interaction(ploidy.class, karyo.class, drop=FALSE)), position='dodge', outlier.shape=NA, show.legend=FALSE) +
  geom_point(aes(col=ploidy.class, shape=Gender), size=3, position=position_jitterdodge(), alpha=0.5) +
  scale_x_discrete(labels=c('No CNAs', 'Amp')) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  scale_shape_manual(values=c('Female'=21, 'Male'=24)) +
  coord_cartesian(ylim=c(0.5, 1)) +
  lemon::facet_rep_wrap(~contig, nrow=4, drop=F, repeat.tick.labels=TRUE) +
  labs(y='ASE value', col='Ploidy', title=tumor.type) +
  guides(col=guide_legend(override.aes=list(alpha=1, shape=1, size=5, order=1)), shape=guide_legend(override.aes=list(alpha=1, size=5, order=2))) +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'FigS6b.png'), dpi=100, width=20, height=14)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'FigS6b.pdf'), width=20, height=14, useDingbats=TRUE)
