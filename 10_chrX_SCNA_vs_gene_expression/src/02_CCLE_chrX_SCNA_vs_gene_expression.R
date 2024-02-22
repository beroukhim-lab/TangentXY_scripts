library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))

## Ploidy
polyploidy.threshold <- 2.5

Mahmoud.supp.file <- here('05_CCLE_data_preparation/data/CCLE_Mahmoud2019Nature', '41586_2019_1186_MOESM4_ESM.xlsx')
annotations <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Cell Line Annotations') %>%
  mutate(DepMapID=sub('-', '.', depMapID))

datasets <- readxl::read_xlsx(Mahmoud.supp.file, sheet='Datasets')

sample.with.wes <- datasets %>%
  filter(WES_CCLE=='TRUE')

sample.with.wes.depmapid <- annotations %>%
  filter(CCLE_ID %in% sample.with.wes$CCLE_ID) %>%
  pull(depMapID) %>%
  sub('-', '.', .)

samples.to.be.analyzed <- sif %>%
  filter(DepMap_ID %in% sample.with.wes.depmapid) %>%
  filter(lineage!='fibroblast') %>%
  pull(DepMap_ID)

cbio.file <- here('07_SCNAs_in_chrX_and_chrY/data', 'ccle_broad_2019_clinical_data.tsv') # Downloaded from cBioPortal (https://www.cbioportal.org/study/clinicalData?id=ccle_broad_2019)
cbio <- read.delim(cbio.file) %>%
  setNames(gsub('\\.', '', colnames(.))) %>%
  mutate(DepMapID=gsub('-', '.', DepMapID)) %>%
  filter(!is.na(DepMapID)) %>%
  filter(AnnotationSource=='CCLE') %>%
  left_join(annotations %>% select(DepMapID, tcga_code), by='DepMapID') %>%
  mutate(ploidy.class=case_when(is.na(Ploidy) & is.na(GenomeDoublings) ~ NA,
                                Ploidy < polyploidy.threshold & GenomeDoublings==0 ~ 'Diploid',
                                TRUE ~ 'Polyploid'))

## Analysis on chrX
## Gene list (Normal gene/Escape gene)
gene.type.file <- here('10_chrX_SCNA_vs_gene_expression/data/XCI_Tukianinen2017Nature', 'Suppl.Table.1.xlsx')
gene.type <- readxl::read_xlsx(gene.type.file, skip=1) %>%
  as.data.frame() %>%
  setNames(gsub(' ', '_', colnames(.))) %>%
  mutate(ensembl_gene_id=sub('\\..*$', '', Gene_ID)) %>%
  mutate(Gene_name=case_when(Gene_name=='42619' ~ 'SEPTIN6', TRUE ~ Gene_name))

## Lift over gene.type object from hg19 to hg38
# library(rtracklayer)
chainObject <- rtracklayer::import.chain('05_CCLE_data_preparation/data/hg19ToHg38.over.chain')

for (i in 1:nrow(gene.type)) {
  gene.type.i <- gene.type[i, ]
  gene.name.i <- gene.type.i$Gene_name
  print(paste(i, gene.name.i))

  grObject <- GenomicRanges::GRanges(seqnames='chrX', ranges=IRanges::IRanges(start=gene.type.i$Start_position, end=gene.type.i$End_position))
  lo <- as.data.frame(rtracklayer::liftOver(grObject, chainObject)) %>%
    filter(seqnames=='chrX') %>%
    arrange(start, end) %>%
    mutate(Chr=sub('chr', '', seqnames)) %>%
    rename(start.hg38=start, end.hg38=end)

  gene.type.i.lo <- gene.type.i %>%
    left_join(lo %>% select(Chr, start.hg38, end.hg38), by='Chr')

  if(nrow(gene.type.i.lo) >= 2) {
    # gene.type.i.lo <- gene.type.i.lo %>%
    #   arrange(start.hg38, end.hg38) %>%
    #   mutate(pre.end=lag(end.hg38)) %>%
    #   mutate(diff=start.hg38-pre.end)

    # print(paste(i, nrow(gene.type.i.lo), max(gene.type.i.lo$diff, na.rm=T)))

    if (gene.name.i=='HTR2C') {
      gene.type.i.lo <- gene.type.i.lo %>%
        filter(end.hg38!=78550462)
    }

    gene.type.i.lo <- gene.type.i.lo %>%
      group_by(Gene_ID) %>%
      mutate(start.min=min(start.hg38), end.max=max(end.hg38)) %>%
      select(-c('start.hg38', 'end.hg38')) %>%
      rename(start.hg38=start.min, end.hg38=end.max) %>%
      distinct()
  }

  gene.type.i.lo <- gene.type.i.lo %>%
    mutate(gene.length.hg19=End_position - Start_position) %>%
    mutate(gene.length.hg38=end.hg38 - start.hg38) %>%
    mutate(length.ratio=gene.length.hg38/gene.length.hg19)

  if (i==1) {
    gene.type.lo <- gene.type.i.lo
  } else {
    gene.type.lo <- gene.type.lo %>%
      bind_rows(gene.type.i.lo)
  }
}
saveRDS(gene.type.lo, file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'gene.type.lo.rds'), compress=FALSE)

## CN
segment.arm <- readRDS(fil=here('07_SCNAs_in_chrX_and_chrY/output/02_CCLE_SCNA_classification', 'segment.arm.rds'))

OmicsDefaultModelProfiles.file <- here('10_chrX_SCNA_vs_gene_expression/data' , 'OmicsOmicsDefaultModelProfiles.csv')
OmicsExpressionAllGenesTPMLogp1Profile.file <- here('10_chrX_SCNA_vs_gene_expression/data', 'OmicsExpressionAllGenesTPMLogp1Profile.csv')

# OmicsProfiles <- read_csv(OmicsProfiles.file)
OmicsDefaultModelProfiles <- read_csv(OmicsDefaultModelProfiles.file)
OmicsExpressionAllGenesTPMLogp1Profile <- read_csv(OmicsExpressionAllGenesTPMLogp1Profile.file)

log2tpm.df <- OmicsExpressionAllGenesTPMLogp1Profile %>%
  rename(ProfileID='...1') %>%
  left_join(OmicsDefaultModelProfiles %>% select(ProfileID, ModelID), by='ProfileID') %>%
  filter(!is.na(ModelID)) %>%
  mutate(DepMapID=sub('-', '.', ModelID)) %>%
  select(-c('ProfileID', 'ModelID')) %>%
  pivot_longer(names_to='gene', values_to='log2tpm', cols=-'DepMapID') %>%
  pivot_wider(names_from='DepMapID', values_from='log2tpm') %>%
  mutate(gene_symbol=gene %>% gsub(' \\(|\\)|ENSG[0-9]{11}', '', .)) %>%
  mutate(ensembl_gene_id=str_extract(gene, 'ENSG[0-9]{11}')) %>%
  select(gene, gene_symbol, ensembl_gene_id, everything())

tpm.df <- OmicsExpressionAllGenesTPMLogp1Profile %>%
  rename(ProfileID='...1') %>%
  left_join(OmicsDefaultModelProfiles %>% select(ProfileID, ModelID), by='ProfileID') %>%
  filter(!is.na(ModelID)) %>%
  mutate(DepMapID=sub('-', '.', ModelID)) %>%
  select(-c('ProfileID', 'ModelID')) %>%
  pivot_longer(names_to='gene', values_to='log2tpm', cols=-'DepMapID') %>%
  mutate(tpm=(2^log2tpm)-1) %>%
  select(-'log2tpm') %>%
  pivot_wider(names_from='DepMapID', values_from='tpm') %>%
  mutate(gene_symbol=gene %>% gsub(' \\(|\\)|ENSG[0-9]{11}', '', .)) %>%
  mutate(ensembl_gene_id=str_extract(gene, 'ENSG[0-9]{11}')) %>%
  select(gene, gene_symbol, ensembl_gene_id, everything())

first.loc.start <- segment.arm %>%
  filter(chr=='X') %>% 
  pull(loc.start) %>%
  min()

last.end <- segment.arm %>%
  filter(chr=='X') %>%
  pull(end) %>%
  max()

gene.type.flt <- gene.type.lo %>%
  filter(start.hg38 >= first.loc.start & end.hg38 <= last.end)

tpm.df.l <- tpm.df %>%
  filter(ensembl_gene_id %in% gene.type.flt$ensembl_gene_id) %>%
  pivot_longer(names_to='DepMapID', values_to='tpm', cols=-c('gene', 'gene_symbol', 'ensembl_gene_id')) %>%
  filter(DepMapID %in% samples.to.be.analyzed)

tpm.medians <- tpm.df.l %>%
  left_join(sif %>% select(DepMap_ID, sex), by=c('DepMapID'='DepMap_ID')) %>%
  group_by(gene_symbol, ensembl_gene_id, sex) %>%
  summarize(median=median(tpm)) %>%
  ungroup() %>%
  mutate(log2median=log2(median + 1))

tpm.threshold <- 1

genes.to.be.analyzed <- tpm.medians %>%
  filter(!is.na(sex)) %>%
  filter(median >= tpm.threshold) %>%
  pull(ensembl_gene_id) %>%
  unique()

tpm.df.nest <- tpm.df.l %>%
  filter(ensembl_gene_id %in% genes.to.be.analyzed) %>%
  left_join(gene.type.flt %>% select(ensembl_gene_id, Combined_XCI_status, start.hg38, end.hg38), by='ensembl_gene_id') %>%
  rename(type=Combined_XCI_status) %>%
  filter(type %in% c('inactive', 'escape')) %>%
  mutate(id=ensembl_gene_id, gene=gene_symbol) %>%
  group_by(id, gene) %>%
  nest()


ensemblids <- tpm.df.nest$id
female.standard <- 0
male.standard <- -0.95
margin <- 0.2
purity.threshold <- 0
sample.num.threshold <- 3
exp.outlier.threshold <-10

df <- tpm.df.nest[11, 'data'] %>% unnest()

## 1 < CN < 2
reg.analysis <- function(df) {
  df.i <- df %>%
    select(gene_symbol, ensembl_gene_id, type) %>%
    distinct()

  i <- which(ensemblids==df.i$ensembl_gene_id)
  print(paste(i, df.i$gene_symbol))

  df.i <- df.i %>%
    mutate(i=i) %>%
    select(c('i', everything()))

  start.position <- df$start.hg38 %>% unique()
  end.position <- df$end.hg38 %>% unique()

  seg.x <- segment.arm %>%
    rename(Gender=sex) %>%
    filter(chr=='X') %>%
    filter(loc.start <= start.position & end >= end.position)

  if (nrow(seg.x)!=0) {
    cn.tpm.female <- seg.x %>%
      inner_join(df, by='DepMapID') %>%
      left_join(cbio %>% select(DepMapID, Purity, Ploidy, GenomeDoublings, ploidy.class), by='DepMapID') %>%
      filter(Gender=='Female') %>%
      filter(!is.na(ploidy.class)) %>%
      filter(ploidy.class=='Diploid') %>%
      filter(tpm!=0) %>%
      mutate(RCN=2^seg.mean) %>%
      mutate(cn=case_when(RCN < 2^(female.standard - margin) ~ 'cn1',
                          RCN >= 2^(female.standard - margin) & RCN <= 2^(female.standard + margin) ~ 'cn2',
                          RCN > 2^(female.standard + margin) ~ 'cn3',
                          TRUE ~ 'other')) %>%
      group_by(lineage) %>%
      mutate(mean.exp.cn2=median(tpm[cn=='cn2'])) %>%
      filter(mean.exp.cn2!=0) %>%
      mutate(rel.exp=tpm/mean.exp.cn2) %>%
      filter(!is.na(rel.exp)) %>%
      filter(rel.exp < exp.outlier.threshold) %>%
      ungroup()

    cn.tpm.male <- seg.x %>%
      inner_join(df, by='DepMapID') %>%
      left_join(cbio %>% select(DepMapID, Purity, Ploidy, GenomeDoublings, ploidy.class), by='DepMapID') %>%
      filter(Gender=='Male') %>%
      filter(!is.na(ploidy.class)) %>%
      filter(ploidy.class=='Diploid') %>%
      filter(tpm!=0) %>%
      mutate(RCN=2^seg.mean) %>%
      mutate(cn=case_when(RCN >= 2^(male.standard - margin) & RCN <= 2^(male.standard + margin) ~ 'cn1',
                          RCN > 2^(male.standard + margin) & RCN <= 2^(female.standard + margin) ~ 'cn2',
                          RCN > 2^(female.standard + margin) ~ 'cn3',
                          TRUE ~ 'other')) %>%
      group_by(lineage) %>%
      mutate(mean.exp.cn1=median(tpm[cn=='cn1'])) %>%
      filter(mean.exp.cn1!=0) %>%
      mutate(rel.exp=tpm/mean.exp.cn1) %>%
      filter(!is.na(rel.exp)) %>%
      filter(rel.exp < exp.outlier.threshold) %>%
      ungroup()

    cn.tpm <- cn.tpm.female %>%
      bind_rows(cn.tpm.male)

    ## 1 < CN < 2
    female.cn1.samples <- cn.tpm.female %>%
      filter(Purity >= purity.threshold) %>%
      filter(cn=='cn1')

    male.cn2.samples <- cn.tpm.male %>%
      filter(Purity >= purity.threshold) %>%
      filter(cn=='cn2')

    cn.tpm.flt <- cn.tpm %>%
      filter(Purity >= purity.threshold)

    if (nrow(female.cn1.samples) >= sample.num.threshold) {
      reg.female1 <- lm(rel.exp ~ RCN, data=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^(female.standard)))
      slope.female1 <- summary(reg.female1)$coefficients['RCN','Estimate']
      pval.female1 <- summary(reg.female1)$coefficients['RCN','Pr(>|t|)']
    } else {
      slope.female1 <- NA
      pval.female1 <- NA
    }

    if (nrow(male.cn2.samples) >= sample.num.threshold) {
      reg.male1 <- lm(rel.exp ~ RCN, data=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard)))
      slope.male1 <- summary(reg.male1)$coefficients['RCN','Estimate']
      pval.male1 <- summary(reg.male1)$coefficients['RCN','Pr(>|t|)']
    } else {
      slope.male1 <- NA
      pval.male1 <- NA
    }

    df.1to2.female <- df.i %>%
      mutate(nrow=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^female.standard) %>% nrow(),
              samples=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^female.standard) %>% pull(DepMapID) %>% unique() %>% length(),
              tpm.median=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^female.standard) %>% pull(tpm) %>% median(),
              slope=slope.female1,
              pval=pval.female1,
              cn='1 < CN < 2',
              Gender='Female')

    df.1to2.male <- df.i %>%
      mutate(nrow=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard)) %>% nrow(),
              samples=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard)) %>% pull(DepMapID) %>% unique() %>% length(),
              tpm.median=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard)) %>% pull(tpm) %>% median(),
              slope=slope.male1,
              pval=pval.male1,
              cn='1 < CN < 2',
              Gender='Male')

    df.1to2 <- df.1to2.female %>%
      bind_rows(df.1to2.male)

    ## 2 < CN
    female.cn3.samples <- cn.tpm.female %>%
      filter(Purity >= purity.threshold) %>%
      filter(cn=='cn3')

    male.cn3.samples <- cn.tpm.male %>%
      filter(Purity >= purity.threshold) %>%
      filter(cn=='cn3')

    if (nrow(female.cn3.samples) >= sample.num.threshold) {
      reg.female2 <- lm(rel.exp ~ RCN, data=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(female.standard)))
      slope.female2 <- summary(reg.female2)$coefficients['RCN','Estimate']
      pval.female2 <- summary(reg.female2)$coefficients['RCN','Pr(>|t|)']
    } else {
      slope.female2 <- NA
      pval.female2 <- NA
    }

    if (nrow(male.cn3.samples) >= sample.num.threshold) {
      reg.male2 <- lm(rel.exp ~ RCN, data=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(female.standard)))
      slope.male2 <- summary(reg.male2)$coefficients['RCN','Estimate']
      pval.male2 <- summary(reg.male2)$coefficients['RCN','Pr(>|t|)']
    } else {
      slope.male2 <- NA
      pval.male2 <- NA
    }

    df.2plus.female <- df.i %>%
      mutate(nrow=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(female.standard)) %>% nrow(),
              samples=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(female.standard)) %>% pull(DepMapID) %>% unique() %>% length(),
              tpm.median=cn.tpm.flt %>% filter(Gender=='Female' & RCN >= 2^(female.standard)) %>% pull(tpm) %>% median(),
              slope=slope.female2,
              pval=pval.female2,
              cn='2 < CN',
              Gender='Female')

    df.2plus.male <- df.i %>%
      mutate(nrow=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(female.standard)) %>% nrow(),
              samples=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(female.standard)) %>% pull(DepMapID) %>% unique() %>% length(),
              tpm.median=cn.tpm.flt %>% filter(Gender=='Male' & RCN >= 2^(female.standard)) %>% pull(tpm) %>% median(),
              slope=slope.male2,
              pval=pval.male2,
              cn='2 < CN',
              Gender='Male')

    df.2plus <- df.2plus.female %>%
      bind_rows(df.2plus.male)

    df.slope <- df.1to2 %>%
      bind_rows(df.2plus)

    ## plot
    female.dotted.lines <- c(2^(female.standard - margin), 2^(female.standard + margin))
    male.dotted.lines <- c(2^(male.standard - margin), 2^(male.standard + margin))

    g1 <- ggplot(cn.tpm %>% left_join(df.1to2, by=c('gene_symbol', 'Gender')), aes(x=RCN, y=rel.exp)) +
      geom_vline(xintercept=2^(female.standard), col='black', linetype='dashed') +
      geom_vline(xintercept=2^(male.standard), col='blue', linetype='dashed') +
      geom_vline(xintercept=1.5, col='red', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Female'), aes(xintercept=female.dotted.lines[1]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Female'), aes(xintercept=female.dotted.lines[2]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Male'), aes(xintercept=male.dotted.lines[1]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Male'), aes(xintercept=male.dotted.lines[2]), col='darkgray', linetype='dashed') +
      geom_point(data=. %>% filter(!((Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^(female.standard) & Purity >= purity.threshold)|(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard) & Purity >= purity.threshold))), col='gray', size=2) +
      geom_point(data=. %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^(female.standard) & Purity >= purity.threshold), aes(col=lineage), size=2) +
      geom_point(data=. %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard) & Purity >= purity.threshold), aes(col=lineage), size=2) +
      geom_smooth(data=. %>% filter(Gender=='Female' & RCN >= 2^(-1.2) & RCN <= 2^(female.standard) & Purity >= purity.threshold), method='lm') +
      geom_smooth(data=. %>% filter(Gender=='Male' & RCN >= 2^(male.standard - margin) & RCN <= 2^(female.standard) & Purity >= purity.threshold), method='lm') +
      geom_text(aes(x=-Inf, y=Inf, hjust=-0.05, vjust=1.5, label=paste0('slope = ', slope)), size=10) +
      scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3), limits=c(0, NA)) +
      facet_wrap(~Gender, nrow=2) +
      {if (nrow(female.cn1.samples) >= sample.num.threshold) geom_rect(data=. %>% filter(Gender=='Female'), fill=NA, color='green', linewidth=3, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)} +
      {if (nrow(male.cn2.samples) >= sample.num.threshold) geom_rect(data=. %>% filter(Gender=='Male'), fill=NA, color='green', linewidth=3, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)} +
      labs(title=paste(df.i$gene_symbol, df.i$ensembl_gene_id, '\n', df.i$type, '\n(1', '\u2264', 'CN', '\u2264', '2)'), y='Relative Gene Exp') +
      theme_bw(base_size=30)

    g2 <- ggplot(cn.tpm %>% left_join(df.2plus, by=c('gene_symbol', 'Gender')), aes(x=RCN, y=rel.exp)) +
      geom_vline(xintercept=2^(female.standard), col='black', linetype='dashed') +
      geom_vline(xintercept=2^(male.standard), col='blue', linetype='dashed') +
      geom_vline(xintercept=1.5, col='red', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Female'), aes(xintercept=female.dotted.lines[1]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Female'), aes(xintercept=female.dotted.lines[2]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Male'), aes(xintercept=male.dotted.lines[1]), col='darkgray', linetype='dashed') +
      geom_vline(data=. %>% filter(Gender=='Male'), aes(xintercept=male.dotted.lines[2]), col='darkgray', linetype='dashed') +
      geom_point(data=. %>% filter(!(RCN >= 2^(female.standard) & Purity >= purity.threshold)), col='gray', size=2) +
      geom_point(data=. %>% filter(RCN >= 2^(female.standard) & Purity >= purity.threshold), aes(col=lineage), size=2) +
      geom_smooth(data=. %>% filter(Gender=='Female' & RCN >= 2^(female.standard) & Purity >= purity.threshold), method='lm') +
      geom_smooth(data=. %>% filter(Gender=='Male' & RCN >= 2^(female.standard) & Purity >= purity.threshold), method='lm') +
      geom_text(aes(x=-Inf, y=Inf, hjust=-0.05, vjust=1.5, label=paste0('slope = ', slope)), size=10) +
      scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3), limits=c(0, NA)) +
      facet_wrap(~Gender, nrow=2) +
      {if (nrow(female.cn3.samples) >= sample.num.threshold) geom_rect(data=. %>% filter(Gender=='Female'), fill=NA, color='green', linewidth=3, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)} +
      {if (nrow(male.cn3.samples) >= sample.num.threshold) geom_rect(data=. %>% filter(Gender=='Male'), fill=NA, color='green', linewidth=3, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)} +
      labs(title=paste(df.i$gene_symbol, df.i$ensembl_gene_id, '\n', df.i$type, '\n(CN', '\u2265', '2)'), y='Relative gene exp') +
      theme_bw(base_size=30)

    g <- ggpubr::ggarrange(g1, g2, nrow=1, common.legend=TRUE, legend='bottom')
    ggsave(g, file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression/plot/chrX', paste0(i, '_', df.i$gene_symbol, '_', df.i$type, '.png')), dpi=100, width=18, height=18)
  }
  return(df.slope)
}

lm.df.list <- parallel::mclapply(tpm.df.nest$data, reg.analysis, mc.cores=parallel::detectCores() - 2)
lm.df <- lm.df.list %>%
  bind_rows()
saveRDS(lm.df, file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'lm.df.rds'), compress=FALSE)

lm.df.l <- lm.df %>%
  filter(tpm.median > 1 | gene_symbol=='XIST') %>%
  filter(nrow >= 10) %>%
  filter(!is.na(slope)) %>%
  mutate(slope=slope/2) %>%
  mutate(type=case_when(type=='inactive' ~ 'Inactive genes', type=='escape' ~ 'Escape genes')) %>%
  mutate(type=factor(.$type, levels=c('Inactive genes', 'Escape genes'))) %>%
  group_by(Gender, type, gene_symbol) %>%
  mutate(n=n()) %>%
  filter((Gender=='Female' & n==2) | Gender=='Male') %>%
  select(-n) %>%
  ungroup() %>%
  group_by(Gender, type, cn) %>%
  mutate(n=n()) %>%
  ungroup()

## Diploid, calculation method1
stat.test.cn <- lm.df.l %>%
  filter(Gender=='Female') %>%
  # filter(gene.symbol!='XIST') %>%
  group_by(Gender, type) %>%
  rstatix::wilcox_test(slope ~ cn, paired=TRUE)

stat.test.ie <- lm.df.l %>%
  filter(Gender=='Female' | (Gender=='Male' & cn=='1 < CN < 2')) %>%
  group_by(Gender, cn) %>%
  rstatix::wilcox_test(slope ~ type)

stat.test.fm <- lm.df.l %>%
  filter(cn=='1 < CN < 2') %>%
  group_by(type, cn) %>%
  rstatix::wilcox_test(slope ~ Gender)

stat.test <- stat.test.cn %>%
  bind_rows(stat.test.ie) %>%
  bind_rows(stat.test.fm) %>%
  mutate(p_bonferroni=p.adjust(p, method='bonferroni')) %>%
  select(type, Gender, cn, everything())

g <- ggplot(lm.df.l, aes(x=cn, y=slope)) +
  geom_hline(yintercept=0, col='blue', linetype='dashed') +
  geom_hline(yintercept=1, col='red', linetype='dashed') +
  geom_violin(position=position_dodge(0.9), show.legend=FALSE) +
  geom_point(shape=21, position=ggbeeswarm::position_beeswarm(dodge.width=0.9), alpha=0.5) +
  geom_boxplot(position=position_dodge(0.9), width=0.1, col='red', alpha=0.1, outlier.shape=NA, show.legend=FALSE) +
  geom_text(aes(x=Inf, y=-Inf, label=paste0(n, ' genes')), vjust=-1, hjust=1, size=5) +
  scale_x_discrete(labels=c('1 < CN < 2', '2 < CN')) +
  coord_cartesian(ylim=c(-2.5, NA)) +
  ggrepel::geom_text_repel(data=. %>% filter(gene_symbol=='XIST'), aes(label=gene_symbol, col=type), nudge_x=0.3, box.padding=1, size=5, show.legend=FALSE) +
  lemon::facet_rep_grid(Gender~type, switch='y') +
  ggpubr::stat_pvalue_manual(stat.test.cn %>% filter(p < 0.05) %>% rstatix::add_xy_position(x='cn'), label='P = {p}', col='red', size=5) +
  ggpubr::stat_pvalue_manual(stat.test.cn %>% filter(p >= 0.05) %>% rstatix::add_xy_position(x='cn'), label='P = {p}', size=5) +
  labs(x='CN', y='Slope') +
  theme_classic(base_size=20) +
  theme(legend.title=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.placement='outside') +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'Fig5c.png'), dpi=100, width=9, height=8)
ggsave(g, file=here('10_chrX_SCNA_vs_gene_expression/output/02_CCLE_chrX_SCNA_vs_gene_expression', 'Fig5c.pdf'), width=9, height=8, useDingbats=TRUE)

## Perform one sample T-test to see if average of slopes is significantly different from 1
female.inactive.cn1to2 <- t.test(lm.df.l %>% filter(Gender=='Female' & type=='Inactive genes' & cn=='1 < CN < 2') %>% pull(slope), mu=0)
female.inactive.cn2 <- t.test(lm.df.l %>% filter(Gender=='Female' & type=='Inactive genes' & cn=='2 < CN') %>% pull(slope), mu=0.5)
female.escape.cn1to2 <- t.test(lm.df.l %>% filter(Gender=='Female' & type=='Escape genes' & cn=='1 < CN < 2') %>% pull(slope), mu=0)
female.escape.cn2 <- t.test(lm.df.l %>% filter(Gender=='Female' & type=='Escape genes' & cn=='2 < CN') %>% pull(slope), mu=0.5)
male.inactive.cn1to2 <- t.test(lm.df.l %>% filter(Gender=='Male' & type=='Inactive genes' & cn=='1 < CN < 2') %>% pull(slope), mu=1)
male.escape.cn1to2 <- t.test(lm.df.l %>% filter(Gender=='Male' & type=='Escape genes' & cn=='1 < CN < 2') %>% pull(slope), mu=1)

one.sample.t.test.table <- data.frame(Gender=c('female', 'female', 'female', 'female', 'male', 'male'),
                                      type=c('inactive', 'inactive', 'escape', 'escape', 'inactive', 'escape'),
                                      cn=c('1 < CN < 2', '2 < CN', '1 < CN < 2', '2 < CN', '1 < CN < 2', '2 < CN'),
                                      mean.target=c(0, 0.5, 0, 0.5, 1, 1),
                                      pval=c(female.inactive.cn1to2$p.value, female.inactive.cn2$p.value, female.escape.cn1to2$p.value, female.escape.cn2$p.value, male.inactive.cn1to2$p.value, male.escape.cn1to2$p.value)) %>%
                            mutate(signif=case_when(pval < 0.05 ~ '*', TRUE ~ '')) %>%
                            mutate(pval_bonnferroni=p.adjust(.$pval, method='bonferroni'))

lm.df.l %>%
  group_by(Gender, type, cn) %>%
  summarize(n=n(), slope.mean=mean(slope), slope.median=median(slope))
