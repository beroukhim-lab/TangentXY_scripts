library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

sample.amp.del <- readRDS(file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'sample.amp.del.rds'))

cosmic.census <- read.delim(here('08_TCGA_SCNA_vs_gene_mutation/data', 'Census_allFri Feb 16 01_41_07 2024.tsv'))

mut <- maftools::read.maf(here('08_TCGA_SCNA_vs_gene_mutation/data', 'mc3.v0.2.8.PUBLIC.maf.gz'))

mut.data <- mut@data %>%
  mutate(tcga.id=str_extract(string=Tumor_Sample_Barcode, pattern='^TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9]{2}')) %>%
  mutate(tcga.id=gsub('-', '.', tcga.id)) %>%
  mutate(tcga.id=stringi::stri_replace_last_regex(str=tcga.id, pattern='\\.', replacement='_')) %>%
  separate(col=tcga.id, into=c('TCGA.ID', 'type.code'), sep='_') %>%
  mutate(type=case_when(type.code=='01' ~ 'TP',
                        type.code=='02' ~ 'TR',
                        type.code=='03' ~ 'TB',
                        type.code=='05' ~ 'TAP',
                        type.code=='06' ~ 'TM')) %>%
  left_join(sif, by=c('TCGA.ID', 'type'))
saveRDS(mut.data, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'mut.data.rds'), compress=FALSE)

chrx.karyo <- sample.amp.del %>%
  filter(chr=='X') %>%
  select(SampleID, karyo.class, TCGA.ID, Gender, project) %>%
  mutate(chrX.class=case_when(karyo.class=='No.Alt' ~ 'WT', karyo.class %in% c('Whole.Amp', 'Arm.Amp') ~ 'Amp', karyo.class %in% c('Whole.Del', 'Arm.Del') ~ 'Del')) %>%
  filter(!is.na(chrX.class)) %>%
  mutate(chrX.class=factor(.$chrX.class, levels=c('WT', 'Amp', 'Del')))

chry.karyo <- sample.amp.del %>%
  filter(chr=='Y') %>%
  select(SampleID, karyo.class, TCGA.ID, Gender, project) %>%
  mutate(chrY.class=case_when(karyo.class=='No.Alt' ~ 'WT', karyo.class %in% c('Whole.Amp', 'Arm.Amp') ~ 'Amp', karyo.class %in% c('Whole.Del', 'Arm.Del') ~ 'Del')) %>%
  filter(!is.na(chrY.class)) %>%
  mutate(chrY.class=factor(.$chrY.class, levels=c('WT', 'Amp', 'Del')))

genes <- mut.data$Hugo_Symbol %>% unique()


## Original functions to run Fisher's exact test on a list of contingency table
fisher.test.on.cont.table <- function(cont.table) {
  if (dim(cont.table)[1]==2 & dim(cont.table)[2]==2) {
    fisher.result <- fisher.test(cont.table)

    fisher.df <- data.frame(wt.wt=cont.table[1,1],
      cna.wt=cont.table[1,2],
      wt.mut=cont.table[2,1],
      cna.mut=cont.table[2,2],
      total.mut.num=cont.table[2,] %>% sum(),
      odds.ratio=fisher.result$estimate,
      pval=fisher.result$p.value)
  } else {
    fisher.df <- data.frame(pval=NA)  
  }
  return(fisher.df)
}


## Fisher's exact test on each tumor type
tumor.types <- sif$project %>% unique()

mut.fisher.each.tumor <- function(gene) {
  num <- match(gene, genes)
  print(paste(num, gene))

  chr <- mut.data %>%
    filter(Hugo_Symbol==gene) %>%
    pull(Chromosome) %>%
    unique()

  mut.pts <- mut.data %>%
    # filter(Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site', 'Translation_Start_Site')) %>%
    filter(Hugo_Symbol==gene) %>%
    filter(!is.na(SampleID)) %>%
    pull(SampleID) %>%
    unique()

  for (i in 1:length(tumor.types)) {
    tumor.type.i <- tumor.types[i]

    chrx.wt.vs.amp.cont.table <- chrx.amp %>%
      filter(project==tumor.type.i) %>%
      mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
      mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut'))) %>%
      droplevels() %>%
      select(mut.class, chrX.class) %>%
      table()

    chrx.wt.vs.del.cont.table <- chrx.del %>%
      filter(project==tumor.type.i) %>%
      filter(Gender=='Female') %>%
      mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
      mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut'))) %>%
      droplevels() %>%
      select(mut.class, chrX.class) %>%
      table()

    chry.wt.vs.amp.cont.table <- chry.amp %>%
      filter(project==tumor.type.i) %>%
      mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
      mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut'))) %>%
      droplevels() %>%
      select(mut.class, chrY.class) %>%
      table()

    chry.wt.vs.del.cont.table <- chry.del %>%
      filter(project==tumor.type.i) %>%
      mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
      mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut'))) %>%
      droplevels() %>%
      select(mut.class, chrY.class) %>%
      table()

    cont.table.list <- vector(mode='list', length=4)
    cont.table.list[[1]] <- chrx.wt.vs.amp.cont.table
    cont.table.list[[2]] <- chrx.wt.vs.del.cont.table
    cont.table.list[[3]] <- chry.wt.vs.amp.cont.table
    cont.table.list[[4]] <- chry.wt.vs.del.cont.table
    names(cont.table.list) <- c('ChrX_Amp', 'ChrX_Del', 'ChrY_Amp', 'ChrY_Del')

    result.i <- lapply(cont.table.list, fisher.test.on.cont.table) %>%
      bind_rows() %>%
      mutate(i=num) %>%
      mutate(tumor.type=tumor.type.i) %>%
      mutate(gene.symbol=gene) %>%
      mutate(chromosome=chr) %>%
      mutate(type=names(cont.table.list)) %>%
      select(c('i', 'tumor.type', 'gene.symbol', 'chromosome', everything()))

    if (i==1) {
      result <- result.i
    } else {
      result <- result %>% bind_rows(result.i)
    }
  }
  return(result)
}

fisher.df.each.tumor.list <- parallel::mclapply(genes, mut.fisher.each.tumor, mc.cores=parallel::detectCores() - 2)
fisher.df.each.tumor <- fisher.df.each.tumor.list %>%
  bind_rows() %>%
  mutate(cosmic.census=case_when(gene.symbol %in% cosmic.census$Gene.Symbol ~ TRUE, TRUE ~ FALSE))
saveRDS(fisher.df.each.tumor, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'fisher.df.each.tumor.rds'), compress=FALSE)

fisher.df.combine.p <- fisher.df.each.tumor %>%
  filter(!is.na(pval)) %>%
  mutate(pval=case_when(pval > 1 ~ 1, TRUE ~ pval)) %>%
  group_by(gene.symbol, chromosome, type) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  group_by(gene.symbol, chromosome, type, n) %>%
  summarize(comb.p=ifelse(n > 1, metap::sumlog(pval)$p, pval)) %>%
  distinct()

fisher.df.or.nest <- fisher.df.each.tumor %>%
  filter(!is.na(odds.ratio)) %>%
  metafor::escalc(measure='OR', ai=wt.wt, bi=cna.wt, ci=wt.mut, di=cna.mut, data=.) %>%
  group_by(gene.symbol, chromosome, type) %>%
  nest()

make.comb.or <- function(df) {
  log.comb.or <- tryCatch(
    metafor::rma(yi=yi, vi=vi, data=df) %>% .$b %>% as.numeric,
    error = function(e) return(NA)
  )
  return(log.comb.or)
}

comb.or.df <- fisher.df.or.nest %>%
  mutate(log.comb.or=map_dbl(data, ~make.comb.or(.))) %>%
  mutate(comb.or=exp(log.comb.or)) %>%
  mutate(log2.comb.or=log2(comb.or)) %>%
  mutate(total.wt.wt=map_dbl(data, ~sum(.$wt.wt))) %>%
  mutate(total.cna.wt=map_dbl(data, ~sum(.$cna.wt))) %>%
  mutate(total.wt.mut=map_dbl(data, ~sum(.$wt.mut))) %>%
  mutate(total.cna.mut=map_dbl(data, ~sum(.$cna.mut))) %>%
  mutate(total.mut.num=map_dbl(data, ~sum(.$total.mut.num))) %>%
  select(-data)
saveRDS(comb.or.df, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'comb.or.df.rds'), compress=FALSE)

fisher.df.anno <- fisher.df.combine.p %>%
  left_join(comb.or.df, by=c('gene.symbol', 'chromosome', 'type')) %>%
  filter(!(type=='ChrX_Del' & chromosome=='X')) %>%
  filter(!(type=='ChrY_Del' & chromosome=='Y')) %>%
  filter(total.mut.num >= 50) %>%
  ungroup() %>%
  mutate(qval=p.adjust(comb.p, method='BH')) %>%
  mutate(signif=case_when(qval < 0.1 ~ 'Q < 0.1')) %>%
  group_by(type) %>%
  arrange(comb.p) %>%
  mutate(rank=1:n()) %>%
  separate(col=type, into=c('chr', 'alt'), sep='_') %>%
  mutate(cosmic.census=case_when(gene.symbol %in% cosmic.census$Gene.Symbol ~ TRUE, TRUE ~ FALSE))

fisher.df.anno %>% filter(chr=='ChrX' & alt=='Amp') %>% arrange(pval) %>% head()
fisher.df.anno %>% filter(chr=='ChrX' & alt=='Del') %>% arrange(pval) %>% head()
fisher.df.anno %>% filter(chr=='ChrY' & alt=='Amp') %>% arrange(pval) %>% head()
fisher.df.anno %>% filter(chr=='ChrY' & alt=='Del') %>% arrange(pval) %>% head()
fisher.df.anno %>% filter(chr=='ChrY' & alt=='Del') %>% arrange(pval) %>% mutate(n=1:n()) %>% filter(gene.symbol=='TP53')

g <- ggplot(fisher.df.anno %>% filter(!is.na(log.comb.or)), aes(x=log2.comb.or, y=comb.p)) +
  geom_point(data=. %>% filter(is.na(signif)), aes(col=signif), size=2) +
  # geom_point(data=. %>% filter(cosmic.census==TRUE), col='limegreen', size=2) +
  geom_point(data=. %>% filter(!is.na(signif)), aes(col=signif), size=2) +
  geom_vline(xintercept=0, linetype='dashed') +
  ggrepel::geom_text_repel(data=. %>% filter(!is.na(signif) & cosmic.census==TRUE), aes(label=gene.symbol), size=5, col='black') +
  ggrepel::geom_text_repel(data=. %>% filter(chr=='ChrY', gene.symbol=='TP53'), aes(label=gene.symbol), size=5, col='grey70') +
  ggh4x::facet_grid2(alt ~ chr, scales='free', independent='all', switch='y') +
  ggh4x::facetted_pos_scales(y=list(
    chr=='ChrX' & alt=='Amp' ~ scale_y_continuous(trans=scales::compose_trans('log10', 'reverse'), labels=scales::label_log(), breaks=c(10^(0), 10^(-5), 10^(-10), 10^(-15)), limits=c(NA, 10^(-15))),
    chr=='ChrX' & alt=='Del' ~ scale_y_continuous(trans=scales::compose_trans('log10', 'reverse'), labels=scales::label_log(), breaks=c(10^(0), 10^(-10), 10^(-20), 10^(-30), 10^(-40), 10^(-50)), limits=c(NA, 10^(-50))),
    chr=='ChrY' & alt=='Amp' ~ scale_y_continuous(trans=scales::compose_trans('log10', 'reverse'), labels=scales::label_log(), breaks=c(10^(0), 10^(-1))),
    chr=='ChrY' & alt=='Del' ~ scale_y_continuous(trans=scales::compose_trans('log10', 'reverse'), labels=scales::label_log(), breaks=c(10^(0), 10^(-1), 10^(-2), 10^(-3))))
  ) +
  # scale_y_continuous(trans=scales::compose_trans('log10', 'reverse'), labels=scales::label_log()) +
  scale_x_continuous(breaks=c(-1, 0, 1, 2)) +
  scale_color_manual(values=c('Q < 0.1'='red'), na.value='grey90', limits='Q < 0.1') +
  labs(x=expression(paste({log[2]}, '[Odds ratio]', sep='')), y='p-value') +
  theme_classic(base_size=20) +
  theme(legend.title=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.placement='outside') +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'Fig3e.png'), dpi=100, width=10, height=8)
ggsave(g, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'Fig3e.pdf'), width=10, height=8, useDingbats=TRUE)

## Make supplementary tables
chrx.amp <- fisher.df.anno %>%
  filter(chr=='ChrX' & alt=='Amp') %>%
  arrange(qval, comb.p) %>%
  select(gene.symbol, chromosome, total.wt.wt, total.cna.wt, total.wt.mut, total.cna.mut, total.mut.num, log2.comb.or, comb.p, qval, cosmic.census) %>%
  rename(combined.pval=comb.p, log2.combined.OR=log2.comb.or)

chrx.del <- fisher.df.anno %>%
  filter(chr=='ChrX' & alt=='Del') %>%
  arrange(qval, comb.p) %>%
  select(gene.symbol, chromosome, total.wt.wt, total.cna.wt, total.wt.mut, total.cna.mut, total.mut.num, log2.comb.or, comb.p, qval, cosmic.census) %>%
  rename(combined.pval=comb.p, log2.combined.OR=log2.comb.or)

chry.amp <- fisher.df.anno %>%
  filter(chr=='ChrY' & alt=='Amp') %>%
  arrange(qval, comb.p) %>%
  select(gene.symbol, chromosome, total.wt.wt, total.cna.wt, total.wt.mut, total.cna.mut, total.mut.num, log2.comb.or, comb.p, qval, cosmic.census) %>%
  rename(combined.pval=comb.p, log2.combined.OR=log2.comb.or)

chry.del <- fisher.df.anno %>%
  filter(chr=='ChrY' & alt=='Del') %>%
  arrange(qval, comb.p) %>%
  select(gene.symbol, chromosome, total.wt.wt, total.cna.wt, total.wt.mut, total.cna.mut, total.mut.num, log2.comb.or, comb.p, qval, cosmic.census) %>%
  rename(combined.pval=comb.p, log2.combined.OR=log2.comb.or)

sheets <- list('ChrX_Amp'=chrx.amp, 'ChrX_Del'=chrx.del, 'ChrY_Amp'=chry.amp, 'ChrY_Del'=chry.del)
openxlsx::write.xlsx(sheets, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'TableS1.xlsx'))
