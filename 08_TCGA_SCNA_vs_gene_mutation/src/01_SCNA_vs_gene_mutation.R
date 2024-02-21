library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

sample.amp.del <- readRDS(file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'sample.amp.del.rds'))

mut <- maftools::read.maf(here('08_TCGA_SCNA_vs_gene_mutation/data', 'mc3.v0.2.8.PUBLIC.maf.gz'))

cosmic.census <- read.delim(here('08_TCGA_SCNA_vs_gene_mutation/data', 'Census_allFri Feb 16 01_41_07 2024.tsv'))

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


## Fisher's exact test
mut.fisher <- function(gene) {
  num <- match(gene, genes)
  print(paste(num, gene))

  chr <- mut.data %>%
    filter(Hugo_Symbol==gene) %>%
    pull(Chromosome) %>%
    unique()

  mut.pts <- mut.data %>%
    filter(Hugo_Symbol==gene) %>%
    filter(!is.na(SampleID)) %>%
    pull(SampleID) %>%
    unique()

  chrx.karyo.mut <- chrx.karyo %>%
    mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
    mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut')))

  chrx.wt.vs.amp.cont.table <- chrx.karyo.mut %>%
    filter(chrX.class %in% c('WT', 'Amp')) %>%
    droplevels() %>%
    select(mut.class, chrX.class) %>%
    table()

  chrx.wt.vs.del.cont.table <- chrx.karyo.mut %>%
    filter(Gender=='Female') %>%
    filter(chrX.class %in% c('WT', 'Del')) %>%
    droplevels() %>%
    select(mut.class, chrX.class) %>%
    table()

  chry.karyo.mut <- chry.karyo %>%
    mutate(mut.class=case_when(SampleID %in% mut.pts ~ 'Mut', TRUE ~ 'WT')) %>%
    mutate(mut.class=factor(.$mut.class, levels=c('WT', 'Mut')))

  chry.wt.vs.amp.cont.table <- chry.karyo.mut %>%
    filter(chrY.class %in% c('WT', 'Amp')) %>%
    droplevels() %>%
    select(mut.class, chrY.class) %>%
    table()

  chry.wt.vs.del.cont.table <- chry.karyo.mut %>%
    filter(chrY.class %in% c('WT', 'Del')) %>%
    droplevels() %>%
    select(mut.class, chrY.class) %>%
    table()

  cont.table.list <- vector(mode='list', length=4)
  cont.table.list[[1]] <- chrx.wt.vs.amp.cont.table
  cont.table.list[[2]] <- chrx.wt.vs.del.cont.table
  cont.table.list[[3]] <- chry.wt.vs.amp.cont.table
  cont.table.list[[4]] <- chry.wt.vs.del.cont.table
  names(cont.table.list) <- c('ChrX_Amp', 'ChrX_Del', 'ChrY_Amp', 'ChrY_Del')

  for (i in 1:length(cont.table.list)) {
    cont.table <- cont.table.list[[i]]
    if (dim(cont.table)[1]==2 & dim(cont.table)[2]==2) {
      fisher.result <- fisher.test(cont.table)

      fisher.df.i <- data.frame(i=num,
        gene.symbol=gene,
        chromosome=chr,
        wt.wt=cont.table[1,1],
        cna.wt=cont.table[1,2],
        wt.mut=cont.table[2,1],
        cna.mut=cont.table[2,2],
        total.mut.num=cont.table[2,] %>% sum(),
        odds.ratio=fisher.result$estimate,
        pval=fisher.result$p.value,
        type=names(cont.table.list)[i])
    } else {
      fisher.df.i <- data.frame(i=num,
        gene.symbol=gene,
        chromosome=chr,
        type=names(cont.table.list)[i])    
    }
    if (i==1) {
      result <- fisher.df.i
    } else {
      result <- result %>% bind_rows(fisher.df.i)
    }
  }
  return(result)
}

fisher.df.list <- parallel::mclapply(genes, mut.fisher, mc.cores=parallel::detectCores() - 2)
fisher.df <- fisher.df.list %>%
  bind_rows() %>%
  mutate(cosmic.census=case_when(gene.symbol %in% cosmic.census$Gene.Symbol ~ TRUE, TRUE ~ FALSE))
saveRDS(fisher.df, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'fisher.df.rds'), compress=FALSE)

fisher.df.anno <- fisher.df %>%
  filter(!(type=='ChrX_Del' & chromosome=='X')) %>%
  filter(!(type=='ChrY_Del' & chromosome=='Y')) %>%
  filter(total.mut.num >= 50) %>%
  filter(!is.na(pval)) %>%
  mutate(qval=p.adjust(pval, method='BH')) %>%
  mutate(signif=case_when(qval < 0.01 ~ 'Q < 0.01')) %>%
  separate(col=type, into=c('chr', 'alt'), sep='_') %>%
  arrange(desc(qval))

g <- ggplot(fisher.df.anno, aes(x=odds.ratio, y=-log10(pval), col=qval)) +
  geom_point(aes(col=signif), size=2) +
  ggrepel::geom_text_repel(data=. %>% filter(!(alt=='Del' & chr=='ChrX') & qval < 0.01 & cosmic.census==TRUE), aes(label=gene.symbol), col='black', size=4, seed=123456) +
  ggrepel::geom_text_repel(data=. %>% filter(alt=='Del' & chr=='ChrX' & qval < 0.01 & cosmic.census==TRUE), aes(label=gene.symbol), col='black', size=4, force=200, max.overlaps=Inf, seed=123456) +
  ggh4x::facet_grid2(alt ~ chr, scales='free', independent='all', switch='y') +
  scale_color_manual(values=c('Q < 0.01'='red'), na.value='grey90') +
  labs(x='Odds ratio', y=expression(paste({-log[10]}, '[p-value]', sep=''))) +
  theme_classic(base_size=20) +
  theme(legend.title=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.placement='outside') +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'FigS5.png'), dpi=100, width=10, height=8)
ggsave(g, file=here('08_TCGA_SCNA_vs_gene_mutation/output/01_SCNA_vs_gene_mutation', 'FigS5.pdf'), width=10, height=8, useDingbats=TRUE)
