library(tidyverse)
library(here)

library(DNAcopy)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

## Run CBS on TangentXY-normalized samples
Tn <- readRDS(file=here('06_CCLE_TangentXY/output/03_TangentXY', 'Tn.rds'))

sample.names <- colnames(Tn)

probes <- Tn %>%
  as.data.frame() %>%
  rownames_to_column('locus') %>%
  select(locus) %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  separate(col=pos, into=c('start', 'end'), sep='-') %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique() %>% gtools::mixedsort()))

segment.smoothed.CNA.obj.pr <- readRDS(file=here('06_CCLE_TangentXY/output/07_CBS', 'CBS_PrototypeTangent.rds'))
segment.smoothed.CNA.obj.xy <- readRDS(file=here('06_CCLE_TangentXY/output/07_CBS', 'CBS_TangentXY.rds'))


## Make figures for publication
samples <- c('ACH.000130', 'ACH.000462')

ggplotSample <- function(sample.id) {
  gender <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(Sex)

  print(paste(sample.id, gender))

  ## Prototype Tangent
  points.pr <- segment.smoothed.CNA.obj.pr$data[[sample.id]] %>%
    as.data.frame() %>%
    setNames('log2rcn') %>%
    mutate(loci=row_number()) %>%
    bind_cols(probes) %>%
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even'))

  segments.pr <- segment.smoothed.CNA.obj.pr$output %>%
    filter(ID==sample.id) %>%
    left_join(points %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.pr2 <- points.pr %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.pr %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.pr2 <- segments.pr %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.pr %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.pr <- ggplot(points.pr2, aes(x=loci, y=log2rcn)) +
    geom_point(aes(col=chr.class), show.legend=FALSE) +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.pr2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    scale_color_manual(values=c('odd'='limegreen', 'even'='black')) +
    coord_cartesian(ylim=c(-2.5, 2.5)) +
    labs(title='Tangent', x='Genomic position', y=expression(paste({log[2]}, '[Relative CN]', sep=''))) +
    theme_bw(base_size=20) +
    theme(panel.grid=element_blank()) +
    theme(axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.line.x=element_line(linewidth=0.5)) +
    theme(axis.line.y=element_line(linewidth=0.5)) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggforce::facet_zoom(x=chr %in% c('chrX', 'chrY'), zoom.size=1, shrink=FALSE) +
    theme(zoom=element_rect(linewidth=0.5, linetype='dashed')) +
    theme(zoom.y = element_blank(), validate = FALSE)

  ## TangentXY
  points.xy <- segment.smoothed.CNA.obj.xy$data[[sample.id]] %>%
    as.data.frame() %>%
    setNames('log2rcn') %>%
    mutate(loci=row_number()) %>%
    bind_cols(probes) %>%
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even'))

  segments.xy <- segment.smoothed.CNA.obj.xy$output %>%
    filter(ID==sample.id) %>%
    left_join(points %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.xy2 <- points.xy %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.xy %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.xy2 <- segments.xy %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.xy %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.xy <- ggplot(points.xy2, aes(x=loci, y=log2rcn)) +
    geom_point(aes(col=chr.class), show.legend=FALSE) +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.xy2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    scale_color_manual(values=c('odd'='limegreen', 'even'='black')) +
    coord_cartesian(ylim=c(-2.5, 2.5)) +
    labs(title='TangentXY', x='Genomic position') +
    theme_bw(base_size=20) +
    theme(panel.grid=element_blank()) +
    theme(axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    theme(axis.line.x=element_line(linewidth=0.5)) +
    theme(axis.line.y=element_line(linewidth=0.5)) +
    theme(plot.title=element_text(hjust=0.5)) +
    ggforce::facet_zoom(x=chr %in% c('chrX', 'chrY'), zoom.size=1, shrink=FALSE) +
    theme(zoom=element_rect(linewidth=0.5, linetype='dashed')) +
    theme(zoom.y = element_blank(), validate = FALSE)

  g <- ggpubr::ggarrange(g.pr, g.xy, nrow=1, align='hv')

  ggsave(g, file=here('06_CCLE_TangentXY/output/08_make_CBS_plots', paste0(sample.id, '_', gender, '.pdf')), width=10, height=8, useDingbats=TRUE)
  ggsave(g, file=here('06_CCLE_TangentXY/output/08_make_CBS_plots', paste0(sample.id, '_', gender, '.png')), dpi=100, width=10, height=8)
}

parallel::mclapply(samples, ggplotSample, mc.cores=parallel::detectCores() - 2)
