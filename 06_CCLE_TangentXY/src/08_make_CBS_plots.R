library(tidyverse)
library(here)

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
female.samples <- c('ACH.000462', 'ACH.000081', 'ACH.000944', 'ACH.000812', 'ACH.000435')
male.samples <- c('ACH.000130', 'ACH.000065', 'ACH.000156', 'ACH.000158', 'ACH.000342', 'ACH.000492')
samples <- c(female.samples, male.samples)
karyotypes <- c('<2n>X, -X', '<2n>XX', '<2n>X, -X', '<3n>XXX', '<3n>XX, -X', '<2n>XY', '<2n>XY', '<2n>X, -Y', '<3n>XXYY', '<4n>XX, -Y, -Y', '<2n>X, -Y')

ggplotCBS <- function(sample.id) {
  gender <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(Sex)

  print(paste(sample.id, gender))

  i <- grep(sample.id, samples)
  karyotype <- karyotypes[i]
  cell.line.name <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(CellLineName)
  primary.disease <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(OncotreePrimaryDisease)

  title <- paste0(sub('\\.', '-', sample.id), ', ', cell.line.name, ' (', gender, ', ', primary.disease, ')\n', 'Ploidy & karyotype: ', karyotype)

  ## Prototype Tangent
  points.pr <- segment.smoothed.CNA.obj.pr$data[[sample.id]] %>%
    as.data.frame() %>%
    setNames('log2rcn') %>%
    mutate(loci=row_number()) %>%
    bind_cols(probes) %>%
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even')) %>%
    mutate(chr.class.col=case_when(chr.class=='odd' ~ '#0CB702', TRUE ~ '#000000'))

  segments.pr <- segment.smoothed.CNA.obj.pr$output %>%
    filter(ID==sample.id) %>%
    left_join(points.pr %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.pr2 <- points.pr %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.pr %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.pr2 <- segments.pr %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.pr %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.pr <- ggplot(points.pr2, aes(x=loci, y=log2rcn)) +
    ggrastr::rasterize(geom_point(aes(col=chr.class), size=1, show.legend=FALSE), dpi=300, dev='ragg_png') +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.pr2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    # ggrastr::rasterize(geom_segment(data=segments.pr2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5), dpi=300, dev='ragg_png') +
    scale_color_manual(values=c('odd'='limegreen', 'even'='black')) +
    coord_cartesian(ylim=c(-2.5, 2.5)) +
    labs(title='Tangent', x='Genomic position', y=expression(paste({log[2]}, '[Relative copy-number]', sep=''))) +
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
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even')) %>%
    mutate(chr.class.col=case_when(chr.class=='odd' ~ '#0CB702', TRUE ~ '#000000'))

  segments.xy <- segment.smoothed.CNA.obj.xy$output %>%
    filter(ID==sample.id) %>%
    left_join(points.xy %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.xy2 <- points.xy %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.xy %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.xy2 <- segments.xy %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.xy %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.xy <- ggplot(points.xy2, aes(x=loci, y=log2rcn)) +
    ggrastr::rasterize(geom_point(aes(col=chr.class), size=1, show.legend=FALSE), dpi=300, dev='ragg_png') +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.xy2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    # ggrastr::rasterize(geom_segment(data=segments.xy2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5), dpi=300, dev='ragg_png') +
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

  g.title <- ggpubr::as_ggplot(ggpubr::text_grob(title, size=20))
  g.main <- ggpubr::ggarrange(g.pr, g.xy, nrow=1, align='hv')
  g <- ggpubr::ggarrange(g.title, g.main, ncol=1, heights=c(1, 10))

  ggsave(g, file=here('06_CCLE_TangentXY/output/08_make_CBS_plots', paste0(sample.id, '_', gender, '.png')), bg='white', dpi=100, width=10, height=10)
  ggsave(g, file=here('06_CCLE_TangentXY/output/08_make_CBS_plots', paste0(sample.id, '_', gender, '.pdf')), width=10, height=10, useDingbats=TRUE)
}

parallel::mclapply(samples, ggplotCBS, mc.cores=parallel::detectCores() - 2)




## For all samples
samples <- segment.smoothed.CNA.obj.pr$output$ID %>% unique()

ggplotCBS <- function(sample.id) {
  gender <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(Sex)

  print(paste(sample.id, gender))

  i <- grep(sample.id, samples)
  cell.line.name <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(CellLineName)
  primary.disease <- sif %>%
    filter(ModelID==sample.id) %>%
    pull(OncotreePrimaryDisease)

  title <- paste0(sub('\\.', '-', sample.id), ', ', cell.line.name, ' (', gender, ', ', primary.disease, ')')

  ## Prototype Tangent
  points.pr <- segment.smoothed.CNA.obj.pr$data[[sample.id]] %>%
    as.data.frame() %>%
    setNames('log2rcn') %>%
    mutate(loci=row_number()) %>%
    bind_cols(probes) %>%
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even')) %>%
    mutate(chr.class.col=case_when(chr.class=='odd' ~ '#0CB702', TRUE ~ '#000000'))

  segments.pr <- segment.smoothed.CNA.obj.pr$output %>%
    filter(ID==sample.id) %>%
    left_join(points.pr %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.pr2 <- points.pr %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.pr %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.pr2 <- segments.pr %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.pr %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.pr <- ggplot(points.pr2, aes(x=loci, y=log2rcn)) +
    ggrastr::rasterize(geom_point(aes(col=chr.class), size=1, show.legend=FALSE), dpi=300, dev='ragg_png') +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.pr2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    # ggrastr::rasterize(geom_segment(data=segments.pr2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5), dpi=300, dev='ragg_png') +
    scale_color_manual(values=c('odd'='limegreen', 'even'='black')) +
    coord_cartesian(ylim=c(-2.5, 2.5)) +
    labs(title='Tangent', x='Genomic position', y=expression(paste({log[2]}, '[Relative copy-number]', sep=''))) +
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
    mutate(chr.class=case_when(chr %in% paste0('chr', c(seq(1, 21, 2), 'X')) ~ 'odd', TRUE ~ 'even')) %>%
    mutate(chr.class.col=case_when(chr.class=='odd' ~ '#0CB702', TRUE ~ '#000000'))

  segments.xy <- segment.smoothed.CNA.obj.xy$output %>%
    filter(ID==sample.id) %>%
    left_join(points.xy %>% select(chr, start, loci), by=c('chrom'='chr', 'loc.start'='start')) %>%
    rename(start.loci=loci) %>%
    mutate(end.loci=start.loci + num.mark -1)

  points.xy2 <- points.xy %>%
    filter(!chr %in% c('chrX', 'chrY')) %>%
    bind_rows(points.xy %>% filter(chr %in% c('chrX', 'chrY')) %>% mutate(loci=loci+500))

  segments.xy2 <- segments.xy %>%
    filter(!chrom %in% c('chrX', 'chrY')) %>%
    bind_rows(segments.xy %>% filter(chrom %in% c('chrX', 'chrY')) %>% mutate(start.loci=start.loci+500, end.loci=end.loci+500))

  g.xy <- ggplot(points.xy2, aes(x=loci, y=log2rcn)) +
    ggrastr::rasterize(geom_point(aes(col=chr.class), size=1, show.legend=FALSE), dpi=300, dev='ragg_png') +
    geom_hline(yintercept=0, col='gray', linetype='dashed', linewidth=1) +
    geom_hline(yintercept=-1, col='blue', linetype='dashed', linewidth=1) +
    geom_segment(data=segments.xy2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5) +
    # ggrastr::rasterize(geom_segment(data=segments.xy2, aes(x=start.loci, xend=end.loci, y=seg.mean, yend=seg.mean), col='red', linewidth=1.5), dpi=300, dev='ragg_png') +
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

  g.title <- ggpubr::as_ggplot(ggpubr::text_grob(title, size=20))
  g.main <- ggpubr::ggarrange(g.pr, g.xy, nrow=1, align='hv')
  g <- ggpubr::ggarrange(g.title, g.main, ncol=1, heights=c(1, 10))

  ggsave(g, file=here('06_CCLE_TangentXY/output/08_make_CBS_plots/all', paste0(sample.id, '_', gender, '.png')), bg='white', dpi=100, width=10, height=10)
}

parallel::mclapply(samples, ggplotCBS, mc.cores=parallel::detectCores() - 2)
