library(tidyverse)
library(here)

library(DNAcopy)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(DepMap_ID=gsub('-', '.', DepMap_ID))

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

CNA.obj.xy <- CNA(Tn, probes$chr, probes$start, data.type='logratio', sampleid=sample.names)

smoothed.CNA.obj.xy <- smooth.CNA(CNA.obj.xy)

segment.smoothed.CNA.obj.xy <- segment(smoothed.CNA.obj.xy, verbose=1)

saveRDS(segment.smoothed.CNA.obj.xy, file=here('06_CCLE_TangentXY/output/07_CBS', 'CBS_TangentXY.rds'), compress=FALSE)


## Make plots
out.dir.whole <- here('06_CCLE_TangentXY/output/07_CBS/plot/TangentXY/whole')
out.dir.xy <- here('06_CCLE_TangentXY/output/07_CBS/plot/TangentXY/XY')

if (!dir.exists(out.dir.whole)) {
  dir.create(out.dir.whole, recursive=TRUE)
}
if (!dir.exists(out.dir.xy)) {
  dir.create(out.dir.xy, recursive=TRUE)
}

for (i in 1:length(sample.names)) {
  sample.i <- sample.names[i]
  gender.i <- sif %>%
    filter(DepMap_ID==sample.i) %>%
    pull(sex)
  print(paste(i, sample.i, gender.i))

  png(file=here(out.dir.whole, paste0(i, '_TangentXY_CBS_', sample.i, '_', gender.i, '.png')))
  plotSample(segment.smoothed.CNA.obj.xy, sampleid=sample.i, plot.type="w", ylim=c(-2.5, 2.5))
  dev.off()

  png(file=here(out.dir.xy, paste0(i, '_TangentXY_CBS_', sample.i, '_', gender.i, '.png')))
  plotSample(segment.smoothed.CNA.obj.xy, sampleid=sample.i, plot.type="w", chromlist=c('chrX', 'chrY'), ylim=c(-2.5, 2.5))
  dev.off()
}


## Make figures for publication
male1 <- 'ACH.000130'
female1 <- 'ACH.000462'

samples <- c(male1, female1)

for (i in 1:length(samples)) {
  sample.i <- samples[i]
  gender.i <- sif %>%
    filter(DepMap_ID==sample.i) %>%
    pull(sex)
  print(paste(i, sample.i, gender.i))

  pdf(file=here('06_CCLE_TangentXY/output/07_CBS', paste0(i, '_', sample.i, '_', gender.i, '_TangentXY_whole.pdf')))
  plotSample(segment.smoothed.CNA.obj.xy, sampleid=sample.i, plot.type="w", ylim=c(-2.5, 2.5))
  dev.off()

  pdf(file=here('06_CCLE_TangentXY/output/07_CBS', paste0(i, '_', sample.i, '_', gender.i, '_TangentXY_XY.pdf')))
  plotSample(segment.smoothed.CNA.obj.xy, sampleid=sample.i, plot.type="w", chromlist=c('chrX', 'chrY'), ylim=c(-2.5, 2.5))
  dev.off()
}
