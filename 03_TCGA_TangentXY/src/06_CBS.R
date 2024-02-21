library(tidyverse)
library(here)

library(DNAcopy)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))
probes <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'probes.rds')) %>%
  mutate(chr=factor(.$chr, levels=.$chr %>% unique()))

## Run CBS on each tumor type
tumor.types <- sif$project %>% unique()
Tn.files <- list.files(here('03_TCGA_TangentXY/output/03_TangentXY/byTumorType'), full.names=TRUE)

for (i in seq_along(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  tumor.type.i.Tn.file <- grep(paste0('Tn_', tumor.type.i), Tn.files, value=TRUE)
  print(paste(i, tumor.type.i, tumor.type.i.Tn.file))

  Tn.i <- readRDS(tumor.type.i.Tn.file)

  sample.names.i <- colnames(Tn.i)

  CNA.obj <- CNA(Tn.i, probes$chr, probes$start, data.type='logratio', sampleid=sample.names.i)

  smoothed.CNA.obj <- smooth.CNA(CNA.obj)

  segment.smoothed.CNA.obj <- segment(smoothed.CNA.obj, verbose=1)

  saveRDS(segment.smoothed.CNA.obj, file=here('03_TCGA_TangentXY/output/06_CBS', paste0('CBS_', tumor.type.i, '.rds')), compress=FALSE)
}

## Merge all tumor types
for (i in seq_along(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))

  segment.smoothed.CNA.obj.i <- readRDS(file=here('03_TCGA_TangentXY/output/06_CBS', paste0('CBS_', tumor.type.i, '.rds')))
  seg.i <- segment.smoothed.CNA.obj.i$output

  if (i==1) {
    segment <- seg.i
  } else {
    segment <- segment %>% bind_rows(seg.i)
  }
}
saveRDS(segment, file=here('03_TCGA_TangentXY/output/06_CBS', 'CBS_all.rds'), compress=FALSE)


## Make plots
for (i in seq_along(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))
  
  out.dir.whole <- here('03_TCGA_TangentXY/output/06_CBS/plot/', tumor.type.i, '/whole')
  out.dir.xy <- here('03_TCGA_TangentXY/output/06_CBS/plot/', tumor.type.i, '/XY')

  if (!dir.exists(out.dir.whole)) {
    dir.create(out.dir.whole, recursive=TRUE)
  }
  if (!dir.exists(out.dir.xy)) {
    dir.create(out.dir.xy, recursive=TRUE)
  }

  segment.smoothed.CNA.obj.i <- readRDS(file=here('03_TCGA_TangentXY/output/06_CBS', paste0('CBS_', tumor.type.i, '.rds')))

  sample.names.i <- segment.smoothed.CNA.obj.i$output$ID %>% unique()

  for (j in seq_along(sample.names.i)) {
    sample.j <- sample.names.i[j]
    gender.j <- sif %>%
      filter(SampleID==sample.j) %>%
      pull(Gender)
    print(paste(j, sample.j, gender.j))

    png(file=here(out.dir.whole, paste0(j, '_TangentXY_CBS_', sample.j, '_', gender.j, '.png')))
    plotSample(segment.smoothed.CNA.obj.i, sampleid=sample.j, plot.type="w", ylim=c(-2.5, 2.5))
    dev.off()

    png(file=here(out.dir.xy, paste0(j, '_TangentXY_CBS_', sample.j, '_', gender.j, '.png')))
    plotSample(segment.smoothed.CNA.obj.i, sampleid=sample.j, plot.type="w", chromlist=c('X', 'Y'), ylim=c(-2.5, 2.5))
    dev.off()
  }
}
