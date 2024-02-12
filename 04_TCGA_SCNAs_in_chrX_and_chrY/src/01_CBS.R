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
for (i in 5:11) {
  tumor.type.i <- tumor.types[i]
  tumor.type.i.Tn.file <- grep(paste0('Tn_', tumor.type.i), Tn.files, value=TRUE)
  print(paste(i, tumor.type.i, tumor.type.i.Tn.file))

  Tn.i <- readRDS(tumor.type.i.Tn.file)

  sample.names.i <- colnames(Tn.i)

  CNA.obj <- CNA(Tn.i, probes$chr, probes$start, data.type='logratio', sampleid=sample.names.i)

  smoothed.CNA.obj <- smooth.CNA(CNA.obj)

  segment.smoothed.CNA.obj <- segment(smoothed.CNA.obj, verbose=1)

  saveRDS(segment.smoothed.CNA.obj, file=here('04_TCGA_SCNAs_in_chrX_and_chrY/output/01_CBS', paste0('CBS_', tumor.type.i, '.rds')), compress=FALSE)
}

## Merge all tumor types
for (i in 1:length(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))

  if (rescaling.after.tangent) {
    segment.smoothed.CNA.obj.y.i <- readRDS(file=here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  } else {
    segment.smoothed.CNA.obj.y.i <- readRDS(file=here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  }
  seg.i <- segment.smoothed.CNA.obj.y.i$output

  if (i==1) {
    segment <- seg.i
  } else {
    segment <- segment %>% bind_rows(seg.i)
  }
}
saveRDS(segment, file=here(paste0('output/15_CBS/rescalingAfterTangent', '/dimension', dim), 'CBS_all.RData'), compress=FALSE)

## Extract segmented data and save as text file for downstream GISTIC
for (i in 1:length(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))

  if (rescaling.after.tangent) {
    segment.smoothed.CNA.obj.y <- readRDS(file=here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  } else {
    segment.smoothed.CNA.obj.y <- readRDS(file=here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  }
  seg <- segment.smoothed.CNA.obj.y$output %>%
    rename(SampleID=ID)

  female.samples <- sif %>%
    filter(SampleID %in% unique(seg$SampleID) & gender=='female') %>%
    pull(SampleID)

  male.samples <- sif %>%
    filter(SampleID %in% unique(seg$SampleID) & gender=='male') %>%
    pull(SampleID)

  ## Adjust chrX and chrY in male so that the values are relative to monosomy
  seg.adj <- seg %>%
    filter(SampleID %in% c(female.samples, male.samples)) %>%
    mutate(seg.mean=case_when(SampleID %in% male.samples & chrom %in% c('X', 'Y') ~ seg.mean + 1, TRUE ~ seg.mean)) %>%
    filter(!(SampleID %in% female.samples & chrom=='Y'))

  ## Pick up one sample per each patient, prioritise TP over TM and TAP
  ## For SKCM, pick up only TM
  if (tumor.type.i=='SKCM') {
    unique.samples <- sif %>%
      filter(type=='TM') %>%
      filter(SampleID %in% unique(seg.adj$SampleID)) %>%
      distinct(ID, .keep_all=TRUE) %>%
      pull(SampleID)
  } else {
    unique.samples <- sif %>%
      filter(NT=='Tumor') %>%
      mutate(type=factor(.$type, levels=c('TP', 'TM', 'TAP', 'TB'))) %>%
      filter(SampleID %in% unique(seg.adj$SampleID)) %>%
      arrange(type) %>%
      distinct(ID, .keep_all=TRUE) %>%
      pull(SampleID)
  }

  if (i==1) {
    samples.used.for.gistic <- unique.samples
  } else {
    samples.used.for.gistic <- samples.used.for.gistic %>% append(unique.samples)
  }

  seg.unique.samples <- seg.adj %>%
    filter(SampleID %in% unique.samples)

  if (rescaling.after.tangent) {
    write.table(seg.unique.samples, file=here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim, '/seg_for_gistic'), paste0('Seg_ZSTangent_10460normals_50dimensions_', tumor.type.i, '.txt')), sep='\t', quote=FALSE, row.names=FALSE)
  } else {
    write.table(seg.unique.samples, file=here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim, '/seg_for_gistic'), paste0('Seg_ZSTangent_10460normals_50dimensions_', tumor.type.i, '.txt')), sep='\t', quote=FALSE, row.names=FALSE)
  }
}

## Show breakdown of sample type for each tumor type
sif %>%
  filter(SampleID %in% samples.used.for.gistic) %>%
  select(project, type) %>%
  table()

## Make plots
for (i in 1:length(tumor.types)) {
  tumor.type.i <- tumor.types[i]
  print(paste(i, tumor.type.i))
  
  if (rescaling.after.tangent) {
    out.dir.whole <- here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim, '/plot/', tumor.type.i, '/whole'))
    out.dir.xy <- here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim, '/plot/', tumor.type.i, '/XY'))
  } else {
    out.dir.whole <- here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim, '/plot/', tumor.type.i, '/whole'))
    out.dir.xy <- here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim, '/plot/', tumor.type.i, '/XY'))
  }

  if (!dir.exists(out.dir.whole)) {
    dir.create(out.dir.whole, recursive=TRUE)
  }
  if (!dir.exists(out.dir.xy)) {
    dir.create(out.dir.xy, recursive=TRUE)
  }

  if (rescaling.after.tangent) {
    segment.smoothed.CNA.obj.y <- readRDS(file=here(paste0('output/15_CBS/rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  } else {
    segment.smoothed.CNA.obj.y <- readRDS(file=here(paste0('output/15_CBS/no_rescalingAfterTangent/dimension', dim), paste0('CBS_', tumor.type.i, '.RData')))
  }
  sample.names.i <- segment.smoothed.CNA.obj.y$output$ID %>% unique()

  for (j in 1:length(sample.names.i)) {
    sample.j <- sample.names.i[j]
    gender.j <- sif %>%
      filter(SampleID==sample.j) %>%
      pull(gender)
    print(paste(j, sample.j, gender.j))

    png(file=here(out.dir.whole, paste0(j, '_ZSTangent_CBS_', sample.j, '_', gender.j, '.png')))
    plotSample(segment.smoothed.CNA.obj.y, sampleid=sample.j, plot.type="w", ylim=c(-2.5, 2.5))
    dev.off()

    png(file=here(out.dir.xy, paste0(j, '_ZSTangent_CBS_', sample.j, '_', gender.j, '.png')))
    plotSample(segment.smoothed.CNA.obj.y, sampleid=sample.j, plot.type="w", chromlist=c('X', 'Y'), ylim=c(-2.5, 2.5))
    dev.off()
  }
}

