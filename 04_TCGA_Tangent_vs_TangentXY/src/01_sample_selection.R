library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.n <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_N_QCed_commonCNVremoved.rds'))
doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

batch.table.female <- sif %>%
  filter(!project %in% c('BRCA', 'COAD', 'KIRC', 'OV', 'READ', 'UCEC')) %>%
  filter(SampleID %in% c(colnames(doc.n), colnames(doc.t))) %>%
  filter(Gender=='Female', type %in% c('NB', 'TP')) %>%
  select(batch, type) %>%
  table() %>%
  as.data.frame.matrix() %>%
  mutate(total=NB+TP) %>%
  arrange(desc(total)) %>%
  rownames_to_column('batch')

batch.table.male <- sif %>%
  filter(!project %in% c('BRCA', 'COAD', 'KIRC', 'OV', 'READ', 'UCEC')) %>%
  filter(SampleID %in% c(colnames(doc.n), colnames(doc.t))) %>%
  filter(Gender=='Male', type %in% c('NB', 'TP')) %>%
  select(batch, type) %>%
  table() %>%
  as.data.frame.matrix() %>%
  mutate(total=NB+TP) %>%
  arrange(desc(total)) %>%
  rownames_to_column('batch')

## Manually pick up batches to be analyzed
## 1. Analyses done at the Broad Institute
## 2. Batches with similar number of NB and TP

## Female
batches.to.be.analyzed.female <- c('BI_A202', 'BI_A14W', 'BI_A13W', 'BI_A22D', 'BI_A10S', 'BI_A17V', 'BI_A25L', 'BI_A21Z', 'BI_A18F', 'BI_1753', 'BI_2393', 'BI_A22Z', 'BI_2036', 'BI_A21A', 'BI_1696', 'BI_1705', 'BI_2394', 'BI_A289', 'BI_2238', 'BI_A23U', 'BI_1486', 'BI_1494', 'BI_2024', 'BI_2253', 'BI_A16O', 'BI_A20C', 'BI_A29Q')
samples.to.be.analyzed.female <- sif %>%
  filter(SampleID %in% c(colnames(doc.n), colnames(doc.t))) %>%
  filter(batch %in% batches.to.be.analyzed.female) %>%
  filter(type %in% c('NB', 'TP')) %>%
  filter(Gender=='Female')

tss.table.female <- samples.to.be.analyzed.female %>%
  select(tss, type) %>%
  table() %>%
  as.data.frame.matrix() %>%
  mutate(total=NB+TP) %>%
  arrange(desc(total)) %>%
  rownames_to_column('tss')

pts.to.be.analyzed.female <- samples.to.be.analyzed.female %>%
  select(TCGA.ID, project, type, Gender, tss, batch) %>%
  mutate(count=1) %>%
  pivot_wider(names_from=type, values_from=count) %>%
  mutate(batch=factor(.$batch, levels=batches.to.be.analyzed.female)) %>%
  mutate(tss=factor(.$tss, levels=tss.table.female$tss)) %>%
  arrange(batch, tss) %>%
  filter(NB==1 & TP==1) %>%
  pull(TCGA.ID) %>%
  head(500)

batch.tss.female.tp <- samples.to.be.analyzed.female %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.female) %>%
  filter(type=='TP') %>%
  select(batch, tss) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(TP=Freq)

batch.tss.female.nb <- samples.to.be.analyzed.female %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.female) %>%
  filter(type=='NB') %>%
  select(batch, tss) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(NB=Freq)

batch.tss.female <- batch.tss.female.tp %>%
  full_join(batch.tss.female.nb, by=c('batch', 'tss')) %>%
  filter(TP!=0) %>%
  arrange(desc(TP)) %>%
  mutate(batch.no=1:n())

normals.order.female <- samples.to.be.analyzed.female %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.female) %>%
  filter(type=='NB') %>%
  left_join(batch.tss.female %>% select(-c('NB', 'TP')), by=c('batch', 'tss')) %>%
  group_by(batch.no) %>%
  mutate(n=1:n()) %>%
  ungroup() %>%
  arrange(n, batch.no)


## Male
batches.to.be.analyzed.male <- c('BI_1494', 'BI_1434', 'BI_1492', 'BI_1912', 'BI_1870', 'BI_2244', 'BI_2394', 'BI_2024', 'BI_A29Q', 'BI_1845', 'BI_A22Z', 'BI_1491', 'BI_1468', 'BI_2253', 'BI_A34A', 'BI_A24P', 'BI_A289', 'BI_A30X', 'BI_1696', 'BI_2078', 'BI_2260', 'BI_A23U')
samples.to.be.analyzed.male <- sif %>%
  filter(SampleID %in% c(colnames(doc.n), colnames(doc.t))) %>%
  filter(batch %in% batches.to.be.analyzed.male) %>%
  filter(type %in% c('NB', 'TP')) %>%
  filter(Gender=='Male')

tss.table.male <- samples.to.be.analyzed.male %>%
  select(tss, type) %>%
  table() %>%
  as.data.frame.matrix() %>%
  mutate(total=NB+TP) %>%
  arrange(desc(total)) %>%
  rownames_to_column('tss')

pts.to.be.analyzed.male <- samples.to.be.analyzed.male %>%
  select(TCGA.ID, project, type, Gender, tss, batch) %>%
  mutate(count=1) %>%
  pivot_wider(names_from=type, values_from=count) %>%
  mutate(batch=factor(.$batch, levels=batches.to.be.analyzed.male)) %>%
  mutate(tss=factor(.$tss, levels=tss.table.male$tss)) %>%
  arrange(batch, tss) %>%
  filter(NB==1 & TP==1) %>%
  pull(TCGA.ID) %>%
  head(500)

batch.tss.male.tp <- samples.to.be.analyzed.male %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.male) %>%
  filter(type=='TP') %>%
  select(batch, tss) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(TP=Freq)

batch.tss.male.nb <- samples.to.be.analyzed.male %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.male) %>%
  filter(type=='NB') %>%
  select(batch, tss) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>%
  rename(NB=Freq)

batch.tss.male <- batch.tss.male.tp %>%
  full_join(batch.tss.male.nb, by=c('batch', 'tss')) %>%
  filter(TP!=0) %>%
  arrange(desc(TP)) %>%
  mutate(batch.no=1:n())

normals.order.male <- samples.to.be.analyzed.male %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.male) %>%
  filter(type=='NB') %>%
  left_join(batch.tss.male %>% select(-c('NB', 'TP')), by=c('batch', 'tss')) %>%
  group_by(batch.no) %>%
  mutate(n=1:n()) %>%
  ungroup() %>%
  arrange(n, batch.no)


normals.female <- sif %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.female & type=='NB') %>%
  mutate(TCGA.ID=factor(.$TCGA.ID, levels=normals.order.female$TCGA.ID)) %>%
  arrange(TCGA.ID) %>%
  pull(SampleID)

normals.male <- sif %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.male & type=='NB') %>%
  mutate(TCGA.ID=factor(.$TCGA.ID, levels=normals.order.male$TCGA.ID)) %>%
  arrange(TCGA.ID) %>%
  pull(SampleID)

tumors.female <- sif %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.female & type=='TP') %>%
  pull(SampleID)

tumors.male <- sif %>%
  filter(TCGA.ID %in% pts.to.be.analyzed.male & type=='TP') %>%
  pull(SampleID)

save(normals.female, normals.male, tumors.female, tumors.male, file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

