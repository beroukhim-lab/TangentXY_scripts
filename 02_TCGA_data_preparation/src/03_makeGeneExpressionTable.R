library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt')) %>%
  mutate(sample.id=paste(TCGA.ID, type, sep='.'))

extract_tpm <- function(df) {
  df.extract <- df %>%
    select(c(1:3, 7)) %>%
    setNames(c('gene_id', 'gene_name', 'gene_type', 'tpm')) %>%
    unite(col='gene_info', c('gene_id', 'gene_name', 'gene_type'), sep='|') %>%
    column_to_rownames('gene_info')
}

rnaseq.files <- list.files(here('02_TCGA_data_preparation/data/RNAseq'), recursive=TRUE, pattern='tsv$', full.names=TRUE)

file.ids <- rnaseq.files %>%
  dirname() %>%
  sub('.*/', '', .)

barcodes <- parallel::mclapply(file.ids, TCGAutils::UUIDtoBarcode, from_type='file_id', mc.cores=parallel::detectCores()-2) %>%
  bind_rows() %>%
  rename(barcode=associated_entities.entity_submitter_id)

biospecs <- parallel::mclapply(barcodes$associated_entities.entity_submitter_id, TCGAutils::TCGAbiospec, mc.cores=parallel::detectCores()-2) %>%
  bind_rows() %>%
  mutate(sampleid=paste(gsub('-', '.', submitter_id), sample, sep='.')) %>%
  mutate(type=case_when(sample=='01' ~ 'TP',
                        sample=='02' ~ 'TR',
                        sample=='03' ~ 'TB',
                        sample=='05' ~ 'TAP',
                        sample=='06' ~ 'TM',
                        sample=='07' ~ 'TAM',
                        sample=='11' ~ 'NT')) %>%
  mutate(sample.id=paste(gsub('-', '.', submitter_id), type, sep='.')) %>%
  left_join(sif, by='sample.id')

bb <- barcodes %>%
  bind_cols(biospecs)

rnaseq.df.list <- parallel::mclapply(rnaseq.files, read_delim, skip=6, col_names=FALSE, mc.cores=parallel::detectCores()-2)
tpm.list <- parallel::mclapply(rnaseq.df.list, extract_tpm, mc.cores=parallel::detectCores()-2)

gene.info <- tpm.list[[1]] %>%
  rownames()

egid.info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=sub('\\..*$', '', gene.info), column='ENTREZID', keytype='ENSEMBL') %>%
  as.matrix() %>%
  as.data.frame() %>%
  setNames('egid')

chr.info <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys=sub('\\..*$', '', gene.info), column='CHR', keytype='ENSEMBL') %>%
  as.matrix() %>%
  as.data.frame() %>%
  setNames('chr')

gene.info2 <- paste(gene.info, egid.info$egid, chr.info$chr, sep='|')

tpm.all <- bind_cols(tpm.list) %>%
  setNames(bb$barcode)

rownames(tpm.all) <- gene.info2

## Pick up the samples that have WES data
wes.samples <- bb %>%
  filter(!is.na(SampleID)) %>%
  distinct(SampleID, .keep_all=TRUE)

tpm.selected <- tpm.all %>%
  select(wes.samples$barcode) %>%
  setNames(wes.samples$SampleID)
saveRDS(tpm.selected, here('output', 'tpm.rds'), compress=TRUE)
