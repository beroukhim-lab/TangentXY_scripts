library(tidyverse)
library(here)

sample.annot <- read.delim(here('01_TCGA_seg_arm_level_average/data', 'clinical.tsv'))

#Segmentation files downloaded from the NIH GDC Data Portal based on a manifest file (../data/gdc_manifest_SNP6_seg_grch38.txt)
seg.files <- list.files(here('01_TCGA_seg_arm_level_average/data/seg'), , pattern='v2\\.txt$', recursive=TRUE, full.names=TRUE)
seg.list <- lapply(seg.files, read_delim)

seg.df <- seg.list %>%
  bind_rows() %>%
  as.data.frame()
# saveRDS(seg.df, here('01_TCGA_seg_arm_level_average/tmp', 'seg.df.rds'), compress=FALSE)
# seg.df <- readRDS(file=here('01_TCGA_seg_arm_level_average/tmp', 'seg.df.rds'))

aliquots <- seg.df$GDC_Aliquot %>% unique()
submitter.id <- TCGAutils::UUIDtoBarcode(aliquots, from_type='aliquot_ids')
biospec.list <- parallel::mclapply(submitter.id$portions.analytes.aliquots.submitter_id, TCGAutils::TCGAbiospec, mc.cores=parallel::detectCores()-2)
biospec <- biospec.list %>%
  bind_rows() %>%
  select(c('submitter_id', 'sample_definition', 'sample', 'vial', 'portion', 'analyte', 'plate', 'center'))
annotation <- submitter.id %>%
  bind_cols(biospec) %>%
  mutate(type=case_when(sample=='01' ~ 'TP',
                        sample=='02' ~ 'TR',
                        sample=='03' ~ 'TB',
                        sample=='05' ~ 'TAP',
                        sample=='06' ~ 'TM',
                        sample=='10' ~ 'NB',
                        sample=='11' ~ 'NT',
                        sample=='12' ~ 'NBC',
                        sample=='14' ~ 'NBM')) %>%
  left_join(sample.annot %>% select(case_submitter_id, gender) %>% distinct(), by=c('submitter_id'='case_submitter_id'))

tumor.aliquot <- annotation %>%
  filter(type %in% c('TP', 'TR', 'TB', 'TAP', 'TM')) %>%
  select(portions.analytes.aliquots.aliquot_id, gender)

hg38.annot <- rCGH::hg38 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom))) %>%
  mutate(centromerCenter=centromerStart+1500000)

seg.tumor <- seg.df %>%
  inner_join(tumor.aliquot, by=c('GDC_Aliquot'='portions.analytes.aliquots.aliquot_id')) %>%
  filter(gender %in% c('female', 'male')) %>%
  left_join(hg38.annot, by=c('Chromosome'='chrom')) %>%
  mutate(Chromosome=factor(.$Chromosome, levels=c(1:22, 'X', 'Y'))) %>%
  mutate(arm=case_when(End < centromerCenter ~ 'p',
                        Start < centromerCenter & End > centromerCenter ~ 'pcq',
                        Start > centromerCenter ~ 'q'))

seg.tumor.newEnd <- seg.tumor %>%
  filter(arm=='pcq') %>%
  mutate(End=centromerCenter) %>%
  mutate(arm='p')

seg.tumor.newStart <- seg.tumor %>%
  filter(arm=='pcq') %>%
  mutate(Start=centromerCenter) %>%
  mutate(arm='q')

seg.tumor2 <- seg.tumor %>%
  filter(arm!='pcq') %>%
  bind_rows(seg.tumor.newEnd) %>%
  bind_rows(seg.tumor.newStart) %>%
  arrange(GDC_Aliquot, Chromosome, Start)

seg.tumor.summary <- seg.tumor2 %>%
  mutate(span=End - Start) %>%
  mutate(segment.span=Segment_Mean * span) %>%
  group_by(GDC_Aliquot, Chromosome, arm, gender) %>%
  summarize(span.total=sum(span), segment.span.total=sum(segment.span)) %>%
  mutate(seg.arm.mean=segment.span.total/span.total) %>%
  ungroup() %>%
  mutate(auto.sex=case_when(Chromosome=='X' ~ 'X', Chromosome=='Y' ~ 'Y', TRUE ~ 'Autosomes'))

g <- ggplot(seg.tumor.summary, aes(x=auto.sex, y=seg.arm.mean)) +
  geom_hline(yintercept=0, col='gray', linetype='dashed') +
  geom_hline(yintercept=log2(0.5), col='blue', linetype='dashed') +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-5, 3, by=1)) +
  facet_wrap(~gender, nrow=1, labeller=as_labeller(c('female'='Female', 'male'='Male'))) +
  labs(y='log2(Relative CN)', title='Arm level average log2(Relative CN)') +
  theme_bw(base_size=30) +
  theme(panel.grid.minor.y=element_blank()) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('output', 'Arm_level_log2RCN.png'), dpi=100, width=16, height=8)
ggsave(g, file=here('output', 'Arm_level_log2RCN.pdf'), width=16, height=8)
