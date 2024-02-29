library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

absolute.file <- here('07_SCNAs_in_chrX_and_chrY/data', 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt')
absolute <- read.delim(absolute.file) %>%
  rename(barcode=sample) %>%
  separate(col=array, into=c('project', 'tss', 'participant', 'sample'), sep='-') %>%
  unite(col=TCGA.ID, c('project', 'tss', 'participant'), sep='.') %>%
  mutate(type=case_when(sample=='01' ~ 'TP',
                        sample=='02' ~ 'TR',
                        sample=='03' ~ 'TB',
                        sample=='05' ~ 'TAP',
                        sample=='06' ~ 'TM')) %>%
  left_join(sif, by=c('TCGA.ID', 'type'))

ase.chrx.annot.major.rate.female <- readRDS(file=here('11_chrX_allele_specific_expression/output/04_TCGA_ASE_analysis', 'ase.chrx.annot.major.rate.female.rds'))

ase.chrx.annot.major.rate.female.purity <- ase.chrx.annot.major.rate.female %>%
  left_join(absolute %>% select(SampleID, purity), by='SampleID')

g <- ggplot(ase.chrx.annot.major.rate.female.purity, aes(x=purity, y=major.rate.median)) +
  geom_point(aes(col=ploidy.class)) +
  geom_smooth(aes(col=ploidy.class), method='lm') +
  lims(x=c(0, 1), y=c(0.5, 1)) +
  coord_fixed() +
  lemon::facet_rep_wrap(~karyo.class, nrow=1, repeat.tick.labels=TRUE) +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'FigS7.png'), dpi=100, width=10, height=5)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'FigS7.pdf'), width=10, height=5, useDingbats=TRUE)

reg.wt <- lm(major.rate.mean ~ purity, data=ase.chrx.annot.major.rate.summary3 %>% filter(karyo.class=='WT'))
reg.amp <- lm(major.rate.mean ~ purity, data=ase.chrx.annot.major.rate.summary3 %>% filter(karyo.class=='Whole.Amp'))
summary(reg.wt)
summary(reg.amp)

reg <- ase.chrx.annot.major.rate.summary3 %>%
  group_by(karyo.class, ploidy.class) %>%
  nest() %>%
  mutate(model=map(data, ~lm(major.rate.mean ~ purity, data=.)),
          intercept=map_dbl(model, ~signif(summary(.)$coef[1,1])),
          slope=map_dbl(model, ~signif(summary(.)$coef[2,1]))) %>%
  select(karyo.class, ploidy.class, intercept, slope)

ase.chrx.annot.major.rate.summary4 <- ase.chrx.annot.major.rate.summary3 %>%
  left_join(reg, by=c('karyo.class', 'ploidy.class')) %>%
  # mutate(intercept=0.5, slope=0.5) %>%
  mutate(predicted.ase.value=intercept + slope * purity) %>%
  mutate(diff=major.rate.mean - predicted.ase.value) %>%
  filter(karyo.class=='WT') %>%
  group_by(project, ploidy.class) %>%
  mutate(diff.median=median(diff), n=n()) %>%
  ungroup() %>%
  arrange(ploidy.class, diff.median) %>%
  mutate(project=factor(.$project, levels=.$project %>% unique()))

g <- ggplot(ase.chrx.annot.major.rate.summary4, aes(x=project, y=diff)) +
  # geom_violin() +
  geom_boxplot(aes(fill=ploidy.class), outlier.shape=NA) +
  # facet_grid(ploidy.class ~ karyo.class) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('output/FiguresForPaper', 'Purity_corrected_ASE_value.png'), dpi=100, width=18, height=8)

