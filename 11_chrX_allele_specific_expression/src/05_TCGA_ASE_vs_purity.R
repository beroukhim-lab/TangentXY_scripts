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
  geom_point(aes(col=ploidy.class), shape=21, alpha=0.5) +
  geom_smooth(aes(col=ploidy.class), method='lm') +
  lims(y=c(0.5, 1)) +
  scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  coord_fixed(ratio=1.5) +
  lemon::facet_rep_wrap(~karyo.class, nrow=1, repeat.tick.labels=TRUE) +
  labs(x='Tumor purity', y='ASE value', col='Ploidy') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'FigS7.png'), dpi=100, width=10, height=5)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'FigS7.pdf'), width=10, height=5, useDingbats=TRUE)

reg <- ase.chrx.annot.major.rate.female.purity %>%
  group_by(karyo.class, ploidy.class) %>%
  nest() %>%
  mutate(model=map(data, ~lm(major.rate.mean ~ purity, data=.)),
          intercept=map_dbl(model, ~signif(summary(.)$coef[1,1])),
          slope=map_dbl(model, ~signif(summary(.)$coef[2,1]))) %>%
  select(karyo.class, ploidy.class, intercept, slope)

# ase.chrx.annot.major.rate.diff <- ase.chrx.annot.major.rate.female.purity %>%
#   left_join(reg, by=c('karyo.class', 'ploidy.class')) %>%
#   # mutate(intercept=0.5, slope=0.5) %>%
#   mutate(expected.ase.value=intercept + slope * purity) %>%
#   mutate(diff=expected.ase.value - major.rate.mean) %>%
#   filter(karyo.class=='No.Arm-level.Alt') %>%
#   group_by(project, ploidy.class) %>%
#   mutate(diff.median=median(diff), n=n(), x=length(diff[diff > 0])) %>%
#   ungroup() %>%
#   arrange(ploidy.class, desc(diff.median)) %>%
#   mutate(project=factor(.$project, levels=.$project %>% unique()))

# one.sample.t.test <- ase.chrx.annot.major.rate.diff %>%
#   filter(n >= 2) %>%
#   group_by(karyo.class, project, ploidy.class) %>%
#   # rstatix::t_test(diff ~ 1, mu=0, alternative='greater') %>%
#   rstatix::t_test(diff ~ 1, mu=0, alternative='two.sided') %>%
#   rstatix::adjust_pvalue(method = "bonferroni") %>%
#   rstatix::add_significance() %>%
#   rstatix::add_xy_position(x='project') %>%
#   mutate(xmin=case_when(ploidy.class=='Diploid' ~ xmin - 0.2, TRUE ~ xmin + 0.2)) %>%
#   mutate(xmax=case_when(ploidy.class=='Diploid' ~ xmax - 0.2, TRUE ~ xmax + 0.2))

# one.sample.wilcox.test <- ase.chrx.annot.major.rate.diff %>%
#   filter(n >= 2) %>%
#   group_by(karyo.class, project, ploidy.class) %>%
#   rstatix::wilcox_test(diff ~ 1, mu=0, alternative = "two.sided") %>%
#   rstatix::adjust_pvalue(method = "bonferroni") %>%
#   rstatix::add_significance() %>%
#   rstatix::add_xy_position(x='project') %>%
#   mutate(xmin=case_when(ploidy.class=='Diploid' ~ xmin - 0.2, TRUE ~ xmin + 0.2)) %>%
#   mutate(xmax=case_when(ploidy.class=='Diploid' ~ xmax - 0.2, TRUE ~ xmax + 0.2))

ase.chrx.annot.major.rate.diff.summary <- ase.chrx.annot.major.rate.female.purity %>%
  left_join(reg, by=c('karyo.class', 'ploidy.class')) %>%
  # mutate(intercept=0.5, slope=0.5) %>%
  mutate(expected.ase.value=intercept + slope * purity) %>%
  mutate(diff=expected.ase.value - major.rate.mean) %>%
  filter(karyo.class=='No.Arm-level.Alt') %>%
  group_by(karyo.class, project, ploidy.class) %>%
  summarize(diff.median=median(diff), N=n(), X=length(diff[diff > 0]), max.diff=max(diff)) %>%
  arrange(ploidy.class, desc(diff.median)) %>%
  ungroup() %>%
  mutate(project=factor(.$project, levels=.$project %>% unique()))

binom.test.df <- mapply(rstatix::binom_test, x=ase.chrx.annot.major.rate.diff.summary$X, n=ase.chrx.annot.major.rate.diff.summary$N, alternative='two.sided') %>%
  t() %>%
  as.data.frame() %>%
  bind_cols(ase.chrx.annot.major.rate.diff.summary %>% select(-N)) %>%
  mutate(y.position=max.diff + 0.02) %>%
  mutate(group1='1', group2='null') %>%
  # mutate(groups=list(c('1', 'null model'))) %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance() %>%
  # rstatix::add_xy_position() %>%
  arrange(project, ploidy.class) %>%
  mutate(x=rep(1:(n()/2), 2) %>% sort) %>%
  mutate(xmin=case_when(ploidy.class=='Diploid' ~ x- 0.2, TRUE ~ x + 0.2)) %>%
  mutate(xmax=case_when(ploidy.class=='Diploid' ~ x - 0.2, TRUE ~ x + 0.2))

g <- ggplot(ase.chrx.annot.major.rate.diff, aes(x=project, y=diff)) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_boxplot(aes(fill=ploidy.class), position=position_dodge(0.8)) +
  # ggpubr::stat_pvalue_manual(one.sample.t.test, label='p.adj.signif', col='ploidy.class', size=5, hide.ns=TRUE, show.legend=FALSE) +
  # ggpubr::stat_pvalue_manual(one.sample.wilcox.test, label='p.adj.signif', col='ploidy.class', size=5, hide.ns=TRUE, show.legend=FALSE) +
  ggpubr::stat_pvalue_manual(binom.test.df, label='p.adj.signif', col='ploidy.class', size=5, hide.ns=TRUE, show.legend=FALSE) +
  scale_fill_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  scale_color_manual(values=c('Diploid'='#4DAF4A', 'Polyploid'='#FF7F00')) +
  scale_y_continuous(breaks=seq(-3, 3, 1)/10, limits=c(-0.27, 0.27)) +
  labs(y='Difference in allelic ratios\nfrom expected', fill='Ploidy') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'Fig4g.png'), dpi=100, width=18, height=6)
ggsave(g, file=here('11_chrX_allele_specific_expression/output/05_TCGA_ASE_vs_purity', 'Fig4g.pdf'), width=18, height=6, useDingbats=TRUE)

female.inactive.cn1to2 <- t.test(lm.df.l %>% filter(Gender=='Female' & type=='Inactive genes' & cn=='1 < CN < 2') %>% pull(slope), mu=0)
