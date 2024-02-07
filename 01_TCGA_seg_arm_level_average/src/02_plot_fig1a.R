library(tidyverse)
library(here)

seg.tumor.summary <- readRDS(file=here('01_TCGA_seg_arm_level_average/tmp/01_arm_level_average_calculation', 'seg.tumor.summary.rds'))

g <- ggplot(seg.tumor.summary, aes(x=auto.sex, y=seg.arm.mean)) +
  geom_hline(yintercept=0, col='gray', linetype='dashed') +
  geom_hline(yintercept=log2(0.5), col='blue', linetype='dashed') +
  geom_boxplot() +
  scale_y_continuous(breaks=seq(-5, 3, by=1)) +
  facet_rep_wrap(~gender, nrow=1, labeller=as_labeller(c('female'='Female', 'male'='Male')), repeat.tick.labels=TRUE) +
  labs(y=expression(paste({log[2]}, '[Relative CN]', sep='')), title=expression(paste('Arm level average ', {log[2]}, '[Relative CN]'))) +
  # coord_capped_cart(bottom='both', left='both') +
  theme_classic(base_size=20) +
  theme(strip.background=element_rect(linewidth=0)) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.title.x=element_blank())
ggsave(g, file=here('01_TCGA_seg_arm_level_average/output/02_plot_fig1a', 'Fig1A.png'), dpi=100, width=10, height=5)
ggsave(g, file=here('01_TCGA_seg_arm_level_average/output/02_plot_fig1a', 'Fig1A.pdf'), width=10, height=5)
