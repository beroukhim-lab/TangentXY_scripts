library(tidyverse)
library(here)

doc.n.transformed <- readRDS(file=here('03_TCGA_TangentXY/output/01_Linear_transformation_on_normals', 'TCGA_WES_hg19_N_Transformed.rds'))

N.autox <- doc.n.transformed[!grepl('Y', rownames(doc.n.transformed)),] %>%
  as.matrix()

## SVD
N.autox.svd <- svd(N.autox)
colnames(N.autox.svd$u) <- colnames(N.autox)
rownames(N.autox.svd$u) <- rownames(N.autox)
rownames(N.autox.svd$v) <- colnames(N.autox)
saveRDS(N.autox.svd, file=here('03_TCGA_TangentXY/output/02_SVD_on_normals', 'N.autox.svd.rds'), compress=FALSE)

d.df <- N.autox.svd$d %>%
  as.data.frame() %>%
  setNames('d') %>%
  mutate(n=1:n())

g <- ggplot(d.df, aes(x=n, y=d)) +
  geom_line() +
  labs(y='r (Importance)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/02_SVD_on_normals', 'r.png'), dpi=100, width=8, height=8)

top <- 100
g <- ggplot(d.df %>% filter(n < top), aes(x=n, y=d)) +
  geom_line() +
  labs(title=paste0('Top ', top, ' important dimensions'), y='r (Importance)') +
  scale_x_continuous(breaks=seq(0, top, by=10)) +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('03_TCGA_TangentXY/output/02_SVD_on_normals', paste0('r_top', top, '.png')), dpi=100, width=10, height=6)

## Construct lower dimensional reference plane
dimensions <- c(10, 30, 50, 100, 200, 500, 5000, 10441)

for (i in 1:length(dimensions)) {
  dim.i <- dimensions[i]
  print(paste(i, dim.i))

  N.lower.i <- N.autox.svd$u[, 1:dim.i] %*% diag(N.autox.svd$d[1:dim.i]) %*% t(N.autox.svd$v[,1:dim.i])
  saveRDS(N.lower.i, file=here('03_TCGA_TangentXY/output/02_SVD_on_normals', paste0('N.lower_', dim.i, 'dimensions.rds')), compress=FALSE)
}
