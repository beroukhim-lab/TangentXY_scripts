library(tidyverse)
library(here)

doc.n.transformed <- readRDS(file=here('06_CCLE_TangentXY/output/01_LinearTransformation', 'CCLE_WES_hg38_N_Transformed.rds'))

N.autox <- doc.n.transformed[!grepl('chrY', rownames(doc.n.transformed)),] %>%
  as.matrix()

## SVD
N.autox.svd <- svd(N.autox)
colnames(N.autox.svd$u) <- colnames(N.autox)
rownames(N.autox.svd$u) <- rownames(N.autox)
rownames(N.autox.svd$v) <- colnames(N.autox)
saveRDS(N.autox.svd, file=here('06_CCLE_TangentXY/output/02_SVD_on_normals', 'N.autox.svd.rds'), compress=FALSE)

d.df <- N.autox.svd$d %>%
  as.data.frame() %>%
  setNames('d') %>%
  mutate(n=1:n())

g <- ggplot(d.df, aes(x=n, y=d)) +
  geom_line() +
  geom_point() +
  labs(y='r (Importance)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('06_CCLE_TangentXY/output/02_SVD_on_normals', 'r.png'), dpi=100, width=16, height=8)

## Construct lower dimensional reference plane
dimension <- 5
N.lower <- N.autox.svd$u[, 1:dimension] %*% diag(N.autox.svd$d[1:dimension]) %*% t(N.autox.svd$v[,1:dimension])
saveRDS(N.lower, file=here('06_CCLE_TangentXY/output/02_SVD_on_normals', paste0('N.lower_23normals_', dimension, 'dimensions.rds')), compress=FALSE)
