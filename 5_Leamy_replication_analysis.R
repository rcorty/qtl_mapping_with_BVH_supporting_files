library(qtl)
library(dplyr)

MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/yi_2007/Yi2007_M16ixCast_B37_Data.csv'
c <- read.cross(file = url(description = MPD_URL), format = 'csv', genotypes = c('M', 'H'))
c <- calc.genoprob(cross = c, step = 1)


c$pheno$father <- factor(c$pheno$father)
c$pheno$SEX <- factor(c$pheno$SEX)
c$pheno$subfam <- factor(c$pheno$subfam)

c$pheno <- c$pheno %>% group_by(subfam) %>% summarise(litter_size = n()) %>% left_join(c$pheno)

standardize <- lm(formula = bw3wk ~ SEX + subfam, data = c$pheno)
c$pheno$bw3wk_std <- residuals(object = standardize)

so <- scanone(cross = c, pheno.col = 'bw3wk_std')
sop <- scanone(cross = c, pheno.col = 'bw3wk_std', n.perm = 1000, n.cluster = 40)

plot(x = so, ylim = c(0, 5))
add.threshold(out = so, perms = sop, alpha = 0.05)
