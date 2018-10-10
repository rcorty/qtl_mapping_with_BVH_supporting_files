library(tidyverse)
library(nlme)
library(car)
library(dglm)

set.seed(27599)

n <- 100
t <- 100


# null locus
Caom_LRs <- Caov_LRs <- Caomv_LRs <- rep(NA, t)
DGLMm_LRs <- DGLMv_LRs <- DGLMmv_LRs <- rep(NA, t)
for (i in 1:t) {
  g <- factor(rbinom(n = n, size = 2, prob = 0.5))
  m <- 0
  s <- 1
  y <- rnorm(n = n, mean = m, sd = s)

  Caom_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'mean')
  DGLMm_LRs[i] <- mean_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caov_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'var')
  DGLMv_LRs[i] <- var_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caomv_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'both')
  DGLMmv_LRs[i] <- joint_test(resp = y, locus = g, covar = NULL)

}

pdf(file = 'Cao_and_DGLM_indistinguishable_null.pdf',
    height = 3,
    width = 8)
par(mfrow = c(1, 3))
plot(x = DGLMm_LRs, Caom_LRs, main = 'mQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMv_LRs, Caov_LRs, main = 'vQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMmv_LRs, Caomv_LRs, main = 'mvQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
dev.off()




# mQTL
Caom_LRs <- Caov_LRs <- Caomv_LRs <- rep(NA, t)
DGLMm_LRs <- DGLMv_LRs <- DGLMmv_LRs <- rep(NA, t)
for (i in 1:t) {
  g <- factor(rbinom(n = n, size = 2, prob = 0.5))
  m <- 0.2*as.numeric(g)
  s <- 1
  y <- rnorm(n = n, mean = m, sd = s)

  Caom_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'mean')
  DGLMm_LRs[i] <- mean_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caov_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'var')
  DGLMv_LRs[i] <- var_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caomv_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'both')
  DGLMmv_LRs[i] <- joint_test(resp = y, locus = g, covar = NULL)

}

pdf(file = 'Cao_and_DGLM_indistinguishable_mqtl.pdf',
    height = 3,
    width = 8)
par(mfrow = c(1, 3))
plot(x = DGLMm_LRs, Caom_LRs, main = 'mQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMv_LRs, Caov_LRs, main = 'vQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMmv_LRs, Caomv_LRs, main = 'mvQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
dev.off()





# vQTL
Caom_LRs <- Caov_LRs <- Caomv_LRs <- rep(NA, t)
DGLMm_LRs <- DGLMv_LRs <- DGLMmv_LRs <- rep(NA, t)
for (i in 1:t) {
  g <- factor(rbinom(n = n, size = 2, prob = 0.5))
  m <- 0
  s <- log(2+as.numeric(g))
  y <- rnorm(n = n, mean = m, sd = s)

  Caom_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'mean')
  DGLMm_LRs[i] <- mean_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caov_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'var')
  DGLMv_LRs[i] <- var_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caomv_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'both')
  DGLMmv_LRs[i] <- joint_test(resp = y, locus = g, covar = NULL)

}

pdf(file = 'Cao_and_DGLM_indistinguishable_vqtl.pdf',
    height = 3,
    width = 8)
par(mfrow = c(1, 3))
plot(x = DGLMm_LRs, Caom_LRs, main = 'mQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMv_LRs, Caov_LRs, main = 'vQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMmv_LRs, Caomv_LRs, main = 'mvQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
dev.off()



# mvQTL
Caom_LRs <- Caov_LRs <- Caomv_LRs <- rep(NA, t)
DGLMm_LRs <- DGLMv_LRs <- DGLMmv_LRs <- rep(NA, t)
for (i in 1:t) {
  g <- factor(rbinom(n = n, size = 2, prob = 0.5))
  m <- 0.1*as.numeric(g)
  s <- log(2+0.5*as.numeric(g))
  y <- rnorm(n = n, mean = m, sd = s)

  Caom_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'mean')
  DGLMm_LRs[i] <- mean_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caov_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'var')
  DGLMv_LRs[i] <- var_test(resp = y, mean_locus = g, var_locus = g, covar = NULL)

  Caomv_LRs[i] <- LRT_mean_var(gt = g, depvar = y, test = 'both')
  DGLMmv_LRs[i] <- joint_test(resp = y, locus = g, covar = NULL)

}

pdf(file = 'Cao_and_DGLM_indistinguishable_mvqtl.pdf',
    height = 3,
    width = 8)
par(mfrow = c(1, 3))
plot(x = DGLMm_LRs, Caom_LRs, main = 'mQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMv_LRs, Caov_LRs, main = 'vQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
plot(x = DGLMmv_LRs, Caomv_LRs, main = 'mvQTL tests', xlab = 'DGLMm LR', ylab = 'Caom LR')
abline(a = 0, b = 1)
dev.off()
