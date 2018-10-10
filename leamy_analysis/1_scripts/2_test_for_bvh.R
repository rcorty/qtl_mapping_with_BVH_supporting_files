library(tidyverse)
library(dglm)
library(qtl)
library(vqtl)
library(evd)

PGEV <- function(q, gev, ...) {

  if ((!'estimate' %in% names(gev))) {
    stop("argument 'gev' must have an element named 'estimate' (all gev objects do)")
  }
  if (length(gev$estimate) != 3) {
    stop("gev$estimate must have three elements")
  }

  return(pgev(q = q,
              loc = gev$estimate[1],
              scale = gev$estimate[2],
              shape = gev$estimate[3], ...))
}

MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/yi_2007/Yi2007_M16ixCast_B37_Data.csv'
c <- read.cross(file = url(description = MPD_URL), format = 'csv', genotypes = c('M', 'H'))
c <- calc.genoprob(cross = c, step = 2)

# predictors
c$pheno$SEX %>% table()         # balanced between two sexes
c$pheno$father %>% table()      # 9 fathers

c$pheno$father <- factor(c$pheno$father)
c$pheno$SEX <- factor(c$pheno$SEX)


# phenotypes
phen_names <- names(c$pheno)[7:14]
phen_names


for (phen_name in phen_names) {

  alt <- dglm(formula = formula(paste(phen_name,  '~ SEX + father')),
              dformula = ~ SEX + father,
              data = c$pheno)

  null <-dglm(formula = formula(paste(phen_name,  '~ SEX + father')),
               dformula = ~ SEX,
               data = c$pheno)

  message(phen_name)
  message(pchisq(q = null$m2loglik - alt$m2loglik, df = length(unique(c$pheno$father)) - 1, lower.tail = FALSE))

}




### dig more into bw3wk effects since that's the one we're focusing on

# baseline fits with SLM
# test sex
alt <- lm(formula = bw3wk ~ SEX + father, data = c$pheno)
null <- lm(formula = bw3wk ~ father, data = c$pheno)
pchisq(q = 2*(logLik(alt) - logLik(null)), df = length(unique(c$pheno$SEX)) - 1, lower.tail = FALSE)

# test father
alt <- lm(formula = bw3wk ~ SEX + father, data = c$pheno)
null <- lm(formula = bw3wk ~ SEX, data = c$pheno)
pchisq(q = 2*(logLik(alt) - logLik(null)), df = length(unique(c$pheno$father)) - 1, lower.tail = FALSE)



# baseline fits with DGLM
# test sex
# mQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ father,
             dformula = ~ SEX + father,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = length(unique(c$pheno$SEX)) - 1, lower.tail = FALSE)


# vQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ SEX + father,
             dformula = ~ father,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = length(unique(c$pheno$SEX)) - 1, lower.tail = FALSE)


# mvQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ father,
             dformula = ~ father,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = 2*length(unique(c$pheno$SEX)) - 2, lower.tail = FALSE)


# test father
# mQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ SEX,
             dformula = ~ SEX + father,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = length(unique(c$pheno$father)) - 1, lower.tail = FALSE)


# vQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ SEX + father,
             dformula = ~ SEX ,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = length(unique(c$pheno$father)) - 1, lower.tail = FALSE)


# mvQTL test
alt <- dglm(formula = bw3wk  ~ SEX + father,
            dformula = ~ SEX + father,
            data = c$pheno)

null <- dglm(formula = bw3wk ~ SEX,
             dformula = ~ SEX,
             data = c$pheno)

pchisq(q = null$m2loglik - alt$m2loglik, df = 2*length(unique(c$pheno$father)) - 2, lower.tail = FALSE)


