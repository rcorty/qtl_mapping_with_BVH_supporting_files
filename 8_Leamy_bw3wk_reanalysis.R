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


# baseline fit
dglm_fit <- dglm(formula = bw3wk ~ SEX + father,
                 dformula = ~ SEX + father,
                 data = c$pheno)
summary(dglm_fit)


# scanones
num_perms <- 1e3
num_cores <- 40

father_mat <-  model.matrix(~ 0 + factor(c$pheno$father))[,-1]
sex_and_father_mat <- model.matrix(~ 0 + factor(c$pheno$SEX) + factor(c$pheno$father))[,-1]


so_no_covar <- scanone(cross = c, pheno.col = 'bw3wk')
sop_no_covar <- scanone(cross = c, pheno.col = 'bw3wk',
                        n.perm = num_perms, n.cluster = 40)
evd_no_covar <- evd::fgev(x = sop_no_covar)
so_no_covar$empir.p <- PGEV(q = so_no_covar$lod, gev = evd_no_covar, lower.tail = FALSE)

so_sex_only <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = as.numeric(c$pheno$SEX))
sop_sex_only <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = as.numeric(c$pheno$SEX),
                        n.perm = num_perms, n.cluster = 40)
evd_sex_only <- evd::fgev(x = sop_sex_only)
so_sex_only$empir.p <- PGEV(q = so_sex_only$lod, gev = evd_sex_only, lower.tail = FALSE)

so_father_only <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = father_mat)
sop_father_only <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = father_mat,
                           n.perm = num_perms, n.cluster = 40)
evd_father_only <- evd::fgev(x = sop_father_only)
so_father_only$empir.p <- PGEV(q = so_father_only$lod, gev = evd_father_only, lower.tail = FALSE)

so_sex_and_father <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = sex_and_father_mat)
sop_sex_and_father <- scanone(cross = c, pheno.col = 'bw3wk', addcovar = sex_and_father_mat,
                              n.perm = num_perms, n.cluster = 40)
evd_sex_and_father <- evd::fgev(x = sop_sex_and_father)
so_sex_and_father$empir.p <- PGEV(q = so_sex_and_father$lod, gev = evd_sex_and_father, lower.tail = FALSE)



# # scanonevars
message('Started sov_no_covar at ', Sys.time())
sov_no_covar <- scanonevar(cross = c,
                           mean.formula = bw3wk ~ father + SEX + mean.QTL.add,
                           var.formula = ~ var.QTL.add)

sov_no_covar <- scanonevar.perm(sov = sov_no_covar,
                                random.seed = 27599,
                                n.perms = num_perms,
                                n.cores = num_cores)
message('Finished sov_no_covar at ', Sys.time())

message('Started sov_sex_and_father at ', Sys.time())
sov_sex_and_father <- scanonevar(cross = c,
                                 mean.formula = bw3wk ~ father + SEX + mean.QTL.add,
                                 var.formula = ~ father + SEX + var.QTL.add)

sov_sex_and_father <- scanonevar.perm(sov = sov_sex_and_father,
                                      random.seed = 27599,
                                      n.perms = num_perms,
                                      n.cores = num_cores)
message('Finished sov_sex_and_father at ', Sys.time())


message('Started sov_sex_only at ', Sys.time())
sov_sex_only <- scanonevar(cross = c,
                           mean.formula = bw3wk ~ father + SEX + mean.QTL.add,
                           var.formula = ~ SEX + var.QTL.add)

sov_sex_only <- scanonevar.perm(sov = sov_sex_only,
                                random.seed = 27599,
                                n.perms = num_perms,
                                n.cores = num_cores)
message('Finished sov_sex_only at ', Sys.time())

message('Started sov_father_only at ', Sys.time())
sov_father_only <- scanonevar(cross = c,
                              mean.formula = bw3wk ~ father + SEX + mean.QTL.add,
                              var.formula = ~ father + var.QTL.add)

sov_father_only <- scanonevar.perm(sov = sov_father_only,
                                   random.seed = 27599,
                                   n.perms = num_perms,
                                   n.cores = num_cores)
message('Finished sov_father_only at ', Sys.time())


saveRDS(object = list(so_no_covar = so_no_covar,
                      so_sex_only = so_sex_only,
                      so_father_only = so_father_only,
                      so_sex_and_father = so_sex_and_father,
                      sov_no_covar = sov_no_covar,
                      sov_sex_only = sov_sex_only,
                      sov_father_only = sov_father_only,
                      sov_sex_and_father = sov_sex_and_father),
        file = 'Leamy_bw3wk_analysis.RDS')


results <- readRDS(file = 'Leamy_bw3wk_analysis.RDS')
attach(results)


plot(x = sov_no_covar, y = so_no_covar)
plot(x = sov_sex_only, y = so_sex_only)
plot(x = sov_father_only, y = so_father_only)
plot(x = sov_sex_and_father, y = so_sex_and_father)



plot(x = sov_no_covar, y = so_no_covar, plot.title = 'Bodyweight at Three Weeks, without Father Effects')
ggsave(filename = 'bw3wk_Caos_tests.pdf', height = 2.5, width = 10)

plot(x = sov_sex_and_father, y = so_sex_and_father, plot.title = 'Bodyweight at Three Weeks, accommodating Father Effects')
ggsave(filename = 'bw3wk_DGLM_tests.pdf', height = 2.5, width = 10)


sov_boot <- scanonevar.boot(sov = sov_sex_and_father, n.resamples = 1000, chr = '11', qtl_type = 'mQTL', random.seed = 27599, n.cores = 40)



results$sov$result %>% filter(mean.empir.p < 0.04)

c$pheno$father <- factor(c$pheno$father, labels = 1:9)

# mean_var_plot_model_free(cross = c,
#                          phenotype.name = 'bw3wk',
#                          grouping.factor.names = c('father', 'D11MIT11'))

mean_var_plot_model_based(cross = c,
                          phenotype.name = 'bw3wk',
                          focal.groups = c('father', 'D11MIT11'),
                          nuisance.groups = 'SEX',
                          genotype.names = c('AA', 'AB'),
                          point_size = 4,
                          title = 'Effects of father and D11MIT11\non Bodyweight at Three Weeks')
ggsave(filename = 'bw3wk_mean_var_plot_father_D11MIT11.pdf', height = 5, width = 8)

c$pheno$D11MIT11 <- c$geno$`11`$data[,'D11MIT11']
c$pheno$father <- factor(x = c$pheno$father, labels = 1:9)

dglm_fit <- dglm(formula = bw3wk ~ SEX + father,
                 dformula = ~ SEX + father,
                 data = c$pheno)

anova(lm(formula = dglm_fit$dispersion.fit$y ~ c$pheno$SEX + c$pheno$father),
      lm(formula = dglm_fit$dispersion.fit$y ~ c$pheno$SEX))


lm_fit <- lm(formula = bw3wk ~ SEX + father,
             data = c$pheno)

c$pheno$lm_residual <- lm_fit$residuals

c$pheno$dglm_residual <- dglm_fit$residuals
c$pheno$dglm_weight <- dglm_fit$weights/mean(dglm_fit$weights)


mean(dglm_fit$weights)
plot(x = dglm_fit$residuals, y = dglm_fit$weights)

ggplot(data = c$pheno,
       mapping = aes(x = father, y = lm_residual)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
              trim = FALSE) +
  geom_jitter(mapping = aes(color = dglm_weight), width = 0.2) +
  scale_color_continuous(name = 'DGLM weight', low = '#56B1F7', high = '#132B43') +
  scale_y_continuous(name = 'SLM residual') +
  ggtitle(label = 'Residuals from Standard Linear Model\nof Bodyweight at Three Weeks') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = 'bw3wk_residuals_by_father.pdf', width = 8, height = 3)



ggplot(data = c$pheno,
       mapping = aes(x = father, y = bw3wk)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),
              trim = FALSE) +
  geom_jitter(width = 0.2) +
  ggtitle(label = 'Bodyweight at Three Weeks, by father') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = 'bw3wk_by_father.pdf', width = 8, height = 4)




sov2 <- scanonevar(cross = c,
                   mean.formula = bw3wk ~ father + SEX + D11MIT11 + mean.QTL.add,
                   var.formula = ~ father + SEX + D11MIT11 + var.QTL.add)

plot(x = sov2)


# drop highest variance father, just to try it
c_drop_father1 <- subset(x = c, ind = c$pheno$father != '9')
so_drop_father1 <- scanone(cross = c_drop_father1, pheno.col = 'bw3wk', addcovar = model.matrix(~ factor(c_drop_father1$pheno$SEX) + factor(c_drop_father1$pheno$father))[,-1])
sop_drop_father1 <- scanone(cross = c_drop_father1, pheno.col = 'bw3wk', addcovar = model.matrix(~ factor(c_drop_father1$pheno$SEX) + factor(c_drop_father1$pheno$father))[,-1],
                            n.perm = num_perms, n.cluster = 40)
evd_drop_father1 <- evd::fgev(x = sop_drop_father1)
so_drop_father1$empir.p <- PGEV(q = so_drop_father1$lod, gev = evd_drop_father1, lower.tail = FALSE)
