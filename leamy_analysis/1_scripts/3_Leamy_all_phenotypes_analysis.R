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

# read in data
MPD_URL <- 'https://phenomedoc.jax.org/QTL_Archive/yi_2007/Yi2007_M16ixCast_B37_Data.csv'
c <- read.cross(file = url(description = MPD_URL), format = 'csv', genotypes = c('M', 'H'))
c <- calc.genoprob(cross = c, step = 5)


c$pheno %>% glimpse()

# predictors
c$pheno$SEX %>% table()         # balanced between two sexes
c$pheno$pgm %>% table()         # all one pgm
c$pheno$father %>% table()      # 9 fathers
c$pheno$mother %>% table()      # too many mothers to estimate their effects well
c$pheno$subfam %>% table()      # too many subfams to estimate their effects well

# phenotypes
phen_names <- names(c$pheno)[c(7, 8, 9, 13, 14)]
phen_names


num_perm <- 1000
num_cores <- 40

# # scanones
# sos <- sops <- list()
# for (phen_name in phen_names) {
#
#   sos[[phen_name]] <- scanone(cross = c,
#                               pheno.col = phen_name,
#                               addcovar = c$pheno$SEX)
#
#   sops[[phen_name]] <- scanone(cross = c,
#                                pheno.col = phen_name,
#                                addcovar = c$pheno$SEX,
#                                n.perm = num_perm,
#                                n.cluster = num_cores)
#
#   evd <- evd::fgev(x = sops[[phen_name]])
#   sos[[phen_name]]$empir.p <- PGEV(q = sos[[phen_name]]$lod, gev = evd, lower.tail = FALSE)
# }
#
# saveRDS(object = sos,
#         file = paste0('Leamy_sos_', num_perm, '_perms.RDS'))
#
#
#
# # vqtl handles factors
# c$pheno$father <- factor(c$pheno$father)
# c$pheno$SEX <- factor(c$pheno$SEX)
#
#
# # scanonevars replicating Cao's test (no BVH)
# sovs_no_bvh <- list()
# for (phen_name in phen_names) {
#
#   message('Starting ', phen_name, ' at ', Sys.time())
#
#   sovs_no_bvh[[phen_name]] <- scanonevar(cross = c,
#                                          mean.formula = formula(paste(phen_name, '~ SEX + father + mean.QTL.add')),
#                                          var.formula = ~ var.QTL.add)
#
#   sovs_no_bvh[[phen_name]] <- scanonevar.perm(sov = sovs_no_bvh[[phen_name]],
#                                               n.perm = num_perm,
#                                               n.cores = num_cores,
#                                               random.seed = 27599)
#
# }
#
# saveRDS(object = sovs_no_bvh,
#         file = paste0('Leamy_sovs_no_BVH_', num_perm, '_perms.RDS'))
#
#
# # scanonevars with BVH
# sovs <- list()
# for (phen_name in phen_names) {
#
#   message('Starting ', phen_name, ' at ', Sys.time())
#
#   sovs[[phen_name]] <- scanonevar(cross = c,
#                                   mean.formula = formula(paste(phen_name, '~ SEX + father + mean.QTL.add')),
#                                   var.formula = ~ SEX + father + var.QTL.add)
#
#   sovs[[phen_name]] <- scanonevar.perm(sov = sovs[[phen_name]],
#                                        n.perm = num_perm,
#                                        n.cores = num_cores,
#                                        random.seed = 27599)
#
# }
#
# saveRDS(object = sovs,
#         file = paste0('Leamy_sovs_', num_perm, '_perms.RDS'))



# plot results
sos <- readRDS(paste0('leamy_analysis/2_results/Leamy_sos_', num_perm, '_perms.RDS'))
caos <- readRDS(paste0('leamy_analysis/2_results/Leamy_sovs_no_BVH_', num_perm, '_perms.RDS'))
sovs <- readRDS(paste0('leamy_analysis/2_results/Leamy_sovs_', num_perm, '_perms.RDS'))


width <- 9
height <- 1.5

plot(x = caos[['bw12d']],     y = sos[['bw12d']], plot.title = 'Bodyweight at 12 days')      # no QTL
ggsave(filename = 'leamy_analysis/3_images/scan_cao_bw12d.pdf', height = height, width = width)

plot(x = caos[['bw3wk']],     y = sos[['bw3wk']], plot.title = 'Bodyweight at 3 weeks')      # no QTL
ggsave(filename = 'leamy_analysis/3_images/scan_cao_bw3wk.pdf', height = height, width = width)

plot(x = caos[['bw6wk']],     y = sos[['bw6wk']], plot.title = 'Bodyweight at 6 weeks')      # concordant, plus one new vQTL
ggsave(filename = 'leamy_analysis/3_images/scan_cao_bw6wk.pdf', height = height, width = width)

plot(x = caos[['scfat12wk']], y = sos[['scfat12wk']], ymax = 4, plot.title = 'subcutaneous fat pad thickness at 12 weeks')  # concordant (three big mQTL/mvQTL)
ggsave(filename = 'leamy_analysis/3_images/scan_cao_scfat12wk.pdf', height = height, width = width)

plot(x = caos[['gfat12wk']],  y = sos[['gfat12wk']], ymax = 4, plot.title = 'gonadal fat pad thickness at 12 weeks')   # concordant (three big mQTL/mvQTL)
ggsave(filename = 'leamy_analysis/3_images/scan_cao_gfat12wk.pdf', height = height, width = width)



plot(x = sovs[['bw12d']],     y = sos[['bw12d']], plot.title = 'Bodyweight at 12 days')      # no QTL
ggsave(filename = 'leamy_analysis/3_images/scan_dglm_bw12d.pdf', height = height, width = width)

plot(x = sovs[['bw3wk']],     y = sos[['bw3wk']], plot.title = 'Bodyweight at 3 weeks')      # new QTL on chr 11
ggsave(filename = 'leamy_analysis/3_images/scan_dglm_bw3wk.pdf', height = height, width = width)

plot(x = sovs[['bw6wk']],     y = sos[['bw6wk']], plot.title = 'Bodyweight at 6 weeks')      # concordant plus two new vQTL
ggsave(filename = 'leamy_analysis/3_images/scan_dglm_bw6wk.pdf', height = height, width = width)

plot(x = sovs[['scfat12wk']], y = sos[['scfat12wk']], ymax = 4, plot.title = 'subcutaneous fat pad thickness at 12 weeks')  # concordant plus more evidence on 11 and 13 and one undiscovery on chr 15
ggsave(filename = 'leamy_analysis/3_images/scan_dglm_scfat12wk.pdf', height = height, width = width)

plot(x = sovs[['gfat12wk']],  y = sos[['gfat12wk']], ymax = 4, plot.title = 'gonadal fat pad thickness at 12 weeks')   # concordant (three QTL)
ggsave(filename = 'leamy_analysis/3_images/scan_dglm_gfat12wk.pdf', height = height, width = width)
