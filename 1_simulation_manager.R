library(rslurm)
library(dplyr)

source('Corty_simulation_function.R')

job_name <- 'fixed_covar_sampled_locus'


num_perms <- 1000 #1000
num_splits <- 100 #100
num_sims_per_split <- 100 #100

num_obs <- 300

if (num_obs == 300) {
  mqtl_pve <- 0.05
  vqtl_effect <- 0.17
  mvqtl_mean_pve <- 0.65*mqtl_pve
  mvqtl_var_effect <-  0.8*vqtl_effect
}

# if (num_obs == 1012) {
#   mqtl_pve <- 0.02
#   vqtl_effect <- 0.2
#   mvqtl_mean_pve <- 0.5*mqtl_pve
#   mvqtl_var_effect <-  0.8*vqtl_effect
# }


scenarios <- expand.grid(qtl = c('none', 'mqtl', 'vqtl', 'mvqtl'),
                         var_covar = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)

scenarios <- do.call(what = bind_rows, args = replicate(n = num_splits, expr = scenarios, simplify = FALSE))

# scenarios that differ only by whether covar is modeled have the same seed
scenarios$seed <- 1:(nrow(scenarios))


sjob <- slurm_apply(f = mvqtl_simulation,
                    params = scenarios,
                    jobname = paste0(job_name, '_', Sys.Date()),
                    nodes = nrow(scenarios),
                    cpus_per_node = 1,
                    slurm_options = list(time = "48:00:00",
                                         share = TRUE,
                                         mem = 2000),
                    add_objects = c('num_sims_per_split',
                                    'num_obs',
                                    'num_perms',
                                    'mqtl_pve',
                                    'vqtl_effect',
                                    'mvqtl_mean_pve',
                                    'mvqtl_var_effect'),
                    pkgs = c('tidyverse', 'nlme', 'car', 'dglm'))

res_raw <- get_slurm_out(sjob, outtype = "raw", wait = TRUE)

saveRDS(object = bind_rows(res_raw),
        file = paste0(job_name, '_',
                      Sys.Date(),
                      'nobs=', num_obs,
                      'nsim=', num_splits*num_sims_per_split,
                      '_nperm=', num_perms,
                      '.RDS'))

quit(save = 'no')
