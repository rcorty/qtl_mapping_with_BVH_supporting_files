library(dplyr)

fns <- list.files(path = '_rslurm_fixed_covar_sampled_locus_20180206',
		full.names = TRUE,
		pattern = 'results')

#fns <- paste0('_rslurm_20170627_sims/results_', 0:1199,'.RDS')


r <- list()
for (fn in fns) {
  try(r[[fn]] <- bind_rows(readRDS(file = fn)))
  message(fn, ' has size ', dim(r[[fn]]))
}

result <- bind_rows(r)

saveRDS(object = result,
        file = paste0('sim_', Sys.time(), '.RDS'))
