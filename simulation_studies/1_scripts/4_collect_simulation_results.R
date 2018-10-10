library(dplyr)

message('Collecting simulation results...')

fns <- list.files(path = '_rslurm_bvh_sims_20180916',
                  full.names = TRUE,
                  pattern = 'results')

message('File names are...', fns)



r <- list()
for (fn in fns) {
  try(r[[fn]] <- bind_rows(readRDS(file = fn)))
  message(fn, ' has size ', nrow(r[[fn]]), ' by ', ncol(r[[fn]]))
}

result <- bind_rows(r)

saveRDS(object = result,
        file = paste0('simulation_studies/3_results/bvh_sim_', Sys.time(), '.RDS'))
