library(rslurm)
library(dplyr)

source('simulation_studies/1_scripts/3_simulation.R')

job_name <- 'bvh_sims'


num_jobs <- 1000 #100
num_sims_per_job <- 10 #100
num_perms <- 1000


num_sims <- num_jobs*num_sims_per_job


# # local test
# library(tidyverse); library(nlme); library(car); library(dglm);
# r <- do.call(what = mvqtl_simulation,
#              args = data_frame(job_idx = 1))



sjob <- slurm_apply(f = mvqtl_simulation,
                    params = data_frame(job_idx = 1:num_jobs),
                    jobname = paste0(job_name, '_', Sys.Date()),
                    nodes = num_jobs,
                    cpus_per_node = 1,
                    slurm_options = list(time = '48:00:00',
                                         share = TRUE,
                                         mem = 2000),
                    add_objects = c('num_jobs',
                                    'num_sims_per_job',
                                    'num_sims',
                                    'num_perms'),
                    pkgs = c('tidyverse', 'nlme', 'car', 'dglm'))

res_raw <- get_slurm_out(slr_job = sjob,
                         outtype = "raw",
                         wait = TRUE)

to_save <- bind_rows(res_raw)

glimpse(to_save)

saveRDS(object = to_save,
        file = paste0('simulation_studies/3_results/',
                      job_name, '_',
                      Sys.Date(),
                      '_num_jobs=', num_jobs,
                      '_num_sims_per_job=', num_sims_per_job,
                      '_num_perms=', num_perms,
                      '.RDS'))

quit(save = 'no')
