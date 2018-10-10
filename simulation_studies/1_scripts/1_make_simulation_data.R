library(tidyverse)

set.seed(27599)
options(scipen = 10)

num_sims <- 1e2
# num_sims <- 1e3
# num_sims <- 1e4

num_obs <- 300

mqtl_pve <- 0.05
vqtl_effect <- 0.17
mvqtl_mean_pve <- 0.65*mqtl_pve
mvqtl_var_effect <-  0.8*vqtl_effect

covar <- rep(x = 1:5, times = num_obs/5)
bvh_effect <- seq(from = -0.4, to = 0.4, length.out = 5)[covar]

scenarios <- crossing(sim_idx = 1:num_sims,
                      bvh = c(FALSE, TRUE),
                      qtl = c('none', 'mqtl', 'vqtl', 'mvqtl'))

scenarios


designs <- data_frame(sim_idx = 1:num_sims,
                      covar = list(covar),
                      locus = list(sample(x = c(-1L, 0L, 1L),
                                          size = num_obs,
                                          replace = TRUE,
                                          prob = c(0.25, 0.5, 0.25))),
                      epsilon = list(rnorm(n = num_obs)))

designs


sim_phenos <- function(bvh, qtl, covar, locus, epsilon) {

  locus_var_linpred <- locus*switch(EXPR = qtl,
                                    'none' = 0,
                                    'mqtl' = 0,
                                    'vqtl' = vqtl_effect,
                                    'mvqtl' = mvqtl_var_effect,
                                    NULL)

  var_linpred <- if (bvh) {
    locus_var_linpred + bvh_effect
  } else {
    locus_var_linpred
  }

  resid_vars <- exp(2*var_linpred)
  average_resid_var <- mean(resid_vars)

  mean_linpred <- locus*switch(EXPR = qtl,
                               'none' = 0,
                               'mqtl' = sqrt(mqtl_pve*average_resid_var/(1-mqtl_pve)),
                               'vqtl' = 0,
                               'mvqtl' = sqrt(mvqtl_mean_pve*average_resid_var/(1-mvqtl_mean_pve)),
                               NULL)

  phen <- mean_linpred + epsilon*sqrt(resid_vars)

  return(phen)
}

inner_join(scenarios, designs) %>%
  group_by(sim_idx, bvh, qtl) %>%
  mutate(phen = pmap(.l = list(covar, locus, epsilon),
                     .f = sim_phenos,
                     bvh = bvh,
                     qtl = qtl)) ->
  simulated_data_nested

simulated_data_nested

simulated_data_unnested <- unnest(simulated_data_nested)

print(simulated_data_unnested)

# simulated_data_unnested %>%
#   filter(sim_idx == 1, qtl == 'none', bvh == FALSE) %>%
#   ggplot(mapping = aes(x = phen, fill = factor(locus))) +
#   geom_density(alpha = 0.3)
#
# simulated_data_unnested %>%
#   filter(sim_idx == 1, qtl == 'mqtl', bvh == FALSE) %>%
#   ggplot(mapping = aes(x = phen, fill = factor(locus))) +
#   geom_density(alpha = 0.3)
#
# simulated_data_unnested %>%
#   filter(sim_idx == 1, qtl == 'vqtl', bvh == FALSE) %>%
#   ggplot(mapping = aes(x = phen, fill = factor(locus))) +
#   geom_density(alpha = 0.3)
#
# simulated_data_unnested %>%
#   filter(sim_idx == 1, qtl == 'mvqtl', bvh == FALSE) %>%
#   ggplot(mapping = aes(x = phen, fill = factor(locus))) +
#   geom_density(alpha = 0.3)
#
# simulated_data_nested %>%
#   group_by(sim_idx) %>%
#   count()
#
# simulated_data_nested %>%
#   group_by(sim_idx, bvh, qtl) %>%
#   count()


write_csv(x = simulated_data_unnested,
          path = paste0('simulation_studies/2_data/',
                        'simulated_data_nsims=',
                        num_sims,
                        '.csv'))

write_rds(x = simulated_data_nested,
          path = paste0('simulation_studies/2_data/',
                        'simulated_data_nsims=',
                        num_sims,
                        '.gz'),
          compress = 'gz')
