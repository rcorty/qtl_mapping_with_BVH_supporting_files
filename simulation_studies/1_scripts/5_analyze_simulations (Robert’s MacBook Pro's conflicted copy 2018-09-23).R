library(dglm)
library(tidyverse)
library(plotROC)
library(RColorBrewer)
library(tables)
library(xtable)

latex_slashes <- function(c) { gsub(replacement = '\\', pattern = '\\textbackslash{}', x = c, fixed = TRUE) }
save_plots <- TRUE
max_alpha <- 0.15
qq_breaks <- 0.05
nominal_alpha <- 0.05

mqtl_test_names <- c('LM', 'Cao[m]', 'DGLM[m]')
vqtl_test_names <- c('Levene', 'Cao[v]', 'DGLM[v]')
mvqtl_test_names <- c('Cao[mv]', 'DGLM[mv]')

test_names <- c(mqtl_test_names, vqtl_test_names, mvqtl_test_names)

r <- readRDS('simulation_studies/3_results/bvh_sim_2018-09-17 22:06:52.RDS')

r %>% glimpse()
r %>%
  select(matches('num_perms')) %>%
  gather(key = test, val = num_perms) %>%
  ggplot(mapping = aes(x = num_perms, fill = test)) +
  geom_histogram(bins = 30, position = position_identity(), alpha = 0.3)

r %>%
  select(-matches('num_perms')) %>%
  glimpse()

r2 <- r %>%
  select(-matches('num_perms')) %>%
  gather(key = test, value = p, LM_asymp_ps:DGLMj_locusperm_ps) %>%
  separate(col = test, into = c('test', 'procedure', 'trash'), sep = '_') %>%
  select(-trash) %>%
  mutate(procedure = factor(procedure,
                            levels = c('asymp', 'rint', 'residperm', 'locusperm'),
                            labels = c('standard', 'RINT', 'residperm', 'locusperm')),
         bvh = factor(bvh,
                      levels = c('FALSE', 'TRUE'),
                      labels = c('BVH absent', 'BVH present')),
         qtl = factor(qtl, levels = c('none', 'mqtl', 'vqtl', 'mvqtl')),
         test = recode_factor(test,
                              LM = 'LM',
                              Caom = 'Cao[m]',
                              DGLMm = 'DGLM[m]',
                              Lev = 'Levene',
                              Caov = 'Cao[v]',
                              DGLMv = 'DGLM[v]',
                              Caoj = 'Cao[mv]',
                              DGLMj = 'DGLM[mv]'))

mqtl_tests <- r2 %>% filter(test %in% mqtl_test_names)
vqtl_tests <- r2 %>% filter(test %in% vqtl_test_names)
mvqtl_tests <- r2 %>% filter(test %in% mvqtl_test_names)

# mqtl_tests %>%
#   filter(qtl == 'none') %>%
#   ggplot(mapping = aes(x = p, fill = procedure)) +
#   geom_histogram(binwidth = 0.02, center = 0.01) +
#   facet_grid(bvh ~ test) +
#   theme_minimal()


#### ROC PLOTS ####
my_rocci_plot <- function(d) {
  d +
    geom_vline(xintercept = nominal_alpha) +
    geom_roc(mapping = aes(m = 1 - p,
                           d = as.numeric(qtl != 'none')),
             cutoffs.at = 1 - nominal_alpha,
             pointsize = 1.5,
             linealpha = 0.5,
             pointalpha = 0.7,
             labels = FALSE) +
    geom_rocci(mapping = aes(m = 1 - p,
                             d = as.numeric(qtl != 'none')),
               ci.at = 1 - nominal_alpha,
               labelsize = 0,
               alpha.box = 0.4) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text.align = 0,
          legend.title.align = 0.5) +
    labs(x = 'false positive rate', y = 'power') +
    coord_cartesian(xlim = c(0, max_alpha)) +
    scale_x_continuous(limits = c(0, max_alpha + 0.01)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}

my_roc_plot <- function(d) {
  d +
    geom_vline(xintercept = nominal_alpha) +
    geom_roc(mapping = aes(m = 1 - p,
                           d = as.numeric(qtl != 'none')),
             cutoffs.at = 1 - nominal_alpha,
             pointsize = 1.5,
             linealpha = 0.5,
             pointalpha = 0.7,
             labels = FALSE) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text.align = 0,
          legend.title.align = 0.5) +
    labs(x = 'false positive rate', y = 'power') +
    coord_cartesian(xlim = c(0, max_alpha)) +
    scale_x_continuous(limits = c(0, max_alpha + 0.01)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2))
}


facet_rocci_by_eval <- function(d) {
  d %>%
    ggplot(mapping = aes(color = test, shape = bvh)) %>%
    my_rocci_plot() +
    facet_grid(facets = ~ procedure) +
    scale_shape(name = 'BVH scenario') +
    scale_color_discrete(labels = scales::parse_format())
}

facet_rocci_by_bvh <- function(d) {
  d %>%
    ggplot(mapping = aes(color = test, shape = procedure)) %>%
    my_rocci_plot() +
    facet_grid(facets = ~ bvh) +
    scale_shape(name = 'procedure') +
    scale_color_discrete(labels = scales::parse_format())
}

facet_rocci_by_test <- function(d) {
  d %>%
    ggplot(mapping = aes(color = procedure, shape = bvh)) %>%
    my_rocci_plot() +
    facet_grid(facets = ~ test, labeller = label_parsed) +
    scale_shape(name = 'BVH scenario')
}

facet_roc_by_bvh <- function(d) {
  d %>%
    ggplot(mapping = aes(color = test, shape = procedure)) %>%
    my_roc_plot() +
    facet_grid(facets = ~ bvh) +
    scale_shape(name = 'procedure') +
    scale_color_discrete(labels = scales::parse_format())
}



# mqtl tests
mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_roc_by_bvh() +
  ggtitle('mQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/rocs_mqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_bvh() +
  ggtitle('mQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }


mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_eval() +
  ggtitle('mQTL discrimination by version')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_test() +
  ggtitle('mQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }



# vqtl tests
vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_roc_by_bvh() +
  ggtitle('vQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/rocs_vqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_bvh() +
  ggtitle('vQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_vqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_eval() +
  ggtitle('vQTL discrimination by version')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_vqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_test() +
  ggtitle('vQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_vqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }



# mvqtl tests
mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_roc_by_bvh() +
  ggtitle('mvQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/rocs_mvqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_bvh() +
  ggtitle('mvQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mvqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_eval() +
  ggtitle('mvQTL discrimination by version')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mvqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_rocci_by_test() +
  ggtitle('mvQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/roccis_mvqtl_all_facet_by_test.pdf', height = 3.5, width = 10) }





#### POSITIVE RATE TABLE ####
options(scipen = 10)
pos_rate_tab <-  r2 %>%
  group_by(bvh, qtl, test, procedure) %>%
  summarise(pos_rate = round(x = mean(p < nominal_alpha, na.rm = TRUE), digits = 3)) %>%
  tabular(table = test*procedure ~ Format(digits=2)*bvh*qtl*pos_rate*identity,
          data = .)
pos_rate_tab

xtable(x = as.matrix(pos_rate_tab)) %>%
  print(include.rownames = FALSE, sanitize.text.function = latex_slashes)


ci_tibble <-  r2 %>%
  mutate(bvh = factor(bvh, labels = c('BVH~absent', 'BVH~present'))) %>%
  group_by(bvh, qtl, test, procedure) %>%
  summarise(pos_count = sum(p < nominal_alpha, na.rm = TRUE),
            pos_rate = pos_count/n(),
            ci_lo = binom.test(x = pos_count,
                               n = n(),
                               p = pos_rate,
                               alternative = 'two.sided',
                               conf.level = 0.95)$conf.int[1],
            ci_hi = binom.test(x = pos_count,
                               n = n(),
                               p = pos_rate,
                               alternative = 'two.sided',
                               conf.level = 0.95)$conf.int[2]) %>%
  mutate(se = (ci_hi - ci_lo)/2,
         ideal_value = case_when(qtl == 'none' ~ 0.05,
                                 qtl == 'mqtl' ~ 0.717,
                                 qtl == 'vqtl' ~ 0.655,
                                 qtl == 'mvqtl' ~ 0.742))

ci_tibble

ci_tibble %>%
  filter(qtl == 'none') %>%
  ggplot() +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi),
               y = 0, yend = 0) +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh, qtl),
             switch = 'y',
             labeller = label_parsed) +
  theme_void() +
  ylim(c(-1, 1)) +
  xlab('false positive rate') +
  theme(panel.spacing = unit(x = 3, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/null_cis.pdf', height = 10, width = 5)

to_ci_zoom <- ci_tibble %>%
  filter(qtl == 'none') %>%
  mutate(range = case_when(ci_lo < 0.03 ~ 'lo',
                           ci_hi > 0.07 ~ 'hi',
                           TRUE ~ 'ok'),
         pos_rate_loc = case_when(ci_lo < 0.03 ~ 0.03,
                                  ci_hi > 0.07 ~ 0.071,
                                  TRUE ~ pos_rate))

to_ci_zoom %>%
  # filter(range == 'ok') %>%
  ggplot() +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi), y = 0, yend = 0) +
  geom_point(mapping = aes(x = pos_rate_loc, color = range), y = 0) +
  geom_text(data = to_ci_zoom %>% filter(range != 'ok'),
            mapping = aes(x = pos_rate_loc,
                          y = 0,
                          # color = range,
                          label = round(x = pos_rate, digits = 2)),
            size = 3) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh),
             switch = 'y',
             labeller = label_parsed) +
  scale_color_manual(values = c('gray95', 'gray95', 'black'),
                     guide = FALSE) +
  scale_x_continuous(name = 'false positive rate',
                     limits = c(0.027, 0.073),
                     breaks = 0.01*(3:7)) +
  theme_void() +
  ylim(c(-1, 1)) +
  theme(panel.spacing = unit(x = 3, units = 'pt'),
        # panel.border = element_rect(fill = NA, colour = 'black'),
        # panel.grid.major.x = element_line(colour = 'white'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        strip.text.x = element_text(size = 12, margin = margin(t = 2, r = 2, b = 4, l = 2, unit = 'pt')),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/null_cis_zoom.pdf', height = 10, width = 5)


ci_tibble %>%
  filter(qtl == 'mqtl', test %in% mqtl_test_names) %>%
  ggplot() +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi), y = 0, yend = 0) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh),
             switch = 'y',
             labeller = label_parsed) +
  theme_void() +
  ylim(c(-1, 1)) +
  scale_x_continuous(name = 'power',
                     limits = 0.01*(6.2:8.8),
                     breaks = 0.01*(6.5:8.5)) +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = 'pt'),
        panel.spacing = unit(x = 3, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        strip.text.x = element_text(size = 12, margin = margin(t = 2, r = 2, b = 4, l = 2, unit = 'pt')),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/mqtl_cis.pdf', height = 4.5, width = 5)


ci_tibble %>%
  filter(qtl == 'vqtl', test %in% vqtl_test_names) %>%
  ggplot() +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi),
               y = 0, yend = 0) +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh),
             switch = 'y',
             labeller = label_parsed) +
  theme_void() +
  ylim(c(-1, 1)) +
  xlim(c(0.5, 0.75)) +
  xlab('positive rate') +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = 'pt'),
        panel.spacing = unit(x = 3, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        strip.text.x = element_text(size = 12, margin = margin(t = 2, r = 2, b = 4, l = 2, unit = 'pt')),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/vqtl_cis.pdf', height = 4.5, width = 5)

ci_tibble %>%
  filter(qtl == 'mvqtl', test %in% mvqtl_test_names) %>%
  ggplot() +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi),
               y = 0, yend = 0) +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh),
             switch = 'y',
             labeller = label_parsed) +
  theme_void() +
  ylim(c(-1, 1)) +
  xlim(c(0.6, 0.85)) +
  xlab('positive rate') +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = 'pt'),
        panel.spacing = unit(x = 3, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        strip.text.x = element_text(size = 12, margin = margin(t = 2, r = 2, b = 4, l = 2, unit = 'pt')),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/mvqtl_cis.pdf', height = 3.5, width = 5)




ci_tibble %>%
  filter(qtl %in% c('none', 'mqtl'), test %in% mqtl_test_names, bvh == 'BVH~absent') %>%
  ggplot() +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi), y = 0, yend = 0) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh, qtl),
             switch = 'y',
             scales = 'free_x') +
  theme_void() +
  ylim(c(-1, 1)) +
  # xlim(c(0.65, 0.85)) +
  xlab('positive rate') +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = 'pt'),
        panel.spacing = unit(x = 5, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/mqtl_FPR_and_power_cis.pdf', height = 5, width = 8)


ci_tibble %>%
  filter(qtl %in% c('none', 'vqtl'), test %in% vqtl_test_names) %>%
  ggplot() +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi), y = 0, yend = 0) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh, qtl),
             switch = 'y',
             scales = 'free_x') +
  theme_void() +
  ylim(c(-1, 1)) +
  # xlim(c(0.65, 0.85)) +
  xlab('positive rate') +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = 'pt'),
        panel.spacing = unit(x = 5, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/vqtl_FPR_and_power_cis.pdf', height = 5, width = 8)



ci_tibble %>%
  filter(qtl %in% c('none', 'mvqtl'), test %in% mvqtl_test_names) %>%
  ggplot() +
  geom_vline(mapping = aes(xintercept = ideal_value), color = 'gray', size = 0.5) +
  geom_segment(mapping = aes(x = ci_lo, xend = ci_hi), y = 0, yend = 0) +
  geom_point(mapping = aes(x = pos_rate), y = 0) +
  facet_grid(rows = vars(test, procedure),
             cols = vars(bvh, qtl),
             switch = 'y',
             scales = 'free_x') +
  theme_void() +
  ylim(c(-1, 1)) +
  # xlim(c(0.65, 0.85)) +
  xlab('positive rate') +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = 'pt'),
        panel.spacing = unit(x = 5, units = 'pt'),
        axis.title.x = element_text(size = unit(x = 12, units = 'pt'), hjust = 0.5),
        axis.text.x = element_text(size = unit(x = 8, units = 'pt')),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 180),
        panel.background = element_rect(fill = 'gray95')) +
  ggsave(filename = 'simulation_studies/4_images/mvqtl_FPR_and_power_cis.pdf', height = 5, width = 8)



ci_lo_tab <- ci_tibble %>%
  tabular(table = test*procedure ~ Format(digits=2)*bvh*qtl*ci_lo*identity,
          data = .)
ci_lo_tab

ci_hi_tab <- ci_tibble %>%
  tabular(table = test*procedure ~ Format(digits=2)*bvh*qtl*ci_hi*identity,
          data = .)
ci_hi_tab

se_tab <- ci_tibble %>%
  tabular(table = test*procedure ~ Format(digits=2)*bvh*qtl*se*identity,
          data = .)
se_tab

xtable(x = as.matrix(se_tab)) %>%
  print(include.rownames = FALSE, sanitize.text.function = latex_slashes)

ci_tibble %>%
  ungroup() %>%
  group_by(bvh, qtl) %>%
  summarise(max(se))



#### AUC TABLE ####
nulls <- r2 %>% filter(qtl == 'none') %>% mutate(null_p = 1-p) %>% select(-p, -qtl)
alts <- r2 %>% filter(qtl != 'none') %>% mutate(alt_p = 1-p) %>% select(-p)

to_auc <- inner_join(x = alts, y = nulls, by = c('job_idx', 'sim_idx', 'bvh', 'test', 'procedure'))

auc_tab <- to_auc %>%
  group_by(bvh, qtl, test, procedure) %>%
  summarise(AUC = MESS::auc(x = null_p, y = alt_p)) %>%
  ungroup() %>%
  mutate(qtl = factor(qtl)) %>%
  tabular(table = test*procedure ~ Format(digits=3)*bvh*qtl*AUC*identity,
          data = .)

auc_tab

xtable(x = as.matrix(auc_tab)) %>%
  print(include.rownames = FALSE, sanitize.text.function = latex_slashes)



#### QQ PLOTS for supplement ####
r3 <- r2 %>% mutate(bvh = factor(bvh, labels = c('BVH~absent', 'BVH~present')))

my_qq <- function(d) {
  d %>%
    ggplot(mapping = aes(sample = p, color = procedure)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_qq(distribution = qunif, geom = 'line', size = 1, alpha = 0.8) +
    coord_flip() +
    facet_grid(test ~ bvh, labeller = label_parsed) +
    scale_x_continuous(name = 'observed false positive rate (FPR)', limits = c(0, max_alpha), breaks = seq(from = 0, to = max_alpha, by = qq_breaks)) +
    scale_y_continuous(name = expression('nominal false positive rate ('~alpha~')'), limits = c(0, max_alpha), breaks = seq(from = 0, to = max_alpha, by = qq_breaks)) +
    scale_color_discrete(name = 'version') +
    theme_minimal() +
    theme(panel.grid.major = element_line(),
          panel.grid.minor = element_blank())
}

# all in one -- too much info for one plot
# r3 %>%
#   filter(qtl == 'none') %>%
#   my_qq()
# if (save_plots) { ggsave(filename = 'null_qqs.pdf', height = 12, width = 6) }


r3 %>%
  filter(test %in% mqtl_test_names, qtl == 'none') %>%
  my_qq()
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/mqtl_null_qqs.pdf', height = 9, width = 7) }


r3 %>%
  filter(test %in% vqtl_test_names, qtl == 'none') %>%
  my_qq()
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/vqtl_null_qqs.pdf', height = 9, width = 7) }


r3 %>%
  filter(test %in% mvqtl_test_names, qtl == 'none') %>%
  my_qq()
if (save_plots) { ggsave(filename = 'simulation_studies/4_images/mvqtl_null_qqs.pdf', height = 7, width = 7) }









# #### simpler ROCs for main paper
# # collapse DGLM down to situations where DGLM matches existence of BVH
#
# mqtl_test_names <- c('LM', 'Cao[m]', 'DGLM[m]')
# vqtl_test_names <- c('Levene', 'Cao[v]', 'DGLM[v]')
# mvqtl_test_names <- c('Cao[mv]', 'DGLM[mv]')
#
# test_names <- c(mqtl_test_names, vqtl_test_names, mvqtl_test_names)
#
# r3 <- r2 %>%
#   filter(!(grepl(pattern = 'no~covar', x = test))) %>%
#   mutate(test = as.character(test),
#          test = case_when(test == 'DGLM[m]~covar' ~ 'DGLM[m]',
#                           test == 'DGLM[v]~covar' ~ 'DGLM[v]',
#                           test == 'DGLM[mv]~covar' ~ 'DGLM[mv]',
#                           TRUE ~ test),
#          test = factor(test, levels = test_names))
#
#
# mqtl_tests <- r3 %>% filter(test %in% mqtl_test_names)
# vqtl_tests <- r3 %>% filter(test %in% vqtl_test_names)
# mvqtl_tests <- r3 %>% filter(test %in% mvqtl_test_names)
#
#
# mqtl_tests %>%
#   filter(qtl == 'mqtl') %>%
#   bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
#   facet_by_bvh() +
#   ggtitle('mQTL discrimination by test')
# if (save_plots) { ggsave(filename = 'simple_rocs_mqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }
#
#
# vqtl_tests %>%
#   filter(qtl == 'vqtl') %>%
#   bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
#   facet_by_bvh() +
#   ggtitle('vQTL discrimination by test')
# if (save_plots) { ggsave(filename = 'simple_rocs_vqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }
#
# mvqtl_tests %>%
#   filter(qtl == 'mvqtl') %>%
#   bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
#   facet_by_bvh() +
#   ggtitle('mvQTL discrimination by test')
# if (save_plots) { ggsave(filename = 'simple_rocs_mvqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }
