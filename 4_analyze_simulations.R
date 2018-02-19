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

mqtl_test_names <- c('LM', 'Cao[m]', 'DGLM[m]~no~covar', 'DGLM[m]~covar')
vqtl_test_names <- c('Levene', 'Cao[v]', 'DGLM[v]~no~covar', 'DGLM[v]~covar')
mvqtl_test_names <- c('Cao[mv]', 'DGLM[mv]~no~covar', 'DGLM[mv]~covar')

test_names <- c(mqtl_test_names, vqtl_test_names, mvqtl_test_names)


# r <- readRDS('sim_2018-01-21 13:27:39.RDS')
# r <- readRDS('fixed_covar_sampled_locus_2018-02-07nobs=300nsim=10000_nperm=0.RDS')
r <- readRDS('sim_2018-02-14 13:38:14.RDS')

r %>% glimpse()
r %>%
  select(matches('num_perms')) %>%
  gather(key = test, val = num_perms) %>%
  ggplot(mapping = aes(x = num_perms, fill = test)) +
  geom_histogram(bins = 30, position = position_identity(), alpha = 0.3)


r2 <- r %>%
  select(-matches('num_perms')) %>%
  gather(key = test, value = p, LM__asymp_ps:DGLMj_ye_covar__locusperm_ps) %>%
  separate(col = test, into = c('test', 'procedure'), sep = '__') %>%
  mutate(procedure = factor(procedure,
                             levels = c('asymp_ps', 'rint_ps', 'residperm_ps', 'locusperm_ps'),
                             labels = c('standard', 'RINT', 'residperm', 'locusperm')),
         bvh = factor(var_covar,
                      levels = c('FALSE', 'TRUE'),
                      labels = c('BVH absent', 'BVH present')),
         qtl = factor(qtl, levels = c('none', 'mqtl', 'vqtl', 'mvqtl')),
         model_covar = ifelse(test = grepl(pattern = 'ye_covar', x = test),
                              yes = TRUE,
                              no = FALSE),
         test = recode_factor(test,
                              LM = 'LM',
                              Caom = 'Cao[m]',
                              DGLMm_no_covar = 'DGLM[m]~no~covar',
                              DGLMm_ye_covar = 'DGLM[m]~covar',
                              Lev = 'Levene',
                              Caov = 'Cao[v]',
                              DGLMv_no_covar = 'DGLM[v]~no~covar',
                              DGLMv_ye_covar = 'DGLM[v]~covar',
                              Caoj = 'Cao[mv]',
                              DGLMj_no_covar = 'DGLM[mv]~no~covar',
                              DGLMj_ye_covar = 'DGLM[mv]~covar')) %>%
  select(-var_covar)


mqtl_tests <- r2 %>% filter(test %in% mqtl_test_names)
vqtl_tests <- r2 %>% filter(test %in% vqtl_test_names)
mvqtl_tests <- r2 %>% filter(test %in% mvqtl_test_names)


#### ROC PLOTS ####
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

facet_by_eval <- function(d) {
  d %>%
    ggplot(mapping = aes(color = test, shape = bvh)) %>%
    my_roc_plot() +
    facet_grid(facets = ~ procedure) +
    scale_shape(name = 'BVH scenario') +
    scale_color_discrete(labels = scales::parse_format())

}

facet_by_bvh <- function(d) {
  d %>%
    ggplot(mapping = aes(color = test, shape = procedure)) %>%
    my_roc_plot() +
    facet_grid(facets = ~ bvh) +
    scale_shape(name = 'procedure') +
    scale_color_discrete(labels = scales::parse_format())
}

facet_by_test <- function(d) {
  d %>%
    ggplot(mapping = aes(color = procedure, shape = bvh)) %>%
    my_roc_plot() +
    facet_grid(facets = ~ test, labeller = label_parsed) +
    scale_shape(name = 'BVH scenario')
}


# mqtl tests
mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_eval() +
  ggtitle('mQTL discrimination by version')
if (save_plots) { ggsave(filename = 'rocs_mqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_bvh() +
  ggtitle('mQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'rocs_mqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('mQTL discrimination by test')
if (save_plots) { ggsave(filename = 'rocs_mqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }



# vqtl tests
vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_eval() +
  ggtitle('vQTL discrimination by version')
if (save_plots) { ggsave(filename = 'rocs_vqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_bvh() +
  ggtitle('vQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'rocs_vqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('vQTL discrimination by test')
if (save_plots) { ggsave(filename = 'rocs_vqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }



# mvqtl tests
mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_eval() +
  ggtitle('mvQTL discrimination by version')
if (save_plots) { ggsave(filename = 'rocs_mvqtl_all_facet_by_eval.pdf', height = 3.5, width = 13) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_bvh() +
  ggtitle('mvQTL discrimination by BVH scenario')
if (save_plots) { ggsave(filename = 'rocs_mvqtl_all_facet_by_bvh.pdf', height = 3.5, width = 7) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('mvQTL discrimination by test')
if (save_plots) { ggsave(filename = 'rocs_mvqtl_all_facet_by_test.pdf', height = 3.5, width = 10) }





#### POSITIVE RATE TABLE ####
big_tab <-  r2 %>%
  group_by(bvh, qtl, test, procedure) %>%
  summarise(pos_rate = mean(p < nominal_alpha, na.rm = TRUE)) %>%
  tabular(table = test*procedure ~ Format(digits=2)*bvh*qtl*pos_rate*identity,
          data = .)
big_tab

xtable(x = as.matrix(big_tab)) %>%
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
if (save_plots) { ggsave(filename = 'mqtl_null_qqs.pdf', height = 7, width = 7) }


r3 %>%
  filter(test %in% vqtl_test_names, qtl == 'none') %>%
  my_qq()
if (save_plots) { ggsave(filename = 'vqtl_null_qqs.pdf', height = 7, width = 7) }


r3 %>%
  filter(test %in% mvqtl_test_names, qtl == 'none') %>%
  my_qq()
if (save_plots) { ggsave(filename = 'mvqtl_null_qqs.pdf', height = 5, width = 7) }









#### simpler ROCs for main paper
# collapse DGLM down to situations where DGLM matches existence of BVH

mqtl_test_names <- c('LM', 'Cao[m]', 'DGLM[m]')
vqtl_test_names <- c('Levene', 'Cao[v]', 'DGLM[v]')
mvqtl_test_names <- c('Cao[mv]', 'DGLM[mv]')

test_names <- c(mqtl_test_names, vqtl_test_names, mvqtl_test_names)

r3 <- r2 %>%
  filter(!(grepl(pattern = 'no~covar', x = test) & bvh == 'BVH present'), !(grepl(pattern = ']~covar', x = test) & bvh == 'BVH absent')) %>%
  mutate(test = as.character(test),
         test = case_when(test == 'DGLM[m]~no~covar' ~ 'DGLM[m]',
                          test == 'DGLM[m]~covar' ~ 'DGLM[m]',
                          test == 'DGLM[v]~no~covar' ~ 'DGLM[v]',
                          test == 'DGLM[v]~covar' ~ 'DGLM[v]',
                          test == 'DGLM[mv]~no~covar' ~ 'DGLM[mv]',
                          test == 'DGLM[mv]~covar' ~ 'DGLM[mv]',
                          TRUE ~ test),
         test = factor(test, levels = test_names))


mqtl_tests <- r3 %>% filter(test %in% mqtl_test_names)
vqtl_tests <- r3 %>% filter(test %in% vqtl_test_names)
mvqtl_tests <- r3 %>% filter(test %in% mvqtl_test_names)


mqtl_tests %>%
  filter(qtl == 'mqtl') %>%
  bind_rows(mqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('mQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simple_rocs_mqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }


vqtl_tests %>%
  filter(qtl == 'vqtl') %>%
  bind_rows(vqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('vQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simple_rocs_vqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }

mvqtl_tests %>%
  filter(qtl == 'mvqtl') %>%
  bind_rows(mvqtl_tests %>% filter(qtl == 'none')) %>%
  facet_by_test() +
  ggtitle('mvQTL discrimination by test')
if (save_plots) { ggsave(filename = 'simple_rocs_mvqtl_all_facet_by_test.pdf', height = 3.5, width = 13) }
