
library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)

source("R/config.R")

perf <- readRDS(sprintf("%s/perfs-experiments.Rds", res_dir))

perf[!train_missing & !test_missing, imputation_method := "None"]
perf[!train_missing & !test_missing, missing_prob := 0]
perf[!train_missing & !test_missing, pattern := "None"]
perf[, missing := "None"]
perf[train_missing & test_missing, missing := "Both"]
perf[train_missing & !test_missing, missing := "Train"]
perf[!train_missing & test_missing, missing := "Test"]
perf[, missing := factor(missing)]

perf_mean <- perf[, .(mse = mean(mse)), by = .(n, problem, algorithm, sampling_strategy, missing,
                                               missing_prob, pattern, imputation_method)]


pars <- expand.grid(n = unique(perf_mean$n), 
                    sampling_strategy = as.character(unique(perf_mean$sampling_strategy)), 
                    missing = setdiff(unique(perf_mean$missing), "None"), 
                    pattern = setdiff(as.character(unique(perf_mean$pattern)), "None"), 
                    #   missing_prob = setdiff(unique(perf_mean$missing_prob), 0), 
                    stringsAsFactors = FALSE)
plots <- mapply(function(xn, xmissing, xpattern, xsampling_strategy) {
  ggplot(perf_mean[n == xn & 
                         missing %in% c("None", xmissing) & 
                         pattern %in% c("None", xpattern) & 
                         #  missing_prob %in% c(0, xmissing_prob) & 
                         sampling_strategy == xsampling_strategy, 
  ],
  aes(x = missing_prob, y = mse, color = imputation_method)) +
    facet_grid(problem ~ algorithm) + 
    geom_line() + 
    geom_point() + 
    #scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
    #scale_x_continuous("Number of Model Refits") +
    scale_color_discrete("Imputation Method") + 
    #geom_hline(yintercept = 0) + 
    ggtitle(sprintf("Performance (n = %s, sampling = %s, missing = %s, pattern = %s)", 
                    xn, xsampling_strategy, xmissing, xpattern))
}, pars$n, pars$missing, pars$pattern, pars$sampling_strategy, SIMPLIFY = FALSE)

ggsave(sprintf("%s/performance.pdf", fig_dir), wrap_plots(plots, ncol = 2), width = 20, height = 20, 
       limitsize = FALSE)
