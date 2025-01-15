
library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)

source("R/config.R")

coverage_pfi_mean = readRDS(sprintf("%s/coverage_pfi_mean.Rds", res_dir))
coverage_pdp_mean = readRDS(sprintf("%s/coverage_pdp_mean.Rds", res_dir))

coverage_pfi_mean[!train_missing & !test_missing, imputation_method := "None"]
coverage_pfi_mean[!train_missing & !test_missing, missing_prob := 0]
coverage_pfi_mean[!train_missing & !test_missing, pattern := "None"]
coverage_pfi_mean[, missing := "None"]
coverage_pfi_mean[train_missing & test_missing, missing := "Both"]
coverage_pfi_mean[train_missing & !test_missing, missing := "Train"]
coverage_pfi_mean[!train_missing & test_missing, missing := "Test"]
coverage_pfi_mean[, missing := factor(missing)]

coverage_pdp_mean[!train_missing & !test_missing, imputation_method := "None"]
coverage_pdp_mean[!train_missing & !test_missing, missing_prob := 0]
coverage_pdp_mean[!train_missing & !test_missing, pattern := "None"]
coverage_pdp_mean[, missing := "None"]
coverage_pdp_mean[train_missing & test_missing, missing := "Both"]
coverage_pdp_mean[train_missing & !test_missing, missing := "Train"]
coverage_pdp_mean[!train_missing & test_missing, missing := "Test"]
coverage_pdp_mean[, missing := factor(missing)]

# Plot PFI
pars <- expand.grid(n = unique(coverage_pfi_mean$n), 
                    missing = setdiff(unique(coverage_pfi_mean$missing), "None"), 
                    pattern = setdiff(as.character(unique(coverage_pfi_mean$pattern)), "None"), 
                    missing_prob = setdiff(unique(coverage_pfi_mean$missing_prob), 0), 
                    sampling_strategy = c("ideal", "subsampling", "bootstrap"), #as.character(unique(coverage_pfi_mean$sampling_strategy)), 
                    stringsAsFactors = FALSE)
plots_pfi <- mapply(function(xn, xmissing, xpattern, xmissing_prob, xsampling_strategy) {
  ggplot(coverage_pfi_mean[n == xn & 
                           missing %in% c("None", xmissing) & 
                           pattern %in% c("None", xpattern) & 
                           missing_prob %in% c(0, xmissing_prob) & 
                           sampling_strategy == xsampling_strategy, 
                           ],
         aes(x = nrefits, y = coverage, color = imputation_method, 
             shape = adjusted, linetype = adjusted)) +
    facet_grid(problem ~ algorithm) + 
    geom_line() + 
    scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
    scale_x_continuous("Number of Model Refits") +
    scale_color_discrete("Imputation Method") + 
    geom_hline(yintercept = .95, linetype = "dashed") + 
    ggtitle(sprintf("PFI (n = %s, sampling = %s, missing = %s, pattern = %s, missing_prob = %s)", 
                    xn, xsampling_strategy, xmissing, xpattern, xmissing_prob))
}, pars$n, pars$missing, pars$pattern, pars$missing_prob, pars$sampling_strategy, SIMPLIFY = FALSE)

# Plot together
ggsave("pfi_coverage_missings.pdf", wrap_plots(plots_pfi, ncol = 2), width = 20, height = 500, 
       limitsize = FALSE)

# Plot PDP
pars <- expand.grid(n = unique(coverage_pdp_mean$n), 
                    missing = setdiff(unique(coverage_pdp_mean$missing), "None"), 
                    pattern = setdiff(as.character(unique(coverage_pdp_mean$pattern)), "None"), 
                    missing_prob = setdiff(unique(coverage_pdp_mean$missing_prob), 0), 
                    sampling_strategy = c("ideal", "subsampling", "bootstrap"), #as.character(unique(coverage_pdp_mean$sampling_strategy)), 
                    stringsAsFactors = FALSE)
plots_pdp <- mapply(function(xn, xmissing, xpattern, xmissing_prob, xsampling_strategy) {
  ggplot(coverage_pdp_mean[n == xn & 
                             missing %in% c("None", xmissing) & 
                             pattern %in% c("None", xpattern) & 
                             missing_prob %in% c(0, xmissing_prob) & 
                             sampling_strategy == xsampling_strategy, 
  ],
  aes(x = nrefits, y = coverage, color = imputation_method, 
      shape = adjusted, linetype = adjusted)) +
    facet_grid(problem ~ algorithm) + 
    geom_line() + 
    scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
    scale_x_continuous("Number of Model Refits") +
    scale_color_discrete("Imputation Method") + 
    geom_hline(yintercept = .95, linetype = "dashed") + 
    ggtitle(sprintf("PDP (n = %s, sampling = %s, missing = %s, pattern = %s, missing_prob = %s)", 
                    xn, xsampling_strategy, xmissing, xpattern, xmissing_prob))
}, pars$n, pars$missing, pars$pattern, pars$missing_prob, pars$sampling_strategy, SIMPLIFY = FALSE)

# Plot together
ggsave("pdp_coverage_missings.pdf", wrap_plots(plots_pdp, ncol = 2), width = 20, height = 500, 
       limitsize = FALSE)
