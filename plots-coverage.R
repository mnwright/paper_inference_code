
library(data.table)
library(ggplot2)
library(patchwork)
library(foreach)

source("R/config.R")

coverage_pfi_mean = readRDS(sprintf("%s/coverage_pfi_mean.Rds", res_dir))
coverage_pdp_mean = readRDS(sprintf("%s/coverage_pdp_mean.Rds", res_dir))
coverage_shap_mean = readRDS(sprintf("%s/coverage_shap_mean.Rds", res_dir))

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

coverage_shap_mean[!train_missing & !test_missing, imputation_method := "None"]
coverage_shap_mean[!train_missing & !test_missing, missing_prob := 0]
coverage_shap_mean[!train_missing & !test_missing, pattern := "None"]
coverage_shap_mean[, missing := "None"]
coverage_shap_mean[train_missing & test_missing, missing := "Both"]
coverage_shap_mean[train_missing & !test_missing, missing := "Train"]
coverage_shap_mean[!train_missing & test_missing, missing := "Test"]
coverage_shap_mean[, missing := factor(missing)]

coverage_pfi_mean <- coverage_pfi_mean[, coverage := mean(coverage), 
                                       by = list(algorithm, problem, sampling_strategy, nrefits, 
                                                 adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)
                                       ]
coverage_pdp_mean <- coverage_pdp_mean[, coverage := mean(coverage), 
                                       by = list(algorithm, problem, sampling_strategy, nrefits, 
                                                 adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)
]
coverage_shap_mean <- coverage_shap_mean[, coverage := mean(coverage), 
                                       by = list(algorithm, problem, sampling_strategy, nrefits, 
                                                 adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)
]

plot_fun <- function(iml_method) {
  if (iml_method == "PFI") {
    coverage_mean = coverage_pfi_mean
  } else if (iml_method == "PDP") {
    coverage_mean = coverage_pdp_mean
  } else if (iml_method == "SHAP") {
    coverage_mean = coverage_shap_mean
  } else {
    stop("Unknown iml method")
  }
  
  pars <- expand.grid(n = unique(coverage_mean$n), 
                      sampling_strategy = as.character(unique(coverage_mean$sampling_strategy)), 
                      missing = setdiff(unique(coverage_mean$missing), "None"), 
                      pattern = setdiff(as.character(unique(coverage_mean$pattern)), "None"), 
                      algorithm = as.character(unique(coverage_mean$algorithm)),
                      #   missing_prob = setdiff(unique(coverage_mean$missing_prob), 0), 
                      stringsAsFactors = FALSE)
  mapply(function(xn, xmissing, xpattern, xalgorithm, xsampling_strategy) {
    ggplot(coverage_mean[n == xn & 
                           missing %in% c("None", xmissing) & 
                           pattern %in% c("None", xpattern) & 
                           algorithm %in% xalgorithm &
                           #  missing_prob %in% c(0, xmissing_prob) & 
                           sampling_strategy == xsampling_strategy, 
    ],
    aes(x = nrefits, y = coverage, color = imputation_method, 
        shape = adjusted, linetype = adjusted)) +
      facet_grid(problem ~ missing_prob) + 
      geom_line() + #geom_point() +
      scale_y_continuous(sprintf("Confidence Interval %s", "Coverage"), limits = c(0, 1)) +
      scale_x_continuous("Number of Model Refits") +
      scale_color_discrete("Imputation Method") + 
      geom_hline(yintercept = .95, linetype = "dashed") + 
      ggtitle(sprintf("%s (learner = %s, n = %s, sampling = %s, missing = %s, pattern = %s)", 
                      iml_method, xalgorithm, xn, xsampling_strategy, xmissing, xpattern))
  }, pars$n, pars$missing, pars$pattern, pars$algorithm, pars$sampling_strategy, SIMPLIFY = FALSE)
}

iml_methods <- c("PFI", "PDP", "SHAP")
pars <- data.table(expand.grid(iml_method = iml_methods, stringsAsFactors = FALSE))
#pars <- pars[!(iml_method == "SHAP" & learner == "randomForest"), ]
mapply(function(iml_method) {
  p <- plot_fun(iml_method)
  ggsave(sprintf("%s/coverage_%s.pdf", fig_dir, iml_method), wrap_plots(p, ncol = 2), width = 20, height = 30, 
         limitsize = FALSE)
}, pars$iml_method, SIMPLIFY = FALSE)

