# =============================================================================
# This script implements the CI coverage experiment
# =============================================================================


# =============================================================================
# Experiment Settings
# =============================================================================
N_EXPERIMENTS = 100 #10000
N_TRAIN = 1000
#N_SAMPLE = 100
MAX_REFITS = 30
# Number of permutations for PFI computation
N_PERM = 5
# Number of refits to estimate true PFI/PDP
# This should be set high
N_TRUE = 10000
NC = 20

# This loads all dependencies and utility functions
devtools::load_all()
set.seed(1)

# Clean up
unlink("registry", recursive = TRUE)

if(file.exists("registry")) {
  reg = loadRegistry("registry", writeable = TRUE)
} else {
  reg = makeExperimentRegistry(file.dir = "registry", source = "source.R")
}
#clearRegistry(reg)

# see file R/ci-experiment.R
addProblem(name = "x12", data = data.frame(), fun = gdata12, seed = 1)
addProblem(name = "x1234", data = data.frame(), fun = gdata, seed = 1)

# see file R/ci-experiment.R
addAlgorithm(name = "lm",  lm_wrapper)
addAlgorithm(name = "rpart",  rpart_wrapper)
addAlgorithm(name = "randomForest",  rf_wrapper)

strgs = c("bootstrap", "ideal") #c("subsampling", "bootstrap", "ideal")
missing_probs <- c(0.1, 0.2, 0.4) # Proportion of missing data
patterns <- c("MCAR", "MAR", "MNAR") # Missing data patterns
train_missing <- c(TRUE, FALSE) # Missing data in training data
test_missing <- c(TRUE, FALSE) # Missing data in test data
imputation_methods <- c("mean", "missForest", "mice")
setting = expand.grid(n = N_TRAIN,
                      max_refits = MAX_REFITS,
                      n_perm = N_PERM,
                      sampling_strategy = strgs, 
                      missing_prob = missing_probs, 
                      pattern = patterns, 
                      train_missing = train_missing,
                      test_missing = test_missing, 
                      imputation_method = imputation_methods
                      )
setting <- data.table(setting)

# Only consider missings in both training and test (or none at all)
setting <- setting[(!train_missing & !test_missing) | (train_missing & test_missing), ]

# Remove duplicated experiments (no missing data)
setting <- setting[!(!train_missing & !test_missing & (missing_prob != 0.1 | pattern != "MCAR" | imputation_method != "mean")), ]

pdes = list(x12 = setting, x1234 = setting)

addExperiments(pdes, repls = N_EXPERIMENTS)
summarizeExperiments()


# =============================================================================
# Run Experiment
# =============================================================================

#testJob(1)
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = NC, fs.latency = 0)
submitJobs()
waitForJobs()

pars = unwrap(getJobPars())



combine_pfis = function(a, b) {append(a, b["pfi"])}
pfis = reduceResults(fun = combine_pfis, init = list())
pfis = rbindlist(pfis)
pfis = ijoin(pars, pfis, by = "job.id")
saveRDS(pfis, sprintf("%s/pfis-experiments.Rds", res_dir))
gc()

combine_pdps = function(a, b) {append(a, b["pdp"])}
pdps = reduceResults(fun = combine_pdps, init = list())
pdps = rbindlist(pdps)
pdps = ijoin(pars, pdps, by = "job.id")
saveRDS(pdps, sprintf("%s/pdps-experiments.Rds", res_dir))
gc()

combine_shaps = function(a, b) {append(a, b["shap"])}
shaps = reduceResults(fun = combine_shaps, init = list())
shaps = rbindlist(shaps)
shaps = ijoin(pars, shaps, by = "job.id")
saveRDS(shaps, sprintf("%s/shaps-experiments.Rds", res_dir))
gc()

combine_perf = function(a, b) {append(a, b["perf"])}
perfs = reduceResults(fun = combine_perf, init = list())
perfs = rbindlist(perfs)
perfs = ijoin(pars, perfs, by = "job.id")
saveRDS(perfs, sprintf("%s/perfs-experiments.Rds", res_dir))
gc()

# =============================================================================
# Compute true importance and PDP
# =============================================================================

tpfis = list()
tpdps = list()
tshaps = list()

tpfi_file = sprintf("%s/tpfis.Rds", res_dir)
tpdp_file = sprintf("%s/tpdps.Rds", res_dir)
tshap_file = sprintf("%s/tshaps.Rds", res_dir)

if (!file.exists(tpfi_file)){
  mod_names = c("lm", "rpart", "randomForest")
  dgp_names = c("x12", "x1234")
  for(i in mod_names) {
    for (j in dgp_names) {
      for (ntrain in N_TRAIN) {
        message(i, j, ntrain)
        tpfi = get_true_pfi(N_TRUE, ntrain = ntrain, dgps[[j]], mods[[i]])
        tpdp = get_true_pdp(N_TRUE, ntrain = ntrain, dgps[[j]], mods[[i]])
        tshap = get_true_shap(N_TRUE, ntrain = ntrain, dgps[[j]], mods[[i]])
        tpfi$algorithm = tpdp$algorithm = tshap$algorithm = i
        tpfi$problem = tpdp$problem = tshap$problem = j
        tpfi$n = tpdp$n = tshap$n = ntrain
        tpfis = append(tpfis, list(tpfi))
        tpdps = append(tpdps, list(tpdp))
        tshaps = append(tshaps, list(tshap))
      }
    }
  }
  tpfis = rbindlist(tpfis)
  saveRDS(tpfis, file = tpfi_file)
  tpdps = rbindlist(tpdps)
  saveRDS(tpdps, file = tpdp_file)
  tshaps = rbindlist(tshaps)
  saveRDS(tshaps, file = tshap_file)
} else {
  tpfis = readRDS(tpfi_file)
  tpdps = readRDS(tpdp_file)
  tshaps = readRDS(tshap_file)
}

cis_pfi = merge(pfis, tpfis, by = c("feature", "algorithm", "problem", "n"))
cis_pdp = merge(pdps, tpdps, by = c("feature", "algorithm", "problem", "feature_value", "n"))
cis_shap = merge(shaps, tshaps, by = c("feature", "algorithm", "problem", "n"))

# =============================================================================
# Compute coverage for PFI
# =============================================================================
cis_pfi = cis_pfi[, in_ci := (lower <= tpfi) & (tpfi <= upper)]
coverage_pfi = cis_pfi[,.(coverage = mean(in_ci),
                  coverage_se = (1/N_EXPERIMENTS) * sd(in_ci),
                  avg_width = mean(upper - lower), 
                  bias = tpfi - pfi),
               by = list(feature, algorithm, problem, max_refits, n_perm, sampling_strategy, nrefits, 
                         adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method, m)]

coverage_pfi_mean = coverage_pfi[, .(coverage = mean(coverage), avg_width = mean(avg_width), coverage_se = mean(coverage_se), bias = mean(bias)),
                                     by = list(algorithm, problem, sampling_strategy, nrefits, 
                                               adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)]
print(coverage_pfi_mean)
saveRDS(coverage_pfi_mean, file = sprintf("%s/coverage_pfi_mean.Rds", res_dir))


# =============================================================================
# Compute coverage  for PDP
# =============================================================================
cis_pdp = cis_pdp[, in_ci := (lower <= tpdp) & (tpdp <= upper)]
coverage_pdp = cis_pdp[,.(coverage = mean(in_ci),
                  coverage_se = (1/N_EXPERIMENTS) * sd(in_ci),
                  avg_width = mean(upper - lower), 
                  bias = tpdp - pdp),
               by = list(feature, feature_value, algorithm, problem, sampling_strategy, max_refits, nrefits, 
                         adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method, m)]

coverage_pdp_mean = coverage_pdp[, .(coverage = mean(coverage), avg_width = mean(avg_width), coverage_se = mean(coverage_se), bias = mean(bias)),
                                 by = list(algorithm, problem, sampling_strategy, nrefits, 
                                           adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)]
saveRDS(coverage_pdp_mean, sprintf("%s/coverage_pdp_mean.Rds", res_dir))
print(coverage_pdp_mean)

# =============================================================================
# Compute coverage for SHAP
# =============================================================================
cis_shap = cis_shap[, in_ci := (lower <= tshap) & (tshap <= upper)]
coverage_shap = cis_shap[,.(coverage = mean(in_ci),
                          coverage_se = (1/N_EXPERIMENTS) * sd(in_ci),
                          avg_width = mean(upper - lower), 
                          bias = tshap - shap),
                       by = list(feature, algorithm, problem, max_refits, n_perm, sampling_strategy, nrefits, 
                                 adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method, m)]

coverage_shap_mean = coverage_shap[, .(coverage = mean(coverage), avg_width = mean(avg_width), coverage_se = mean(coverage_se), bias = mean(bias)),
                                 by = list(algorithm, problem, sampling_strategy, nrefits, 
                                           adjusted, n, missing_prob, pattern, train_missing, test_missing, imputation_method)]
print(coverage_shap_mean)
saveRDS(coverage_shap_mean, file = sprintf("%s/coverage_shap_mean.Rds", res_dir))


