# =============================================================================
# This script implements the CI coverage experiment
# =============================================================================


# =============================================================================
# Experiment Settings
# =============================================================================
N_EXPERIMENTS = 100 #10000
N_TRAIN = c(100, 1000)
#N_SAMPLE = 100
MAX_REFITS = 30
# Number of permutations for PFI computation
N_PERM = 5
# Number of refits to estimate true PFI/PDP
# This should be set high
N_TRUE = 10000
NC = 10

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
#addAlgorithm(name = "randomForest",  rf_wrapper)

strgs = "ideal" #c("subsampling", "bootstrap", "ideal")
impute = c("none", "mean", "mice")
setting = expand.grid(n = N_TRAIN,
                      max_refits = MAX_REFITS,
                      n_perm = N_PERM,
                      sampling_strategy = strgs, 
                      impute = impute)

pdes = list(x12 = setting, x1234 = setting)

addExperiments(pdes, repls = N_EXPERIMENTS)
summarizeExperiments()


# =============================================================================
# Run Experiment
# =============================================================================

#testJob(10)
if (grepl("node\\d{2}|bipscluster", system("hostname", intern = TRUE))) {
  ids <- findNotStarted()
  ids[, chunk := chunk(job.id, chunk.size = 200)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = "iml", chunks.as.arrayjobs = TRUE,
  				              ncpus = 1, memory = 6000, walltime = 10*24*3600,
  							        max.concurrent.jobs = 200))
} else if (grepl("glogin\\d+", system("hostname", intern = TRUE))) {
  ids <- findNotSubmitted()#[1:400, ]
  ids[, chunk := chunk(job.id, chunk.size = 800)]
  submitJobs(ids = ids, # walltime in seconds, 10 days max, memory in MB
             resources = list(name = "iml", chunks.as.arrayjobs = FALSE, 
                              partition = "medium40:test", ntasks = 40, # ntaskspernode depends on partition
                              max.concurrent.jobs = 8000,
                              #partition = "large40:shared", 
  				              ncpus = 1, walltime = 3600))
} else {
  reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = NC, fs.latency = 0)
  submitJobs()
}
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

# =============================================================================
# Compute true importance and PDP
# =============================================================================

tpfis = list()
tpdps = list()

tpfi_file = sprintf("%s/tpfis.Rds", res_dir)
tpdp_file = sprintf("%s/tpdps.Rds", res_dir)

if (!file.exists(tpfi_file)){
  mod_names = c("lm", "rpart", "randomForest")
  dgp_names = c("x12", "x1234")
  for(i in mod_names) {
    for (j in dgp_names) {
      for (ntrain in N_TRAIN) {
        message(i, j, ntrain)
        tpfi = get_true_pfi(N_TRUE, ntrain = ntrain, dgps[[j]], mods[[i]])
        tpdp = get_true_pdp(N_TRUE, ntrain = ntrain, dgps[[j]], mods[[i]])
        tpfi$algorithm = tpdp$algorithm = i
        tpfi$problem = tpdp$problem = j
        tpfi$n = tpdp$n = ntrain
        tpfis = append(tpfis, list(tpfi))
        tpdps = append(tpdps, list(tpdp))
      }
    }
  }
  tpfis = rbindlist(tpfis)
  saveRDS(tpfis, file = tpfi_file)
  tpdps = rbindlist(tpdps)
  saveRDS(tpdps, file = tpdp_file)
} else {
  tpfis = readRDS(tpfi_file)
  tpdps = readRDS(tpdp_file)
}

cis_pfi = merge(pfis, tpfis, by = c("feature", "algorithm", "problem", "n"))
cis_pdp = merge(pdps, tpdps, by = c("feature", "algorithm", "problem", "feature_value", "n"))

# =============================================================================
# Compute coverage for PFI
# =============================================================================
cis_pfi = cis_pfi[, in_ci := (lower <= tpfi) & (tpfi <= upper)]
coverage_pfi = cis_pfi[,.(coverage = mean(in_ci),
                  coverage_se = (1/N_EXPERIMENTS) * sd(in_ci),
                  avg_width = mean(upper - lower)),
               by = list(feature, algorithm, problem, max_refits, n_perm, sampling_strategy, nrefits, adjusted, n, impute)]

coverage_pfi_mean = coverage_pfi[, .(coverage = mean(coverage), avg_width = mean(avg_width), coverage_se = mean(coverage_se)),
                                     by = list(algorithm, problem, sampling_strategy, nrefits, adjusted, n, impute)]
print(coverage_pfi_mean)
saveRDS(coverage_pfi_mean, file = sprintf("%s/coverage_pfi_mean.Rds", res_dir))


# =============================================================================
# Compute coverage  for PDP
# =============================================================================
cis_pdp = cis_pdp[, in_ci := (lower <= tpdp) & (tpdp <= upper)]
coverage_pdp = cis_pdp[,.(coverage = mean(in_ci),
                  coverage_se = (1/N_EXPERIMENTS) * sd(in_ci),
                  avg_width = mean(upper - lower)),
               by = list(feature, feature_value, algorithm, problem, sampling_strategy, max_refits, nrefits, adjusted, n, impute)]

coverage_pdp_mean = coverage_pdp[, .(coverage = mean(coverage), avg_width = mean(avg_width), coverage_se = mean(coverage_se)),
                                 by = list(algorithm, problem, sampling_strategy, nrefits, adjusted, n, impute)]
saveRDS(coverage_pdp_mean, sprintf("%s/coverage_pdp_mean.Rds", res_dir))
print(coverage_pdp_mean)


