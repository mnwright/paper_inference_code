# =============================================================================
# Functions used in CI experiments (code/ci-experiment.R)
# =============================================================================

# How big the training data should be (n1)
# At 0.632 subsampling and bootstrapping sizes coincide (in expectation)
SAMPLING_FRACTION = 0.632

# =============================================================================
# For data generating process
# =============================================================================

#' Data scenario x1234
#'
#' @param nsample Number of samples to generate
#' @return data.frame
generate_data1234 = function(nsample){
  x1 = runif(nsample, min = 0, max = 1)
  x2 = runif(nsample, min = 0, max = 1)
  x3 = runif(nsample, min = 0, max = 1)
  x4 = runif(nsample, min = 0, max = 1)
  y = x1 - sqrt(1 + x2) + x3 * x4 + (x4/10)^2 + rnorm(nsample, sd = 1)
  data.frame(x1, x2, x3, x4, y)
}

#' Wrapper for scenario x1234
#'
#' @param data data object from batchtools experiments
#' @param job batchtools job
#' @param n number of data points to sample
#' @return list with sampled data, scenario name and sample size
gdata = function(data, job, n, ...){
  dat = generate_data1234(n)
  list("dat" = dat, name = "x1234", n = n)
}


#' Data scenario x12
#'
#' @param nsample Number of samples to generate
#' @return data.frame
generate_data12 = function(nsample){
  x1 = runif(nsample, min = 0, max = 1)
  x2 = runif(nsample, min = 0, max = 1)
  y = x1 - x2 + rnorm(nsample, sd = 1)
  data.frame(x1, x2, y)
}

#' Wrapper for scenario x12
#'
#' @param data data object from batchtools experiments
#' @param job batchtools job
#' @param n number of data points to sample
#' @return list with sampled data, scenario name and sample size
gdata12 = function(data, job, n, ...){
 dat = generate_data12(n)
 list("dat" = dat, name = "x12", n = n)
}


# =============================================================================
# Lists of DGPs and models
# =============================================================================

rf100 = function(...) randomForest(..., ntree = 100)

dgps = list("x12" = generate_data12, "x1234" = generate_data1234)
mods = list("lm" = lm, "rpart" = rpart, "randomForest" = rf100)

# =============================================================================
# Model wrappers for batchtools
# =============================================================================

#' Funcion to produce Model Wrapper for batchtools
#' 
#' @param model Character name of model, see mods list
#' @return function, model wrapper
get_model_wrapper = function(model){
  # Retrieve model from list
  train_mod = mods[[model]]
  model_wrapper = function(data, job, instance, ...){
    gen_data = dgps[[instance$name]]
    #dat = instance$dat
    dat = gen_data(instance$n)
    nrefits = 1:job$prob.pars$max_refits
    samp_strategy = job$prob.pars$sampling_strategy
    gd = NA
    results = lapply(nrefits, function(m){
      if (samp_strategy == "subsampling"){
        train_ids = sample(1:nrow(dat), size = SAMPLING_FRACTION * nrow(dat), replace = FALSE)
        train_dat = dat[train_ids,]
        test_dat = dat[setdiff(1:nrow(dat), train_ids), ]
      } else if (samp_strategy == "bootstrap") {
        # Size here is n and not n1. Due to replace = TRUE, n_unique(train_dat) =~ 0.632 * n
        train_ids = sample(1:nrow(dat), size = nrow(dat), replace = TRUE)
        train_dat = dat[train_ids,]
        test_dat = dat[setdiff(1:nrow(dat), train_ids), ]
      } else if (samp_strategy == "ideal"){
        # Completely fresh data for both training and test
        train_dat = gen_data(SAMPLING_FRACTION * nrow(dat))
        test_dat  = gen_data(SAMPLING_FRACTION * nrow(dat))
        #gd = gen_data
      } else {
        print(sprintf("Strategy %s not implemented", job$prob.pars$sampling_strategy))
      }
      
      # Introduce missing data and impute
      pattern = job$prob.pars$pattern
      missing_prob = job$prob.pars$missing_prob
      if (pattern == "MCAR") {
        miss_fun <- function(data) missMethods::delete_MCAR(data, p = missing_prob)
      } else if (pattern == "MAR") {
        miss_fun <- function(data) {
          # Random half (rounded down) of columns are missing, the other half are used as control variables
          cols_mis <- sample(1:ncol(data), floor(ncol(data)/2))
          cols_ctrl <- sample(setdiff(1:ncol(data), cols_mis), floor(ncol(data)/2))
          missMethods::delete_MAR_rank(data, p = missing_prob, 
                                       cols_mis = cols_mis, cols_ctrl = cols_ctrl) 
        }
      } else if (pattern == "MNAR") {
        miss_fun <- function(data) {
          cols_mis <- 1:ncol(data)
          missMethods::delete_MNAR_rank(data, p = missing_prob, cols_mis = cols_mis)
        } 
      } else {
        stop("Unknown missing data pattern")
      }
      
      # Impute missing data
      imputation_method = job$prob.pars$imputation_method
      if (imputation_method == "mean") {
        impute_fun <- function(data) list(missMethods::impute_mean(data))
      } else if (imputation_method == "mice") {
        impute_fun <- function(data) {
          imp <- mice::mice(data, m = job$prob.pars$missing_prob * 100, print = FALSE)
          complete(imp, "all")
        }
      } else if (imputation_method == "missForest") {
        impute_fun <- function(data) list(missRanger(data, verbose = 0, num.threads = 1))
      } else {
        stop("Unknown imputation method")
      }
      
      train_missing = job$prob.pars$train_missing
      test_missing = job$prob.pars$test_missing
      if (train_missing) {
        train_dat <- miss_fun(train_dat)
        train_dat <- impute_fun(train_dat)
      } else {
        train_dat <- list(train_dat)
      }
      if (test_missing) {
        test_dat <- miss_fun(test_dat)
        test_dat <- impute_fun(test_dat)
      } else {
        test_dat <- list(test_dat)
      }
      
      # Creates the prediction function
      mods <- lapply(1:length(train_dat), function(i) {
        train_mod(y ~ ., data = train_dat[[i]])
      })
      
      imps <- max(length(train_dat), length(test_dat))
      
      # PFI
      pfis = rbindlist(lapply(1:imps, function(i) {
        if (length(train_dat) > 1) {
          fh = function(x) predict(mods[[i]], newdata = x)
        } else {
          fh = function(x) predict(mods[[1]], newdata = x)
        }
        if (length(test_dat) > 1) {
          td = test_dat[[i]]
        } else {
          td = test_dat[[1]]
        }
        
        pfis = compute_pfis(td, job$prob.pars$n_perm, fh, gen_data = gd)
        pfis$refit_id = m
        pfis$imp_id = i
        pfis
      }))
      
      # PDP
      pdps = rbindlist(lapply(1:imps, function(i) {
        if (length(train_dat) > 1) {
          fh = function(x) predict(mods[[i]], newdata = x)
        } else {
          fh = function(x) predict(mods[[1]], newdata = x)
        }
        if (length(test_dat) > 1) {
          td = test_dat[[i]]
        } else {
          td = test_dat[[1]]
        }
        
        pdps = compute_pdps(td, fh)
        pdps$refit_id = m
        pdps$imp_id = i
        pdps
      }))
      
      # SHAP
      shaps = rbindlist(lapply(1:imps, function(i) {
        if (length(train_dat) > 1) {
          fh = function(x) predict(mods[[i]], newdata = x)
          fit = mods[[i]]
          trdat = train_dat[[i]] 
        } else {
          fh = function(x) predict(mods[[1]], newdata = x)
          fit = mods[[1]]
          trdat = train_dat[[1]] 
        }
        if (length(test_dat) > 1) {
          tedat = test_dat[[i]]
        } else {
          tedat = test_dat[[1]]
        }
        
        shaps = compute_shaps(tedat, trdat, fh, fit)
        shaps$refit_id = m
        shaps$imp_id = i
        shaps
      }))
      
      list("pfis" = pfis, "pdps" = pdps, shaps = shaps)
    })
    pfis = rbindlist(lapply(results, function(x) x[["pfis"]]))
    pdps = rbindlist(lapply(results, function(x) x[["pdps"]]))
    shaps = rbindlist(lapply(results, function(x) x[["shaps"]]))

    # Loop over refit_id to simulate different number of models
    res_pfi = rbindlist(lapply(2:max(pfis$refit_id), function(i) {
      t_alpha = qt(1 - 0.05/2, df = i - 1)
      res = compute_pfi_cis(pfis[refit_id <= i, ], t_alpha, adjust = FALSE, type = samp_strategy)
      res$adjusted = FALSE
      if (samp_strategy != "ideal"){
        res_adjusted = compute_pfi_cis(pfis[refit_id <= i, ], t_alpha, adjust = TRUE, type = samp_strategy)
        res = rbind(data.table(res),
                    data.table(res_adjusted, adjusted = TRUE))
      }
      res[, nrefits := i]
      res
    }))

    # Loop over refit_id to simulate different number of models
    res_pdp = rbindlist(lapply(2:max(pdps$refit_id), function(i) {
      t_alpha = qt(1 - 0.05/2, df = i - 1)
      res = compute_pdp_cis(pdps[refit_id <= i, ], t_alpha, adjust = FALSE, type = samp_strategy)
      res$adjusted = FALSE
      if (samp_strategy != "ideal") {
        res_adjusted = compute_pdp_cis(pdps[refit_id <= i, ], t_alpha, adjust = TRUE, type = samp_strategy)
        res = rbind(data.table(res),
                    data.table(res_adjusted, adjusted = TRUE))
      }
      res[, nrefits := i]
      res
    }))
    
    # Loop over refit_id to simulate different number of models
    res_shap = rbindlist(lapply(2:max(shaps$refit_id), function(i) {
      t_alpha = qt(1 - 0.05/2, df = i - 1)
      res = compute_shap_cis(shaps[refit_id <= i, ], t_alpha, adjust = FALSE, type = samp_strategy)
      res$adjusted = FALSE
      if (samp_strategy != "ideal") {
        res_adjusted = compute_shap_cis(shaps[refit_id <= i, ], t_alpha, adjust = TRUE, type = samp_strategy)
        res = rbind(data.table(res),
                    data.table(res_adjusted, adjusted = TRUE))
      }
      res[, nrefits := i]
      res
    }))

    res_shap$job.id = res_pdp$job.id = res_pfi$job.id = job$id
    list("pdp" = res_pdp, "pfi" = res_pfi, "shap" = res_shap)
  }
}

# Generate the wrappers 
lm_wrapper = get_model_wrapper("lm")
rpart_wrapper = get_model_wrapper("rpart")
rf_wrapper = get_model_wrapper("randomForest")

# =============================================================================
# Adjustment term for both PD and PFI
# =============================================================================


#' Calculate the adjustment term based on Nadeau and Bengio
#' 
#' When SAMPLING_FRACTION is set to 0.632, same for both
#' @param type Either 'bootstrap' or 'subsampling'
  #' @return Variance adjustment term (n2/n1)
get_adjustment_term = function(type){
  if (type == "bootstrap") {
    fraction = 0.632
  } else if (type == "subsampling") {
    fraction = SAMPLING_FRACTION
  }  else {
    stop("not implemented")
  }
  # same as n2/n1 from Nadeau & Bengio paper
  (1 - fraction) / fraction
}


# =============================================================================
# PFI-specific functions
# =============================================================================
#' Permute data
#'
#' @param dat data.frame
#' @param fname Name of feature to permute
#' @param gen_data If NA, permute the data. If function, then used to generate
#'        new data to replace dat[fname]. Set NA for resampling, and use fun.
#'        for the infinite/ideal scenario.
#' @return dat, but with permuted dat[fname[]
permute = function(dat, fname, gen_data = NA){
  dat2 = dat
  if(is.na(gen_data)) {
    dat2[,fname] = sample(dat[,fname], size = nrow(dat), replace = FALSE)
  } else {
    dat2[,fname] = gen_data(nrow(dat))[,fname]
  }
  dat2
}

#' Compute PFI for 1 feature (L2 based)
#'
#' @param dat data.frame with unseen data
#' @param fh prediction function
#' @param nperm Number of permutations to be used
#' @param fname Feature for which to compute PFI
#' @param gen_data NA for shuffling, a dat gen. function for "ideal"
#' @return PFI for feature, averaged over permutations 
pfi = function(dat, fh, nperm, fname, gen_data = NA){
  R = mean((fh(dat) - dat$y)^2)
  Rt = lapply(1:nperm, function(k) {
    dat_p = permute(dat, fname, gen_data)
    mean((fh(dat_p) - dat$y)^2)
  })
  Rt = mean(unlist(Rt))
  data.table(pfi = Rt - R)
}

#' Compute PFI for all features (L2-based)
#'
#' @param test_dat data.frame with unseen data
#' @param fh prediction function
#' @param nperm Number of permutations
#' @param gen_data NA for permutation. set to function for specific data generation.
#' @return data.frame with PFIs
compute_pfis = function(test_dat, nperm, fh, gen_data = NA){
  fnames = setdiff(colnames(test_dat), "y")
  pfis = lapply(fnames, function(fname) {
    pfis_x = pfi(test_dat, fh, nperm, fname, gen_data = gen_data)
    pfis_x$feature = fname
    pfis_x
  })
  rbindlist(pfis)
}

#' Compute confidence intervals for PFI
#' 
#' @param resx data.frame with the PFI results
#' @param t_alpha 1-alpha/2 quantile of t-distribution
#' @param adjust TRUE if variance adjustment term by Nadeau/Bengio should be used
#' @param type Either "bootstrap" or "subsampling". Ignored if adjust is FALSE. 
#' @return data.frame with PFIs and their lower and upper CI boundaries.
compute_pfi_cis = function(resx, t_alpha, adjust = FALSE, type = NULL){
  nrefits = length(unique(resx$refit_id))
  if (resx[, any(imp_id) > 1]) {
    aa = resx[, .(var = var(pfi),
                  pfi = mean(pfi)),
              by = list(feature, imp_id)]
    myfun = function(...) {
      mice::pool.scalar(...)$t
    }
    resx = aa[, .(var3 = myfun(pfi, var), 
                  pfi = mean(pfi)), 
              by = list(feature)]
  } else {
    resx = resx[, .(var3 = var(pfi),
                    pfi = mean(pfi)),
                by = list(feature)]
  }
  
  m = (1/nrefits)
  if (adjust) m = m + get_adjustment_term(type)
  resx$se2 = m * resx$var3
  resx$lower = resx$pfi - t_alpha * sqrt(resx$se2)
  resx$upper = resx$pfi + t_alpha * sqrt(resx$se2)
  resx
}


#' Compute the expected learner-PFIs for scenario
#'
#' Used as groundtruth in the experiment. Repeatedly draws
#' new data, fits model and computes PFIs.
#' Results are averaged over these PFIs.
#'
#' @param ntrue Number of times new data is sampled
#' @param ntrain Size of sampled data
#' @param gen_data function for data generation
#' @param train_mod function to train model
#' @param nperm Number of permutations
#' @return data.frame with true PFIs
get_true_pfi = function(ntrue, ntrain, gen_data, train_mod, nperm = 5){
  true_pfis = mclapply(1:ntrue,  mc.cores = NC, function(m){
    train_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    mod = train_mod(y ~ ., data = train_dat)
    fh = function(x) predict(mod, newdata = x)
    test_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    pfis = compute_pfis(test_dat, nperm, fh, gen_data = gen_data)
    pfis[,.(pfi = mean(pfi)), by = list(feature)]
  })
  true_pfis = rbindlist(true_pfis)
  tpfi = true_pfis[,.(tpfi = mean(pfi), tse = sqrt(1/ntrue) * sd(pfi)), by = list(feature)]
}

# =============================================================================
# PDP-specific functions
# =============================================================================

#' Compute PDP for 1 feature
#'
#' @param data.frame with data for MC integration
#' @param fh prediction function
#' @param fname Feature for which to computed PDP
#' @param xgrid Grid values at which to computed the PDP
#' @return data.frame with PDP values
pdp = function(dat, fh, fname, xgrid = c(0.1, 0.3, 0.5, 0.7, 0.9)){
  dat2 = dat
  res = lapply(xgrid, function(x_g){
    dat2[,fname] = x_g
    data.frame(pdp = mean(fh(dat2)),
               feature_value = x_g)
  })
  rbindlist(res)
}


#' Compute PDPs for all features
#' 
#' @param test_dat data.frame for MC integration
#' @param fh prediction function
#' @return data.frame with PDPs
compute_pdps = function(test_dat, fh){
  fnames = setdiff(colnames(test_dat), "y")
  pdps = lapply(fnames, function(fname) {
    pdp_dat = pdp(test_dat, fh, fname)
    pdp_dat$feature = fname
    pdp_dat
  })
  rbindlist(pdps)
}

#' Compute confidence intervals for PDP
#' 
#' @param resx data.frame with the PDP results
#' @param t_alpha 1-alpha/2 quantile of t-distribution
#' @param adjust TRUE if variance adjustment term by Nadeau/Bengio should be used
#' @param type Either "bootstrap" or "subsampling". Ignored if adjust is FALSE. 
#' @return data.frame with PDPs and their lower and upper CI boundaries.
compute_pdp_cis = function(resx, t_alpha, adjust = FALSE, type = NULL){
  nrefits = length(unique(resx$refit_id))
  if (resx[, any(imp_id) > 1]) {
    aa = resx[, .(var = var(pdp),
                  pdp = mean(pdp)),
              by = list(feature, feature_value, imp_id)]
    myfun = function(...) {
      mice::pool.scalar(...)$t
    }
    resx = aa[, .(var2 = myfun(pdp, var), 
                  mpdp = mean(pdp)), 
              by = list(feature, feature_value)]
  } else {
    resx = resx[, .(var2 = var(pdp),
                    mpdp =  mean(pdp)),
                by = list(feature, feature_value)]
  }
  
  m = (1/nrefits)
  if (adjust) m = m + get_adjustment_term(type)
  resx$se2 = m * resx$var2
  resx$lower = resx$mpdp - t_alpha * sqrt(resx$se2)
  resx$upper = resx$mpdp + t_alpha * sqrt(resx$se2)
  resx
}


#' Compute the expected learner-PDPs for scenario
#'
#' Used as groundtruth in the experiment. Repeatedly draws
#' new data, fits model and computes PDPs.
#' Results are averaged over these PDPs.
#'
#' @param ntrue Number of times new data is sampled
#' @param ntrain Size of sampled data
#' @param gen_data function for data generation
#' @param train_mod function to train model
#' @return data.frame with true PDPs.
get_true_pdp = function(ntrue, ntrain, gen_data, train_mod){
  true_pdps = mclapply(1:ntrue,  mc.cores = NC, function(m){
    # Make sure nsample is the same as for resampled models
    train_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    mod = train_mod(y ~ ., train_dat)
    fh = function(x) predict(mod, newdata = x)
    # Here sample size does not matter. Only tradeoff: Accuracy and computation time
    test_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    pdps = compute_pdps(test_dat, fh)
    pdps[,.(pdp = mean(pdp)), by = list(feature, feature_value)]
})
  true_pdps = rbindlist(true_pdps)
  true_pdps[,.(tpdp = mean(pdp), tse = sqrt(1/ntrue) * sd(pdp)), by = list(feature, feature_value)]
}

# =============================================================================
# SHAP-specific functions
# =============================================================================

#' Compute SHAP for all features 
#'
#' @param test_dat data.frame with unseen data
#' @param fh prediction function
#' @return data.frame with SHAP values
compute_shaps = function(test_dat, train_dat, fh, fit){
  fnames = setdiff(colnames(test_dat), "y")
  explanation  <- shapr::explain(fit, test_dat[, fnames], train_dat[, fnames], 
                                 approach = "empirical", phi0 = mean(train_dat$y), 
                                 predict_model = stats::predict, verbose = NULL)
  shaps <- colMeans(abs(explanation$shapley_values_est[, ..fnames]))
  data.table(shap = shaps, feature = names(shaps))
  #browser()
  #melt(explanation$shapley_values_est[, ..fnames], measure.vars = fnames, 
  #     variable.name = "feature", value.name = "shap")
}

#' Compute confidence intervals for SHAP
#' 
#' @param resx data.frame with the SHAP results
#' @param t_alpha 1-alpha/2 quantile of t-distribution
#' @param adjust TRUE if variance adjustment term by Nadeau/Bengio should be used
#' @param type Either "bootstrap" or "subsampling". Ignored if adjust is FALSE. 
#' @return data.frame with SHAP values and their lower and upper CI boundaries.
compute_shap_cis = function(resx, t_alpha, adjust = FALSE, type = NULL){
  nrefits = length(unique(resx$refit_id))
  if (resx[, any(imp_id) > 1]) {
    aa = resx[, .(var = var(shap),
                  shap = mean(shap)),
              by = list(feature, imp_id)]
    myfun = function(...) {
      mice::pool.scalar(...)$t
    }
    resx = aa[, .(var3 = myfun(shap, var), 
                  shap = mean(shap)), 
              by = list(feature)]
  } else {
    resx = resx[, .(var3 = var(shap),
                    shap = mean(shap)),
                by = list(feature)]
  }
  
  m = (1/nrefits)
  if (adjust) m = m + get_adjustment_term(type)
  resx$se2 = m * resx$var3
  resx$lower = resx$shap - t_alpha * sqrt(resx$se2)
  resx$upper = resx$shap + t_alpha * sqrt(resx$se2)
  resx
}

#' Compute the expected learner SHAP values for scenario
#'
#' Used as groundtruth in the experiment. Repeatedly draws
#' new data, fits model and computes SHAP values
#' Results are averaged over these SHAP values.
#'
#' @param ntrue Number of times new data is sampled
#' @param ntrain Size of sampled data
#' @param gen_data function for data generation
#' @param train_mod function to train model
#' @return data.frame with true SHAP values
get_true_shap = function(ntrue, ntrain, gen_data, train_mod){
  true_shapss = mclapply(1:ntrue,  mc.cores = NC, function(m){
    # Make sure nsample is the same as for resampled models
    train_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    mod = train_mod(y ~ ., train_dat)
    fh = function(x) predict(mod, newdata = x)
    test_dat = gen_data(nsample = SAMPLING_FRACTION * ntrain)
    shaps = compute_shaps(test_dat, train_dat, fh, mod)
    shaps[,.(shaps = mean(shaps)), by = list(feature)]
  })
  true_pdps = rbindlist(true_pdps)
  true_pdps[,.(tshap = mean(shaps), tse = sqrt(1/ntrue) * sd(shaps)), by = list(feature)]
}
