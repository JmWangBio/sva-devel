#' Adjust for batch effects using a beta-binomial regression framework in bisulfite sequencing count data
#' 
#' ComBat_biseq is a model adapted from ComBat-seq using beta-binomial regression,
#' which specifically targets bisulfite sequencing count data.  
#'
#' @inheritParams ComBat_met
#' @param numCs matrix of number of methylated cytosines
#' @param coverage matrix of coverage
#' @param shrink Boolean variable indicating whether to apply EB-shrinkage on parameter estimation
#' 
#' @return \code{ComBat_biseq} returns a batch-adjusted feature x sample count matrix, representing number
#' of methylated cytosines
#' @export
#' 
#' @examples
#' library(methylKit)
#' # Simulate bisulfite sequencing data
#' my.methylBase <- dataSim(replicates = 8, sites = 50, treatment = c(1,1,1,1,0,0,0,0),
#'                          percentage = 10, effect = 25, add.info = TRUE)
#' my.methylData <- getData(my.methylBase[[1]])
#' coverage_mat <- as.matrix(my.methylData[, grep("coverage", colnames(my.methylData))])
#' class(coverage_mat) <- "numeric"
#' numCs_mat <- as.matrix(my.methylData[, grep("numCs", colnames(my.methylData))])
#' class(numCs_mat) <- "numeric"
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#' 
#' # Adjust for batch effects including biological conditions
#' adj_numCs_mat <- ComBat_biseq(numCs = numCs_mat, coverage = coverage_mat, batch = batch, 
#' group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_numCs_mat <- ComBat_biseq(numCs = numCs_mat, coverage = coverage_mat, batch = batch, 
#' group = group, full_mod = FALSE)
#' 

ComBat_biseq <- function(numCs, coverage, batch, group = NULL, covar_mod = NULL, full_mod = TRUE, 
                         shrink = FALSE, mean.only = FALSE, feature.subset.n = NULL) {  
  ########  Preparation  ########
  ## Check if coverage and numCs are the same size
  if (!all(dim(coverage) == dim(numCs))) {
    stop("Coverage matrix and numCs matrix are not the same size.")
  }
  
  ## check if coverage and batch have matching sizes
  if (ncol(coverage) != length(batch)) {
    stop("Coverage matrix and batch vector do not have matching sizes.")
  }
  
  ## check if coverage and group have matching sizes
  if (!is.null(group)) {
    if (ncol(coverage) != length(group)) {
      stop("Coverage matrix and group vector do not have matching sizes.")
    }
  }
  
  ## Does not allow numCs to be larger than coverage
  if (sum(numCs > coverage, na.rm = TRUE) > 0) {
    stop("Number of methylated Cs cannot exceed total number of Cs.")
  }
  
  ## Does not support more than 1 batch variable
  if (length(dim(batch)) > 1) {
    stop("ComBat-biseq does not allow more than one batch variable!")
  }
  
  ## Does not support 1 sample across all batches
  batch <- as.factor(batch)
  if (all(table(batch) <= 1)) {
    stop("ComBat-biseq doesn't support only 1 sample across all batches!")
  }
  
  ## Does not support only 1 batch level
  if (length(levels(batch)) <= 1) {
    stop("Found only one batch. No need to adjust for batch effects!")
  }
  
  ## Correct for mean batch effects only if any batch has only 1 sample
  if (any(table(batch) == 1)) {
    cat("At least one batch contains only 1 sample. Only mean batch effects will be corrected.\n")
    mean.only <- TRUE
  }
  
  ## Remove features with zero variance across all batches
  zero.var.rows.lst <- lapply(levels(batch), function(b) {
    intersect(which(apply(coverage[, batch == b], 1, function(x){var(x, na.rm = TRUE) == 0})), 
              which(apply(numCs[, batch == b], 1, function(x){var(x, na.rm = TRUE) == 0})))
  })
  all.zero.var.rows <- Reduce(intersect, zero.var.rows.lst)
  
  if (length(all.zero.var.rows) > 0) {
    cat(sprintf("Found %s features with uniform numCs and coverage values across all batches; 
                these features will not be adjusted for batch effects.\n", 
                length(all.zero.var.rows)))
  }
  
  keep <- setdiff(1:nrow(coverage), all.zero.var.rows)
  numCsOri <- numCs
  coverageOri <- coverage
  numCs <- numCsOri[keep, ]
  coverage <- coverageOri[keep, ]
  
  ## Create a vector for correction types
  if (mean.only) {
    mean.only.vec <- rep(TRUE, length(keep))
  } else {
    mean.only.vec <- rep(FALSE, length(keep))
  }
  
  ## Calculate batch-related statistics
  # number of batches
  n_batch <- nlevels(batch)
  # list of samples in each batch  
  batches_ind <- lapply(1:n_batch, function(i) {
    which(batch == levels(batch)[i])
  })
  # number of samples in each batch
  n_batches <- sapply(batches_ind, length)
  # total number of samples
  n_sample <- sum(n_batches)  
  cat(sprintf("Found %s batches and %s samples.\n", n_batch, n_sample))
  
  ## Make design matrix
  # biological condition matrix
  group <- as.factor(group)
  if (full_mod & nlevels(group) > 1) {
    cat("Using full model in ComBat-biseq.\n")
    mod <- model.matrix(~group)
  } else {
    cat("Using null model in ComBat-biseq.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(coverage)))
  }
  # covariate matrix
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i) {
        model.matrix(~covar_mod[,i])
      }))
    }
  }
  # combine covariate matrix with biological condition matrix
  mod <- cbind(mod, covar_mod)
  # combine with batch matrix
  batchmod <- model.matrix(~-1 + batch)
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), 'covariate(s) or covariate level(s).\n')
  
  ## Check if the design is confounded
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-biseq.")
    }
    if (ncol(design) > (n_batch+1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, -c(1:n_batch)]))) {
        stop('The covariates are confounded! Please remove one or more of the covariates 
             so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded 
             covariates and rerun ComBat-biseq.")
      }
    }
  }
  
  ## Check for missing values in count matrices
  NAs = any(is.na(coverage) | is.na(numCs))
  if (NAs) {
    cat(c('Found', sum(is.na(coverage)), 
          'missing data values in the coverage matrix.\n'), sep = ' ')
    cat(c('Found', sum(is.na(numCs)),
          'missing data values in the numCs matrix.\n'), sep = ' ')
  }
  
  ########  Estimate parameters from beta-binomial GLM  ########
  cat("Fitting the GLM model\n")
  gamma_hat_lst <- vector("list", length = nrow(coverage))
  mu_hat_lst <- vector("list", length = nrow(coverage))
  phi_hat_lst <- vector("list", length = nrow(coverage))
  delta_hat_lst <- vector("list", length = nrow(coverage))
  n_zero_modvar <- 0
  
  for (k in 1:nrow(coverage)) {
    if (k %% 1000 == 0) {
      cat(sprintf("%s features processed\n", k))
    }
    # mark rows with NA values and zero coverage
    full_mat <- as.data.frame(cbind(design, 
                                    coverage = coverage[k, ], 
                                    numCs = numCs[k, ]))
    full_mat$coverage[full_mat$coverage == 0] <- NA
    full_mat$numCs[full_mat$coverage == 0] <- NA
    nona <- which(stats::complete.cases(full_mat))
    
    # check if the data are all NAs
    if (length(nona) == 0) {
      n_zero_modvar <- n_zero_modvar + 1
      next
    }
    
    # check if the model has zero model variance
    full_ratio <- full_mat$numCs / full_mat$coverage
    if (qr(cbind(full_mat[, 1:(ncol(full_mat) - 2)], full_ratio)[nona, ])$rank < ncol(full_mat) - 1)  {
      n_zero_modvar <- n_zero_modvar + 1
      next
    }
    
    # model fit
    full_df <- as.data.frame(full_mat)
    
    if (mean.only.vec[k]) {
      glm_f <- aod::betabin(cbind(numCs, coverage - numCs) ~ -1 + ., 
                            ~ 1, 
                            data = full_df[nona, ])
    } else {
      batch_nona <- batch[nona]
      glm_f <- aod::betabin(cbind(numCs, coverage - numCs) ~ -1 + ., 
                            ~ batch_nona, 
                            data = full_df[nona, ])
    }
    
    # compute mean and precision intercepts as batch-size-weighted average from batches
    alpha_x <- as.numeric(aod::coef(glm_f)[1:n_batch] %*% 
                            as.matrix(colSums(batchmod[nona, ]) / length(nona)))
    if (mean.only.vec[k]) {
      alpha_z <- as.numeric(log(1 / glm_f@random.param - 1))
    } else {
      alpha_z <- as.numeric(log(1 / glm_f@random.param - 1) %*%
                              as.matrix(colSums(batchmod[nona, ]) / length(nona)))
    }
    
    # estimate parameters
    gamma_hat <- aod::coef(glm_f)[1:n_batch] - alpha_x
    mu_hat <- rep(NA, nrow(full_mat))
    mu_hat[nona] <- aod::fitted(glm_f)
    phi_hat <- as.numeric(exp(alpha_z)) * rep(1, nrow(full_mat))
    if (mean.only.vec[k]) {
      delta_hat <- rep(0, n_batch)
    } else {
      delta_hat <- log(1 / glm_f@random.param - 1) - alpha_z
    }
    
    # store result
    gamma_hat_lst[[k]] <- gamma_hat
    mu_hat_lst[[k]] <- mu_hat
    phi_hat_lst[[k]] <- phi_hat
    delta_hat_lst[[k]] <- delta_hat
  }
  
  cat(sprintf("Found %s features with zero model variance;
              these features won't be adjusted for batch effects.\n",
              n_zero_modvar))
  
  # convert NULLs to NAs
  gamma_hat_lst[sapply(gamma_hat_lst, is.null)] <- NA
  mu_hat_lst[sapply(mu_hat_lst, is.null)] <- NA
  phi_hat_lst[sapply(phi_hat_lst, is.null)] <- NA
  delta_hat_lst[sapply(delta_hat_lst, is.null)] <- NA
  
  # reformat lists as matrices
  gamma_hat_mat <- do.call('rbind', gamma_hat_lst)
  mu_hat_mat <- do.call('rbind', mu_hat_lst)
  phi_hat_mat <- do.call('rbind', phi_hat_lst)
  delta_hat_mat <- do.call('rbind', delta_hat_lst)
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########
  if (shrink) {
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_betabin
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(numCs.dat = numCs[, batches_ind[[ii]]],
                           coverage.dat = coverage[, batches_ind[[ii]]],
                           mu = mu_hat_mat[, batches_ind[[ii]]],
                           gamma = gamma_hat_mat[, ii],
                           phi = phi_hat_mat[, batches_ind[[ii]]],
                           delta = delta_hat_mat[, ii],
                           feature.subset.n = feature.subset.n)
      } else {
        invisible(capture.output(mcres <- mcint_fun(numCs.dat = numCs[, batches_ind[[ii]]],
                                                    coverage.dat = coverage[, batches_ind[[ii]]],
                                                    mu = mu_hat_mat[, batches_ind[[ii]]],
                                                    gamma = gamma_hat_mat[, ii],
                                                    phi = phi_hat_mat[, batches_ind[[ii]]],
                                                    delta = delta_hat_mat[, ii],
                                                    feature.subset.n = feature.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res) {res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    delta_star_mat <- lapply(monte_carlo_res, function(res) {res$delta_star})
    delta_star_mat <- do.call(cbind, delta_star_mat)
    
    if (mean.only) {
      cat("Apply shrinkage to mean only\n")
      delta_star_mat <- delta_hat_mat
    }
  } else {
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat_mat
    delta_star_mat <- delta_hat_mat
  }
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star_mat <- matrix(NA, nrow = nrow(coverage), ncol = ncol(coverage))
  phi_star_mat <- phi_hat_mat
  
  for (jj in 1:n_batch) {
    logit_mu_star_subset <- log(mu_hat_mat[, batches_ind[[jj]]] / 
                                  (1 - mu_hat_mat[, batches_ind[[jj]]])) - 
      vec2mat(gamma_star_mat[, jj], n_batches[jj])
    mu_star_mat[, batches_ind[[jj]]] <- exp(logit_mu_star_subset) / 
      (1 + exp(logit_mu_star_subset))
    if (!mean.only) {
      log_phi_star_subset <- log(phi_hat_mat[, batches_ind[[jj]]]) +
        vec2mat(delta_hat_mat[, jj], n_batches[jj]) -
        vec2mat(delta_star_mat[, jj], n_batches[jj])
      phi_star_mat[, batches_ind[[jj]]] <- exp(log_phi_star_subset)
    }
  }
  
  ########  Adjust the data  ########
  cat("Adjusting the data\n")
  adj_numCs_raw <- matrix(NA, nrow = nrow(coverage), ncol = ncol(coverage))
  for (kk in 1:n_batch) {
    numCs_sub <- numCs[, batches_ind[[kk]]]
    coverage_sub <- coverage[, batches_ind[[kk]]]
    old_mu <- mu_hat_mat[, batches_ind[[kk]]]
    old_phi <- phi_hat_mat[, batches_ind[[kk]]] * exp(delta_hat_mat[, kk])
    new_mu <- mu_star_mat[, batches_ind[[kk]]]
    new_phi <- phi_star_mat[, batches_ind[[kk]]]
    adj_numCs_raw[, batches_ind[[kk]]] <- match_quantiles_betabin(numCs_sub = numCs_sub,
                                                                  coverage_sub = coverage_sub,
                                                                  old_mu = old_mu, 
                                                                  old_phi = old_phi, 
                                                                  new_mu = new_mu, 
                                                                  new_phi = new_phi)
  }
  
  # Negative counts might occur for numeric reasons; convert negative counts to zero
  adj_numCs_raw[adj_numCs_raw < 0] <- 0
  
  ## Add back features with uniform coverage and numCs values so that dimensions won't change
  adj_numCs <- numCsOri
  adj_numCs[keep, ] <- adj_numCs_raw
  return(adj_numCs)
}
