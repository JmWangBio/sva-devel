#' Adjust for batch effects using a beta regression framework in DNA methylation data
#'
#' ComBat_met is a model adapted from ComBat and ComBat-seq using beta regression,
#' which specifically targets DNA methylation data.
#'
#' @param bv matrix of beta-values
#' @param batch vector for batch
#' @param group optional vector for biological condition of interest
#' @param covar_mod optional model matrix representing co-variates to be included in the model
#' @param full_mod Boolean variable indicating whether to include biological condition of interest in the model
#' @param shrink Boolean variable indicating whether to apply EB-shrinkage on parameter estimation
#' @param mean.only Boolean variable indicating whether to apply EB-shrinkage on the estimation of dispersion effects
#' @param feature.subset.n number of features to use in non-parametric EB estimation, only useful when shrink equals TRUE
#' @param pseudo_beta pseudo beta-values to be used for replacing extreme 0 and 1 beta-values. 
#' Value needs to be between 0 and 0.5.
#' @param ref.batch NULL by default. If given, that batch will be selected as a reference for batch correction.
#' 
#' @return \code{ComBat_met} returns a feature x sample beta-value matrix, adjusted for batch effects.
#' @export
#' 
#' @examples
#' # Generate a random beta-value matrix
#' bv_mat <- matrix(runif(n = 400, min = 0, max = 1), nrow = 50, ncol = 8)
#' batch <- c(rep(1, 4), rep(2, 4))
#' group <- rep(c(0, 1), 4)
#'
#' # Adjust for batch effects including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, batch = batch, group = group, full_mod = TRUE)
#' # Adjust for batch effects without including biological conditions
#' adj_bv_mat <- ComBat_met(bv_mat, batch = batch, group = group, full_mod = FALSE)
#'

ComBat_met <- function(bv, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                       shrink = FALSE, mean.only = FALSE, feature.subset.n = NULL,
                       pseudo_beta = 1e-4, ref.batch = NULL) {
  ########  Preparation  ########
  ## convert extreme 0 or 1 values to pseudo-beta
  if (pseudo_beta <= 0 | pseudo_beta >= 0.5) {
    stop("Invalid pseudo beta-values.")
  }
  bv[bv == 0] <- pseudo_beta
  bv[bv == 1] <- 1 - pseudo_beta
  
  ## check if beta values and batch have matching sizes
  if (ncol(bv) != length(batch)) {
    stop("Coverage matrix and batch vector do not have matching sizes.")
  }
  
  ## check if beta values and group have matching sizes
  if (!is.null(group)) {
    if (ncol(bv) != length(group)) {
      stop("Coverage matrix and group vector do not have matching sizes.")
    }
  }
  
  ## Does not allow beta values outside the (0, 1) range
  if (sum(bv >= 1, na.rm = TRUE) > 0 | sum(bv <= 0, na.rm = TRUE) > 0) {
    stop("All beta values must be between 0 and 1.")
  }
  
  ## Does not support more than 1 batch variable
  if (length(dim(batch)) > 1) {
    stop("ComBat-met does not allow more than one batch variable!")
  }
  
  ## Does not support 1 sample across all batches
  batch <- as.factor(batch)
  if (all(table(batch) <= 1)) {
    stop("ComBat-met doesn't support only 1 sample across all batches!")
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
    which(apply(bv[, batch == b], 1, function(x) {var(x, na.rm = TRUE) == 0}))
  })
  all.zero.var.rows <- Reduce(intersect, zero.var.rows.lst)
  
  if (length(all.zero.var.rows) > 0) {
    cat(sprintf("Found %s features with uniform beta values across all batches; 
                these features will not be adjusted for batch effects.\n",
                length(all.zero.var.rows)))
  }
  
  keep <- setdiff(1:nrow(bv), all.zero.var.rows)
  bvOri <- bv
  bv <- bvOri[keep, ]
  
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
    cat("Using full model in ComBat-met.\n")
    mod <- model.matrix(~group)
  } else {
    cat("Using null model in ComBat-met.\n")
    mod <- model.matrix(~1, data = as.data.frame(t(bv)))
  }
  # covariate matrix
  if (!is.null(covar_mod)) {
    if (is.data.frame(covar_mod)) {
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i) {
        model.matrix(~covar_mod[, i])
      }))
    }
  }
  # combine covariate matrix with biological condition matrix
  mod <- cbind(mod, covar_mod)
  # combine with batch matrix
  batchmod <- model.matrix(~-1 + batch)
  if (!is.null(ref.batch)) {
    ## check for reference batch and make appropriate changes
    if (!(ref.batch %in% levels(batch))) {
      stop("Reference level ref. batch is not one of the levels of the batch variable.")
    }
    cat(sprintf("Using batch %s as the reference batch\n", ref.batch))
    ref <- which(levels(batch) == ref.batch)
  } else {
    ref <- NULL
  }
  
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  cat("Adjusting for", ncol(design) - ncol(batchmod), 'covariate(s) or covariate level(s).\n')
  
  ## Check if the design is confounded
  if (qr(design)$rank < ncol(design)) {
    if (ncol(design) == (n_batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-met.")
    }
    if (ncol(design) > (n_batch+1)) {
      if ((qr(design[, -c(1:n_batch)])$rank < ncol(design[, -c(1:n_batch)]))) {
        stop('The covariates are confounded! Please remove one or more of the covariates 
             so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded 
             covariates and rerun ComBat-met.")
      }
    }
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(bv))
  if (NAs) {cat(c('Found', sum(is.na(bv)), 'missing data values\n'), sep = ' ')}
  
  ########  Estimate parameters from beta GLM  ########
  cat("Fitting the GLM model\n")
  gamma_hat_lst <- vector("list", length = nrow(bv))
  mu_hat_lst <- vector("list", length = nrow(bv))
  phi_hat_lst <- vector("list", length = nrow(bv))
  delta_hat_lst <- vector("list", length = nrow(bv))
  n_zero_modvar <- 0
  n_zero_modvar_batch <- 0
  n_moderr <- 0
  
  for (k in 1:nrow(bv)) {
    if (k %% 1000 == 0) {
      cat(sprintf("%s features processed\n", k))
    }
    # mark rows with NA values
    full_mat <- cbind(design, bv[k, ])
    nona <- which(stats::complete.cases(full_mat))
    
    # check if the data are all NAs
    if (length(nona) == 0) {
      n_zero_modvar <- n_zero_modvar + 1
      next
    }
    
    # check if the model has zero model variance
    if (qr(full_mat[nona, ])$rank < ncol(full_mat)) {
      n_zero_modvar <- n_zero_modvar + 1
      next
    }
    
    # if dispersion correction enabled, check whether the model has zero model 
    # variance within any batch
    if (!mean.only.vec[k]) {
      for (i in 1:length(batches_ind)) {
        if (qr(full_mat[intersect(batches_ind[[i]], nona), c(i, (n_batch+1):ncol(full_mat))])$rank < 
            ncol(full_mat) - n_batch + 1) {
          n_zero_modvar_batch <- n_zero_modvar_batch + 1
          next
        }
      }
    }
    
    # model fit
    if (mean.only.vec[k]) {
      glm_f <- tryCatch({
        betareg::betareg.fit(x = design[nona, ], y = bv[k, ][nona])
      }, error = function(e) {
        e
      })
    } else {
      glm_f <- tryCatch({
        betareg::betareg.fit(x = design[nona, ], y = bv[k, ][nona], 
                             z = batchmod[nona, ])
      }, error = function(e) {
        e
      })
    }
    
    # if error with model fitting
    if (inherits(glm_f, "error")) {
      n_moderr <- n_moderr + 1
      next
    }
    
    # compute mean and precision intercepts as batch-size-weighted average from batches
    if (!is.null(ref.batch)) {
      alpha_x <- glm_f$coefficients$mean[ref]
    } else {
      alpha_x <- glm_f$coefficients$mean[1:n_batch] %*% 
        as.matrix(colSums(batchmod[nona, ]) / length(nona))
    }
    
    if (mean.only.vec[k]) {
      alpha_z <- glm_f$coefficients$precision
    } else {
      if (!is.null(ref.batch)) {
        alpha_z <- glm_f$coefficients$precision[ref]
      } else {
        alpha_z <- glm_f$coefficients$precision %*%
          as.matrix(colSums(batchmod[nona, ]) / length(nona))
      }
    }
    
    # estimate parameters
    gamma_hat <- glm_f$coefficients$mean[1:n_batch] - as.numeric(alpha_x)
    mu_hat <- rep(NA, nrow(full_mat))
    mu_hat[nona] <- glm_f$fitted.values
    phi_hat <- as.numeric(exp(alpha_z)) * rep(1, nrow(full_mat))
    if (mean.only.vec[k]) {
      delta_hat <- rep(0, n_batch)
    } else {
      delta_hat <- glm_f$coefficients$precision - as.numeric(alpha_z)
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
  cat(sprintf("Errors encountered in %s features with model fitting; 
              these features won't be adjusted for batch effects.\n",
              n_moderr))
  if (!mean.only) {
    cat(sprintf("Found %s features with zero model variance within at least one batch; 
                these features won't be adjusted for batch effects.\n",
                n_zero_modvar_batch))
  }
  
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
    mcint_fun <- monte_carlo_int_beta
    monte_carlo_res <- lapply(1:n_batch, function(ii) {
      if (ii == 1) {
        mcres <- mcint_fun(dat = bv[, batches_ind[[ii]]], 
                           mu = mu_hat_mat[, batches_ind[[ii]]], 
                           gamma = gamma_hat_mat[, ii], 
                           phi = phi_hat_mat[, batches_ind[[ii]]],
                           delta = delta_hat_mat[, ii],
                           feature.subset.n = feature.subset.n)
      } else {
        invisible(capture.output(mcres <- mcint_fun(dat = bv[, batches_ind[[ii]]], 
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
    
    ## set gamma and delta equal to 0 for reference batch (probably unnecessary, 
    ## but just to make sure)
    if (!is.null(ref.batch)) {
      gamma_star_mat[, ref] <- 0
      delta_star_mat[, ref] <- 0
    }
    
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
  mu_star_mat <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
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
  adj_bv_raw <- matrix(NA, nrow = nrow(bv), ncol = ncol(bv))
  for (kk in 1:n_batch) {
    bv_sub <- bv[, batches_ind[[kk]]]
    old_mu <- mu_hat_mat[, batches_ind[[kk]]]
    old_phi <- phi_hat_mat[, batches_ind[[kk]]] * exp(delta_hat_mat[, kk])
    new_mu <- mu_star_mat[, batches_ind[[kk]]]
    new_phi <- phi_star_mat[, batches_ind[[kk]]]
    adj_bv_raw[, batches_ind[[kk]]] <- match_quantiles_beta(bv_sub = bv_sub,
                                                            old_mu = old_mu, 
                                                            old_phi = old_phi, 
                                                            new_mu = new_mu, 
                                                            new_phi = new_phi)
  }
  
  ## Add back features with beta-values unqualified for beta regression  
  ## so that dimensions won't change
  adj_bv <- bvOri
  adj_bv[keep, ] <- adj_bv_raw
  return(adj_bv)
}
