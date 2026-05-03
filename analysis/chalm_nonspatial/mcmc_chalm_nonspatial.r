## Helper Functions
rCMLG <- function(H, alpha, kappa, tol = 1e-12){
  repeat{
    m <- length(kappa)
    w <- log(rgamma(m, shape = alpha, rate = kappa))
    cand <- tryCatch(
      as.numeric(qr.solve(as.matrix(H), w, tol = tol)),
      error = function(e) rep(NA_real_, ncol(as.matrix(H)))
    )
    if(all(is.finite(cand))) return(cand)
  }
}

.is_identity <- function(M, tol = 1e-12) {
  if (is.null(dim(M))) return(FALSE)
  if (nrow(M) != ncol(M)) return(FALSE)
  n <- nrow(M)
  if (inherits(M, "Matrix")) M <- as.matrix(M)
  if (max(abs(diag(M) - 1)) > tol) return(FALSE)
  M_off <- M
  diag(M_off) <- 0
  max(abs(M_off)) <= tol
}

rmvn1 <- function(mu, Sigma){
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  eps <- 1e-10 * max(1, mean(diag(Sigma)))
  Sigma <- Sigma + diag(eps, nrow(Sigma))
  L <- chol(Sigma)
  as.numeric(mu + t(L) %*% rnorm(length(mu)))
}

psd_sqrt <- function(Q, tol = 1e-10){
  Qm <- as.matrix(Q)
  Qm <- (Qm + t(Qm))/2
  ev <- eigen(Qm, symmetric = TRUE)
  vals <- ev$values
  mx <- max(1, max(vals))
  vals[vals < tol * mx] <- 0
  vals <- pmax(vals, 0)
  U <- ev$vectors
  U %*% diag(sqrt(vals), nrow = length(vals)) %*% t(U)
}

dlMLG <- function(q, V, a, k){
  V <- Matrix::Matrix(V)
  Vinv <- solve(V)
  det <- as.numeric(determinant(V %*% t(V), logarithm = TRUE))
  det <- det[1] * det[2]
  ll <- -0.5 * det + sum(a * log(k)) - sum(lgamma(a)) +
    t(a) %*% Vinv %*% q - t(k) %*% exp(Vinv %*% q)
  return(as.numeric(ll))
}

rTruncCMLG <- function(H=matrix(rnorm(6),3), alpha=c(1,1,1), kappa=c(1,1,1), cut=0 ){
  repeat {
    emp <- rCMLG(H=H, alpha=alpha, kappa=kappa)
    if(all(emp > cut)) break
  }
  return(emp)
}

## KNN Imputation
knn_constrained_allocation_by_state <- function(counties_sf, value_col, state_col = "STUPS", all_other_totals, k = 5, crs_projected = 5070, distance_power = 1) {
  stopifnot(inherits(counties_sf, "sf"))
  all_other_totals <- unlist(all_other_totals)
  if (is.null(names(all_other_totals))) {stop("all_other_totals must be a named list or named vector, e.g. list(IL = ..., IA = ...).")}
  original_crs <- st_crs(counties_sf)
  counties_proj <- st_transform(counties_sf, crs_projected)
  pts <- st_point_on_surface(counties_proj)
  coords <- st_coordinates(pts)
  values <- counties_proj[[value_col]]
  states <- counties_proj[[state_col]]
  counties_proj$value_imputed <- values
  counties_proj$allocation_status <- ifelse(is.na(values), "missing_not_allocated", "observed")
  state_summary <- list()
  for (s in unique(states)) {
    state_idx <- which(states == s)
    state_values <- values[state_idx]
    observed_local <- which(!is.na(state_values))
    missing_local  <- which(is.na(state_values))
    if (length(missing_local) == 0) {
      state_summary[[s]] <- data.frame(state = s, n_missing = 0, target_total = NA_real_, allocated_total = NA_real_, difference = NA_real_)
      next
    }
    if (!(s %in% names(all_other_totals))) {stop(paste0("No all-other total supplied for state: ", s))}
    state_target <- as.numeric(all_other_totals[s])
    if (is.na(state_target) || state_target < 0) {stop(paste0("Invalid all-other total for state: ", s))}
    if (length(observed_local) == 0) {stop(paste0("No observed counties available for KNN in state: ", s))}
    observed_idx <- state_idx[observed_local]
    missing_idx  <- state_idx[missing_local]
    observed_coords <- coords[observed_idx, , drop = FALSE]
    missing_coords  <- coords[missing_idx, , drop = FALSE]
    observed_values <- values[observed_idx]
    k_use <- min(k, length(observed_idx))
    knn_result <- FNN::get.knnx(data = observed_coords, query = missing_coords, k = k_use)
    neighbor_ids  <- knn_result$nn.index
    neighbor_dist <- knn_result$nn.dist
    # Avoid division by zero
    neighbor_dist[neighbor_dist == 0] <- 1e-6
    weights <- 1 / (neighbor_dist ^ distance_power)
    weights <- weights / rowSums(weights)
    neighbor_values <- matrix(observed_values[neighbor_ids], nrow = length(missing_idx), ncol = k_use)
    raw_imputed <- rowSums(weights * neighbor_values)
    raw_sum <- sum(raw_imputed, na.rm = TRUE)
    if (raw_sum <= 0) {stop(paste0("Raw KNN predictions sum to zero or less in state: ", s))}
    scale_factor <- state_target / raw_sum
    scaled_imputed <- raw_imputed * scale_factor
    counties_proj$value_imputed[missing_idx] <- scaled_imputed
    counties_proj$allocation_status[missing_idx] <- "knn_constrained_state"
    allocated_total <- sum(counties_proj$value_imputed[missing_idx], na.rm = TRUE)
    state_summary[[s]] <- data.frame(state = s, n_missing = length(missing_idx), target_total = state_target, allocated_total = allocated_total, difference = allocated_total - state_target)
  }
  state_summary <- bind_rows(state_summary)
  counties_out <- st_transform(counties_proj, original_crs)
  return(counties_out)
}


## MCMC SAE
saeFit_con_multistate <- function(
    Y, X, Psi, S2, SS,
    a = 0.5, b = 0.5, w = 1000, rho = 1000, sig2b = 1000,
    iter = 3000, burn = 1000, alpha = 1, propSD = 10,
    cVec = NULL,
    aTot_state = NULL,   # ADDED: vector of state-specific benchmark totals
    state_id = NULL      # ADDED: state label for each county/area
){
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  df <- SS - 1
  XP <- cbind(X, Psi)
  
  rmvn1 <- function(mu, Sigma){
    mu <- as.numeric(mu)
    Sigma <- as.matrix(Sigma)
    Sigma <- (Sigma + t(Sigma)) / 2
    eps <- 1e-10 * mean(diag(Sigma))
    Sigma <- Sigma + diag(eps, nrow(Sigma))
    L <- chol(Sigma)
    as.numeric(mu + t(L) %*% rnorm(length(mu)))
  }
  
  # Checks
  has_any_constraint_inputs <- !is.null(cVec) || !is.null(aTot_state) || !is.null(state_id)
  constrained <- !(is.null(cVec) || is.null(aTot_state) || is.null(state_id))
  if(has_any_constraint_inputs && !constrained){stop("For the constrained multi-state model, provide cVec, state_id, and aTot_state together.")}
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, n)
    if(length(cVec) != n) stop("cVec is not length n.")
    if(length(state_id) != n) stop("state_id must have length n.")
    state_id <- as.character(state_id)
    state_levels <- unique(state_id)
    G <- length(state_levels)
    if(length(aTot_state) == 1 && G > 1){stop("aTot_state must contain one benchmark value for each state.")}
    if(!is.null(names(aTot_state))){
      if(!all(state_levels %in% names(aTot_state))){
        stop("Names of aTot_state must include all states appearing in state_id.")
      }
      aTot_lookup <- as.numeric(aTot_state[state_levels])
      names(aTot_lookup) <- state_levels
    } else {
      if(length(aTot_state) != G){
        stop("If aTot_state is unnamed, it must have length equal to the number of states.")
      }
      aTot_lookup <- as.numeric(aTot_state)
      names(aTot_lookup) <- state_levels
    }
    if(any(!is.finite(aTot_lookup))){
      stop("All state-specific benchmark totals in aTot_state must be finite.")
    }
    state_index <- split(seq_len(n), factor(state_id, levels = state_levels))
    
    # feasibility check is now done each state
    for(st in state_levels){
      idx <- state_index[[st]]
      if(sum(cVec[idx]) > aTot_lookup[st]){
        stop(
          paste0(
            "Infeasible constraints in state ", st,
            ": sum(cVec[state]) > aTot_state[state]."
          )
        )
      }
    }
  }
  
  ## initial values ----
  beta1 <- rep(0, p)
  sig2e1 <- var(as.numeric(Y - X %*% beta1)) #median(S2)
  beta2 <- rep(0, p)
  eta2  <- rep(0, r)
  sig2e2 <- 10
  sig2 <- as.numeric(S2)
  
  if(constrained){
    theta <- numeric(n)
    for(st in state_levels){
      idx <- state_index[[st]]
      n_st <- length(idx)
      aTot_st <- aTot_lookup[st]
      theta[idx] <- cVec[idx] + (aTot_st - sum(cVec[idx])) / n_st
    }
  } else {
    theta <- as.numeric(Y)
  }
  
  ## store values
  thetaOut  <- matrix(NA_real_, n, iter)
  beta1Out  <- matrix(NA_real_, p, iter)
  sig2e1Out <- numeric(iter)
  beta2Out  <- matrix(NA_real_, p, iter)
  eta2Out   <- matrix(NA_real_, r, iter)
  sig2Out   <- matrix(NA_real_, n, iter)
  sig2e2Out <- numeric(iter)
  ll <- numeric(iter)
  pb <- txtProgressBar(min = 0, max = iter, style = 3)
  
  for(i in 1:iter){
    ## Variance model
    H <- rbind(XP, XP, Matrix::bdiag(alpha^(-0.5)/sqrt(sig2b) * Matrix::Diagonal(p), alpha^(-0.5)/sqrt(sig2e2) * Matrix::Diagonal(r)))
    aVec <- c(rep(1/2, n), rep(df/2, n), rep(alpha, p + r))
    kVec <- c((Y - theta)^2 / 2, 0.5 * (SS - 1) * S2, rep(alpha, p + r))
    tmp <- rCMLG(H, aVec, kVec)
    beta2 <- beta2Out[, i] <- tmp[1:p]
    eta2  <- eta2Out[, i]  <- tmp[-c(1:p)]
    g <- as.numeric(XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <- log(S2) + g
    sig2 <- sig2Out[, i] <- exp(log_sig2)
    
    ## Constrained theta model
    mu_theta <- as.numeric(X %*% beta1)
    if(!constrained){
      v <- 1 / (1/sig2 + 1/sig2e1)
      m <- v * (Y/sig2 + mu_theta/sig2e1)
      theta <- rnorm(n, mean = m, sd = sqrt(v))
    } else {
      for(st in state_levels){
        idx <- state_index[[st]]
        m_st <- length(idx)
        aTot_st <- aTot_lookup[st]
        if(m_st == 1){
          theta[idx] <- aTot_st
          if(theta[idx] < cVec[idx]){
            stop(paste0("Single-area state ", st, " violates lower bound."))
          }
          next
        }
        dth_local_pos <- m_st
        dth_idx <- idx[dth_local_pos]
        v_dth <- 1 / (1/sig2[dth_idx] + 1/sig2e1)
        m_dth <- v_dth * (Y[dth_idx]/sig2[dth_idx] + mu_theta[dth_idx]/sig2e1)
        max_reject <- 2000
        accept <- FALSE
        for(rr in 1:max_reject){
          total_first <- sum(theta[idx[1:(m_st - 1)]])
          bad <- FALSE
          for(jj in 1:(m_st - 1)){
            ii <- idx[jj]
            v_i <- 1 / (1/sig2[ii] + 1/sig2e1)
            m_i <- v_i * (Y[ii]/sig2[ii] + mu_theta[ii]/sig2e1)
            S_i <- total_first - theta[ii]
            v_tilde <- 1 / (1/v_i + 1/v_dth)
            m_tilde <- v_tilde * (m_i / v_i + (aTot_st - S_i - m_dth) / v_dth)
            
            # upper truncation bound uses the state's benchmark and its corresponding balancing county
            U_i <- aTot_st - cVec[dth_idx] - S_i
            if(!is.finite(U_i) || U_i <= cVec[ii]){
              bad <- TRUE
              break
            }
            theta_new <- truncnorm::rtruncnorm(1, a = cVec[ii], b = U_i, mean = m_tilde, sd = sqrt(v_tilde))
            total_first <- S_i + theta_new
            theta[ii] <- theta_new
          }
          if(bad) next
          
          # enforce the state-specific benchmark exactly
          theta[dth_idx] <- aTot_st - sum(theta[idx[1:(m_st - 1)]])
          if(theta[dth_idx] < cVec[dth_idx]) next
          accept <- TRUE
          break
        }
        
        if(!accept){stop(paste0("theta-update failed to find a feasible draw in state ", st, ". Check cVec/aTot_state."))}
      } # end state loop
    }
    
    if(constrained){
      # verify constraints state by state
      for(st in state_levels){
        idx <- state_index[[st]]
        if(abs(sum(theta[idx]) - aTot_lookup[st]) > 1e-8){stop(paste0("sum(theta[state]) != aTot_state[state] for state ", st))}
      }
      if(any(theta < cVec - 1e-10)){stop("theta violates lower bounds")}
    }
    thetaOut[, i] <- theta
    
    ## Update beta_1 | theta, sig2e1
    Qb1 <- (crossprod(X) / sig2e1) + diag(1/sig2b, p)
    Qb1 <- (Qb1 + t(Qb1)) / 2
    Vb1 <- solve(Qb1)
    mb1 <- Vb1 %*% (crossprod(X, theta) / sig2e1)
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, i] <- beta1
    
    ## Update sigma_theta^2
    shape <- a + n/2
    rate  <- b + 0.5 * sum((theta - as.numeric(X %*% beta1))^2)
    sig2e1 <- 1 / rgamma(1, shape = shape, rate = rate)
    sig2e1Out[i] <- sig2e1
    
    ## sample RE variances
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    log_dtr <- function(x, a, b, mean = 0, sd = 1){log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))}
    log_target_prop <- dlMLG(q = eta2, V = sqrt(alpha) * sqrt(prop) * diag(r), a = rep(alpha, r), k = rep(alpha, r))
    log_target_curr <- dlMLG(q = eta2, V = sqrt(alpha) * sqrt(sig2e2) * diag(r), a = rep(alpha, r), k = rep(alpha, r))
    log_prior_prop <- log_dtr(sqrt(prop), a = 0, b = Inf, mean = 0, sd = 10)
    log_prior_curr <- log_dtr(sqrt(sig2e2), a = 0, b = Inf, mean = 0, sd = 10)
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    ratio <- (log_target_prop - log_target_curr) + (log_prior_prop  - log_prior_curr) + (log_q_backward  - log_q_forward)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2

    ll[i] <- -2 * sum(dnorm(Y, mean = theta, sd = sqrt(sig2), log = TRUE))
    setTxtProgressBar(pb, i)
  }
  
  keep <- (burn + 1):iter
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(Y, mean = rowMeans(thetaOut[, keep, drop = FALSE]), sd = sqrt(rowMeans(sig2Out[, keep, drop = FALSE])), log = TRUE))
  
  list(
    Y_hat = Y,
    Preds = thetaOut,
    DIC = DIC,
    sig2 = sig2Out,
    beta1 = beta1Out,
    sig2e1 = sig2e1Out,
    beta2 = beta2Out,
    eta2 = eta2Out,
    sig2e2 = sig2e2Out
  )
}
