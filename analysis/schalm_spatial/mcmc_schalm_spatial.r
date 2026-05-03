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


## Incidence Matrix
make_incidence <- function(W){
  W <- Matrix::Matrix(W, sparse = TRUE)                                
  W <- (W + t(W)) / 2                                                 
  Matrix::diag(W) <- 0                                                
  U  <- Matrix::triu(W, k = 1)                                    
  ij <- Matrix::summary(U)                                             
  m  <- nrow(ij)                                                       
  n  <- nrow(W)                                                        
  E_inc <- Matrix::sparseMatrix(i = rep(seq_len(m), 2), j = c(ij$i, ij$j), x = c(rep(1, m), rep(-1, m)), dims = c(m, n))
  E_inc
}


## Theta Constraint Function
theta_constraint <- function(theta, A, B, cVec, aTot){
  n <- length(theta)
  d <- n
  idx <- 1:(d-1)
  A11 <- A[idx, idx, drop = FALSE]
  a12 <- as.numeric(A[idx, d])
  a22 <- as.numeric(A[d, d])
  b1  <- as.numeric(B[idx])
  b2  <- as.numeric(B[d])
  qz <- b1 - aTot * a12 + rep(aTot * a22 - b2, d - 1)
  z <- theta[idx]
  z <- pmax(z, cVec[idx])
  cap <- aTot - cVec[d]
  if(sum(z) > cap){
    z <- z - (sum(z) - cap) / (d - 1)
    z <- pmax(z, cVec[idx])
  }
  for(i in idx){
    S_minus <- sum(z[-i])
    U_i <- aTot - cVec[d] - S_minus
    if(!is.finite(U_i) || U_i <= cVec[i]) stop("Infeasible bounds")
    s <- sum(z)
    t12 <- sum(a12 * z)
    Az_i <- as.numeric(A11[i, ] %*% z)
    Qdiag_i <- as.numeric(A11[i, i]) - 2 * a12[i] + a22
    if(!is.finite(Qdiag_i) || Qdiag_i <= 0) stop("Nonpositive conditional precision in theta-update.")
    Qzz_i <- Az_i - t12 - a12[i] * s + a22 * s
    sum_off <- Qzz_i - Qdiag_i * z[i]
    mu_i <- (qz[i] - sum_off) / Qdiag_i
    sd_i <- sqrt(1 / Qdiag_i)
    z[i] <- truncnorm::rtruncnorm(1, a = cVec[i], b = U_i, mean = mu_i, sd = sd_i)
  }
  theta_new <- numeric(d)
  theta_new[idx] <- z
  theta_new[d] <- aTot - sum(z)
  if(theta_new[d] < cVec[d] - 1e-10) stop("theta_d < c_d after constrained update.")
  if(any(theta_new < cVec - 1e-10)) stop("theta violates lower bounds after constrained update.")
  if(abs(sum(theta_new) - aTot) > 1e-8) stop("sum(theta) != aTot after constrained update.")
  theta_new
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


## MCMC SAEICAR2 
saeFit_con_icar2_incidence_obs <- function(
    Y, X, Psi = NULL, S2, SS,        
    W,                               
    obs_idx,                         
    a = 0.5, b = 0.5,
    sig2b1 = 1000,
    sig2b2 = 1000,
    alpha = 1,
    propSD = 10,
    iter = 3000, burn = 1000,
    cVec = NULL, aTot = NULL,
    ridge_var = 1e-6,
    progress = TRUE
){
  ## Dimensions and indexing
  W <- Matrix::Matrix(W, sparse = TRUE)
  W <- (W + t(W)) / 2
  Matrix::diag(W) <- 0
  N_all <- nrow(W)                                                     
  X <- as.matrix(X)                                                    
  if(nrow(X) != N_all) stop("X must have one row for each county in W.")
  p <- ncol(X)
  obs_idx <- as.integer(obs_idx)                                       
  n_obs <- length(obs_idx)                                             
  if(length(Y) != n_obs) stop("length(Y) must equal length(obs_idx).")
  if(anyDuplicated(obs_idx)) stop("obs_idx contains duplicates.")
  if(any(obs_idx < 1 | obs_idx > N_all)) stop("obs_idx out of range.")
  obs <- obs_idx                                                       
  mis <- setdiff(seq_len(N_all), obs)                                  
  n_mis <- length(mis)                                                 
  X_obs <- X[obs, , drop = FALSE]                                      
  
  ## Psi handling
  if(!is.null(Psi)){warning("observed-county selector P_obs is used instead.")}
  P_obs <- Matrix::sparseMatrix(i = seq_len(n_obs), j = obs, x = 1, dims = c(n_obs, N_all))
  
  ## Data checks
  constrained <- !(is.null(cVec) || is.null(aTot))
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, N_all)                     
    if(length(cVec) != N_all) stop("cVec must be length N_all (or scalar).")
    if(!is.finite(aTot)) stop("aTot must be a finite scalar.")
    if(sum(cVec[obs]) > aTot){stop("Infeasible constraints: sum(cVec[obs_idx]) > aTot.")}
  }
  df_vec <- as.numeric(SS - 1)
  if(length(df_vec) == 1) df_vec <- rep(df_vec, n_obs)                 
  if(length(df_vec) != n_obs) stop("SS must be scalar or length n_obs.")
  S2 <- as.numeric(S2)
  if(length(S2) == 1) S2 <- rep(S2, n_obs)                             
  if(length(S2) != n_obs) stop("S2 must be scalar or length n_obs.")
  
  ## Incidence Matrix Graph
  di <- Matrix::rowSums(W)
  E_inc <- make_incidence(W)                                          
  Q <- Matrix::crossprod(E_inc)                                       
  nb <- lapply(1:N_all, function(i) which(W[i, ] != 0))
  count_components <- function(nb){
    n0 <- length(nb)
    visited <- rep(FALSE, n0)
    comp <- 0
    for(i in 1:n0){
      if(!visited[i]){
        comp <- comp + 1
        stack <- i
        visited[i] <- TRUE
        while(length(stack)){
          v <- stack[[length(stack)]]
          stack <- stack[-length(stack)]
          for(u in nb[[v]]){
            if(!visited[u]){
              visited[u] <- TRUE
              stack <- c(stack, u)
            }
          }
        }
      }
    }
    comp
  }
  nComp <- count_components(nb)
  rankQ <- N_all - nComp                                               
  
  ## Variance-side precision construction
  D_var <- rbind(E_inc, sqrt(ridge_var) * Matrix::Diagonal(N_all))
  Q_var <- Matrix::crossprod(D_var)                                    
  Qsqrt <- psd_sqrt(Q_var)
  Qsqrt <- Matrix::Matrix(Qsqrt)
  
  ## Initial values
  beta1  <- rep(0, p)
  sig2e1 <- median(S2)
  beta2  <- rep(0, p)
  eta2   <- rep(0, N_all)                                              
  sig2e2 <- 10
  if(constrained){
    theta <- rep(mean(Y), N_all)
    if(n_obs == 1){
      theta_obs_init <- aTot
      if(theta_obs_init < cVec[obs] - 1e-10){
        stop("Initial observed theta not feasible for n_obs = 1.")
      }
    } else {
      theta_obs_init <- cVec[obs] + (aTot - sum(cVec[obs])) / n_obs
      theta_obs_init[n_obs] <- aTot - sum(theta_obs_init[-n_obs])
    }
    if(any(theta_obs_init < cVec[obs] - 1e-10)){stop("Initial observed theta violates cVec on obs_idx.")}
    if(abs(sum(theta_obs_init) - aTot) > 1e-8){stop("Initial observed theta does not sum to aTot.")}
    theta[obs] <- theta_obs_init                                       
    if(n_mis > 0) theta[mis] <- mean(theta_obs_init)                   
  } else {
    theta <- rep(mean(Y), N_all)                                       
    theta[obs] <- as.numeric(Y)                                        
  }
  
  ## Storing Values
  thetaOut  <- matrix(NA_real_, N_all, iter)                           
  beta1Out  <- matrix(NA_real_, p, iter)
  sig2e1Out <- numeric(iter)
  beta2Out  <- matrix(NA_real_, p, iter)
  eta2Out   <- matrix(NA_real_, N_all, iter)                          
  sig2Out   <- matrix(NA_real_, N_all, iter)                           
  sig2e2Out <- numeric(iter)
  ll <- numeric(iter)
  
  if(progress) pb <- txtProgressBar(min = 0, max = iter, style = 3)
  
  eps_k <- 1e-12
  EX_all <- E_inc %*% X                                                
  XtQX   <- Matrix::crossprod(EX_all)                                  
  
  ## MCMC
  for(s in 1:iter){
    ## Variance model
    theta_obs <- theta[obs]                                            
    eta2_obs  <- as.numeric(P_obs %*% eta2)                            
    expEta_obs <- exp(eta2_obs)                                        
    H_beta2 <- rbind(X_obs, X_obs, alpha^(-0.5)/sqrt(sig2b2) * Matrix::Diagonal(p))
    a_beta2 <- c(rep(1/2, n_obs), df_vec/2, rep(alpha, p))             
    k_beta2 <- c(0.5 * (Y - theta_obs)^2 * expEta_obs, (df_vec/2) * S2 * expEta_obs, rep(alpha, p))
    k_beta2 <- pmax(k_beta2, eps_k)
    beta2 <- rCMLG(H_beta2, a_beta2, k_beta2)
    beta2Out[, s] <- beta2
    expXBeta2_obs <- exp(as.numeric(X_obs %*% beta2))                  
    H_eta2 <- rbind(P_obs, P_obs, alpha^(-0.5)/sqrt(sig2e2) * Qsqrt)
    a_eta2 <- c(rep(1/2, n_obs), df_vec/2, rep(alpha, N_all))          
    k_eta2 <- c(0.5 * (Y - theta_obs)^2 * expXBeta2_obs, (df_vec/2) * S2 * expXBeta2_obs, rep(alpha, N_all))
    k_eta2 <- pmax(k_eta2, eps_k)
    eta2 <- rCMLG(H_eta2, a_eta2, k_eta2)
    eta2Out[, s] <- eta2
    g_all <- as.numeric(X %*% beta2 + eta2)                            
    g_all <- g_all - mean(g_all)                                       
    sig2_obs <- exp(log(S2) + g_all[obs])                        
    sig2_full <- rep(NA_real_, N_all)                                  
    sig2_full[obs] <- sig2_obs                                         
    sig2Out[, s] <- sig2_full                                          
    
    ## Constrained theta model
    prec_y <- numeric(N_all)
    prec_y[obs] <- 1 / sig2_obs
    b_y <- numeric(N_all)
    b_y[obs] <- Y / sig2_obs
    A <- Matrix::Diagonal(x = prec_y) + (1/sig2e1) * Q
    B <- b_y + (1/sig2e1) * as.numeric(Q %*% (X %*% beta1))
    
    if(!constrained){
      A_dense <- as.matrix(A)
      R <- chol((A_dense + t(A_dense))/2)
      mu <- backsolve(R, forwardsolve(t(R), as.numeric(B)))
      z <- rnorm(N_all)
      v <- backsolve(R, z)
      theta <- as.numeric(mu + v)
    } else {
      Aoo <- as.matrix(A[obs, obs, drop = FALSE])                      
      if(n_mis > 0){
        Aom <- as.matrix(A[obs, mis, drop = FALSE])                    
        Amo <- as.matrix(A[mis, obs, drop = FALSE])                    
        Amm <- as.matrix(A[mis, mis, drop = FALSE])                    
        Bobs_cond <- as.numeric(B[obs] - Aom %*% theta[mis])
      } else {
        ## no unobserved counties
        Bobs_cond <- as.numeric(B[obs])
      }
      ## Constrained update on observed counties.
      if(n_obs == 1){
        theta[obs] <- aTot
        if(theta[obs] < cVec[obs] - 1e-10){
          stop("Observed benchmark violates lower bound for n_obs = 1.")
        }
      } else {
        theta[obs] <- theta_constraint(theta = theta[obs], A = Aoo, B = Bobs_cond, cVec = cVec[obs], aTot = aTot)
      }
      
      ## Update unobserved counties conditionally with no constraints but just to carry the ICAR spatial structure.
      if(n_mis > 0){
        Bmis_cond <- as.numeric(B[mis] - Amo %*% theta[obs])
        Amm <- (Amm + t(Amm)) / 2                                      
        Rmis <- chol(Amm)
        mu_mis <- backsolve(Rmis, forwardsolve(t(Rmis), Bmis_cond))
        zmis <- rnorm(n_mis)
        vmis <- backsolve(Rmis, zmis)
        theta[mis] <- as.numeric(mu_mis + vmis)
      }
    }
    thetaOut[, s] <- theta
    
    ## Update beta_1 | theta, sig2e1
    Etheta <- E_inc %*% theta                                          
    XtQtheta <- Matrix::crossprod(EX_all, Etheta)                      
    Qb1 <- (XtQX / sig2e1) + diag(1/sig2b1, p)
    Qb1 <- (Qb1 + t(Qb1))/2
    Vb1 <- solve(as.matrix(Qb1))
    mb1 <- Vb1 %*% (as.numeric(XtQtheta) / sig2e1)
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, s] <- beta1
    
    ## Update sigma_theta^2
    resid1  <- theta - as.numeric(X %*% beta1)
    Eresid1 <- E_inc %*% resid1                                        
    quad1   <- sum(Eresid1^2)                                          
    shape1 <- a + rankQ/2
    rate1  <- b + 0.5 * quad1
    sig2e1 <- 1 / rgamma(1, shape = shape1, rate = rate1)
    sig2e1Out[s] <- sig2e1
    
    ## Sample RE variance parameter
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    log_dtr <- function(x, a, b, mean = 0, sd = 1){log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))}
    log_target <- function(sig2_val){
      sigma <- sqrt(sig2_val)
      if(!is.finite(sigma) || sigma <= 0) return(-Inf)
      u  <- as.numeric((alpha^(-0.5)/sigma) * (Qsqrt %*% eta2))
      ll <- -N_all * log(sigma) + alpha * sum(u) - alpha * sum(exp(u)) 
      lp <- log_dtr(sigma, a = 0, b = Inf, mean = 0, sd = 10) - log(2 * sigma) #added - log(2 * sigma)
      ll + lp
    }
    
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    ratio <- (log_target(prop) - log_target(sig2e2)) + (log_q_backward - log_q_forward)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[s] <- sig2e2
    
    ## Fit criterion
    ll[s] <- -2 * sum(dnorm(Y, mean = theta[obs], sd = sqrt(sig2_obs), log = TRUE))  
    if(progress) setTxtProgressBar(pb, s)
  }
  if(progress) close(pb)
  keep <- (burn + 1):iter
  theta_bar_obs <- rowMeans(thetaOut[obs, keep, drop = FALSE])         
  sig2_bar_obs  <- rowMeans(sig2Out[obs, keep, drop = FALSE], na.rm = TRUE)  
  
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(Y, mean = theta_bar_obs, sd = sqrt(sig2_bar_obs), log = TRUE))
  
  list(
    Y_hat      = Y,
    obs_idx    = obs,                                                 
    Preds      = thetaOut, #[obs, , drop = FALSE],                      
    DIC        = DIC,
    sig2       = sig2Out,                                              
    beta1      = beta1Out,
    sig2e1     = sig2e1Out,
    beta2      = beta2Out,
    eta2       = eta2Out,
    sig2e2     = sig2e2Out,
    W          = W,
    E_inc      = E_inc,                                                
    D_var      = D_var,                                                
    Q          = Q,
    Q_var      = Q_var,                                                
    di         = di,
    rankQ      = rankQ,
    ridge_var  = ridge_var,
    N_all      = N_all,                                                
    n_obs      = n_obs                                                 
  )
}



## MCMC SAEICAR2 - Multistate 
saeFit_con_icar2_incidence_obs_multistate <- function(
    Y, X, Psi = NULL, S2, SS,        
    W,                               
    obs_idx,                         
    a = 0.5, b = 0.5,
    sig2b1 = 1000,
    sig2b2 = 1000,
    alpha = 1,
    propSD = 10,
    iter = 3000, burn = 1000,
    cVec = NULL,
    aTot_state = NULL,               
    state_id = NULL,                 
    ridge_var = 1e-6,
    progress = TRUE
){
  ## Dimensions and indexing
  W <- Matrix::Matrix(W, sparse = TRUE)
  W <- (W + t(W)) / 2
  Matrix::diag(W) <- 0
  N_all <- nrow(W)                                                     
  X <- as.matrix(X)                                                    
  if(nrow(X) != N_all) stop("X must have one row for each county in W.")
  p <- ncol(X)
  obs_idx <- as.integer(obs_idx)                                       
  n_obs <- length(obs_idx)                                             
  if(length(Y) != n_obs) stop("length(Y) must equal length(obs_idx).")
  if(anyDuplicated(obs_idx)) stop("obs_idx contains duplicates.")
  if(any(obs_idx < 1 | obs_idx > N_all)) stop("obs_idx out of range.")
  obs <- obs_idx                                                       
  mis <- setdiff(seq_len(N_all), obs)                                  
  n_mis <- length(mis)                                                 
  X_obs <- X[obs, , drop = FALSE]                                      
  
  ## Psi handling
  if(!is.null(Psi)){warning("observed-county selector P_obs is used instead.")}
  P_obs <- Matrix::sparseMatrix(i = seq_len(n_obs), j = obs, x = 1, dims = c(n_obs, N_all))
  
  ## Data checks: multistate constraint inputs
  has_any_constraint_inputs <- !is.null(cVec) || !is.null(aTot_state) || !is.null(state_id)
  constrained <- !(is.null(cVec) || is.null(aTot_state) || is.null(state_id))
  if(has_any_constraint_inputs && !constrained){stop("For the constrained multi-state model, provide cVec, state_id, and aTot_state together.")}
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, N_all)                     
    if(length(cVec) != N_all) stop("cVec must be length N_all (or scalar).")
    if(length(state_id) != N_all) stop("state_id must have length N_all.")
    state_id <- as.character(state_id)
    state_id_obs <- state_id[obs]
    state_levels_obs <- unique(state_id_obs)
    G_obs <- length(state_levels_obs)
    if(length(aTot_state) == 1 && G_obs > 1){stop("aTot_state must contain one benchmark value for each observed state.")}
    if(!is.null(names(aTot_state))){
      if(!all(state_levels_obs %in% names(aTot_state))){
        stop("Names of aTot_state must include all states appearing in state_id[obs_idx].")
      }
      aTot_lookup <- as.numeric(aTot_state[state_levels_obs])
      names(aTot_lookup) <- state_levels_obs
    } else {
      if(length(aTot_state) != G_obs){
        stop("If aTot_state is unnamed, it must have length equal to the number of observed states.")
      }
      aTot_lookup <- as.numeric(aTot_state)
      names(aTot_lookup) <- state_levels_obs
    }
    if(any(!is.finite(aTot_lookup))){stop("All state-specific benchmark totals in aTot_state must be finite.")}
    ## Split observed counties by state
    state_index_obs <- split(obs, factor(state_id_obs, levels = state_levels_obs))
    ## Feasibility check for each state
    for(st in state_levels_obs){
      idx_st <- state_index_obs[[st]]
      if(sum(cVec[idx_st]) > aTot_lookup[st]){
        stop(
          paste0(
            "Infeasible constraints in observed state ", st,
            ": sum(cVec[state ∩ obs]) > aTot_state[state]."
          )
        )
      }
    }
  }
  df_vec <- as.numeric(SS - 1)
  if(length(df_vec) == 1) df_vec <- rep(df_vec, n_obs)                 
  if(length(df_vec) != n_obs) stop("SS must be scalar or length n_obs.")
  S2 <- as.numeric(S2)
  if(length(S2) == 1) S2 <- rep(S2, n_obs)                             
  if(length(S2) != n_obs) stop("S2 must be scalar or length n_obs.")
  
  ## ----------------------------------------------------------
  ## Incidence Matrix Graph objects
  ## ----------------------------------------------------------
  di <- Matrix::rowSums(W)
  E_inc <- make_incidence(W)                                           
  Q <- Matrix::crossprod(E_inc)                                        
  nb <- lapply(1:N_all, function(i) which(W[i, ] != 0))
  count_components <- function(nb){
    n0 <- length(nb)
    visited <- rep(FALSE, n0)
    comp <- 0
    for(i in 1:n0){
      if(!visited[i]){
        comp <- comp + 1
        stack <- i
        visited[i] <- TRUE
        while(length(stack)){
          v <- stack[[length(stack)]]
          stack <- stack[-length(stack)]
          for(u in nb[[v]]){
            if(!visited[u]){
              visited[u] <- TRUE
              stack <- c(stack, u)
            }
          }
        }
      }
    }
    comp
  }
  nComp <- count_components(nb)
  rankQ <- N_all - nComp                                               
  
  ## Variance-side precision construction
  D_var <- rbind(E_inc, sqrt(ridge_var) * Matrix::Diagonal(N_all))
  Q_var <- Matrix::crossprod(D_var)                                    
  Qsqrt <- psd_sqrt(Q_var)
  Qsqrt <- Matrix::Matrix(Qsqrt)
  
  ## Initial values
  beta1  <- rep(0, p)
  sig2e1 <- median(S2)
  beta2  <- rep(0, p)
  eta2   <- rep(0, N_all)                                              
  sig2e2 <- 10
  if(constrained){
    theta <- rep(mean(Y), N_all)
    for(st in state_levels_obs){
      idx_st <- state_index_obs[[st]]
      n_st <- length(idx_st)
      aTot_st <- aTot_lookup[st]
      if(n_st == 1){
        theta_st_init <- aTot_st
        if(theta_st_init < cVec[idx_st] - 1e-10){
          stop(paste0("Initial observed theta not feasible for single-county state ", st, "."))
        }
      } else {
        theta_st_init <- cVec[idx_st] + (aTot_st - sum(cVec[idx_st])) / n_st
        theta_st_init[n_st] <- aTot_st - sum(theta_st_init[-n_st])
      }
      if(any(theta_st_init < cVec[idx_st] - 1e-10)){stop(paste0("Initial observed theta violates cVec in state ", st, "."))}
      if(abs(sum(theta_st_init) - aTot_st) > 1e-8){stop(paste0("Initial observed theta does not sum to aTot_state in state ", st, "."))}
      theta[idx_st] <- theta_st_init
    }
    if(n_mis > 0) theta[mis] <- mean(theta[obs])                       
  } else {
    theta <- rep(mean(Y), N_all)                                       
    theta[obs] <- as.numeric(Y)                                        
  }
  
  ## Storage
  thetaOut  <- matrix(NA_real_, N_all, iter)                           
  beta1Out  <- matrix(NA_real_, p, iter)
  sig2e1Out <- numeric(iter)
  beta2Out  <- matrix(NA_real_, p, iter)
  eta2Out   <- matrix(NA_real_, N_all, iter)                           
  sig2Out   <- matrix(NA_real_, N_all, iter)                           
  sig2e2Out <- numeric(iter)
  ll <- numeric(iter)
  if(progress) pb <- txtProgressBar(min = 0, max = iter, style = 3)
  
  eps_k <- 1e-12
  EX_all <- E_inc %*% X                                                
  XtQX   <- Matrix::crossprod(EX_all)                                  
  
  ## MCMC
  for(s in 1:iter){
    ## Variance model
    theta_obs <- theta[obs]                                            
    eta2_obs  <- as.numeric(P_obs %*% eta2)                            
    expEta_obs <- exp(eta2_obs)
    H_beta2 <- rbind(X_obs, X_obs, alpha^(-0.5)/sqrt(sig2b2) * Matrix::Diagonal(p))
    a_beta2 <- c(rep(1/2, n_obs), df_vec/2, rep(alpha, p))
    k_beta2 <- c(0.5 * (Y - theta_obs)^2 * expEta_obs, (df_vec/2) * S2 * expEta_obs, rep(alpha, p))
    k_beta2 <- pmax(k_beta2, eps_k)
    beta2 <- rCMLG(H_beta2, a_beta2, k_beta2)
    beta2Out[, s] <- beta2
    expXBeta2_obs <- exp(as.numeric(X_obs %*% beta2))
    H_eta2 <- rbind(P_obs, P_obs, alpha^(-0.5)/sqrt(sig2e2) * Qsqrt)
    a_eta2 <- c(rep(1/2, n_obs), df_vec/2, rep(alpha, N_all))
    k_eta2 <- c(0.5 * (Y - theta_obs)^2 * expXBeta2_obs, (df_vec/2) * S2 * expXBeta2_obs, rep(alpha, N_all))
    k_eta2 <- pmax(k_eta2, eps_k)
    eta2 <- rCMLG(H_eta2, a_eta2, k_eta2)
    eta2Out[, s] <- eta2
    g_all <- as.numeric(X %*% beta2 + eta2)                            
    g_all <- g_all - mean(g_all)
    sig2_obs <- exp(log(S2) + g_all[obs])                              
    sig2_full <- rep(NA_real_, N_all)
    sig2_full[obs] <- sig2_obs
    sig2Out[, s] <- sig2_full
    
    ## Constrained theta model: Only observed counties contribute direct-data precision, but the ICAR precision is still full-graph Q.
    prec_y <- numeric(N_all)
    prec_y[obs] <- 1 / sig2_obs
    b_y <- numeric(N_all)
    b_y[obs] <- Y / sig2_obs
    A <- Matrix::Diagonal(x = prec_y) + (1/sig2e1) * Q
    B <- b_y + (1/sig2e1) * as.numeric(Q %*% (X %*% beta1))
    if(!constrained){
      A_dense <- as.matrix(A)
      R <- chol((A_dense + t(A_dense))/2)
      mu <- backsolve(R, forwardsolve(t(R), as.numeric(B)))
      z <- rnorm(N_all)
      v <- backsolve(R, z)
      theta <- as.numeric(mu + v)
    } else {
      for(st in state_levels_obs){
        idx_st <- state_index_obs[[st]]
        n_st <- length(idx_st)
        aTot_st <- aTot_lookup[st]
        if(n_st == 1){
          theta[idx_st] <- aTot_st
          if(theta[idx_st] < cVec[idx_st] - 1e-10){
            stop(paste0("Observed state benchmark violates lower bound in state ", st, "."))
          }
          next
        }
        rest_st <- setdiff(seq_len(N_all), idx_st)
        A_st <- as.matrix(A[idx_st, idx_st, drop = FALSE])            
        if(length(rest_st) > 0){
          A_st_rest <- as.matrix(A[idx_st, rest_st, drop = FALSE])     
          B_st_cond <- as.numeric(B[idx_st] - A_st_rest %*% theta[rest_st])
        } else {
          B_st_cond <- as.numeric(B[idx_st])
        }
        theta[idx_st] <- theta_constraint(theta = theta[idx_st], A = A_st, B = B_st_cond, cVec = cVec[idx_st], aTot = aTot_st)
      }
      
      ## Update unobserved counties conditionally with no benchmar but only to carry the ICAR spatial structure.
      if(n_mis > 0){
        Amo <- as.matrix(A[mis, obs, drop = FALSE])                    
        Amm <- as.matrix(A[mis, mis, drop = FALSE])
        Bmis_cond <- as.numeric(B[mis] - Amo %*% theta[obs])
        Amm <- (Amm + t(Amm)) / 2
        Rmis <- chol(Amm)
        mu_mis <- backsolve(Rmis, forwardsolve(t(Rmis), Bmis_cond))
        zmis <- rnorm(n_mis)
        vmis <- backsolve(Rmis, zmis)
        theta[mis] <- as.numeric(mu_mis + vmis)
      }
      
      ## Verify state-specific observed benchmarks after the update.
      for(st in state_levels_obs){
        idx_st <- state_index_obs[[st]]
        if(abs(sum(theta[idx_st]) - aTot_lookup[st]) > 1e-8){
          stop(paste0("sum(theta[state ∩ obs]) != aTot_state[state] for state ", st))
        }
      }
      if(any(theta[obs] < cVec[obs] - 1e-10)){
        stop("Observed theta violates lower bounds after multistate constrained update.")
      }
    }
    thetaOut[, s] <- theta
    
    ## Update beta_1 | theta, sig2e1
    Etheta <- E_inc %*% theta                                          
    XtQtheta <- Matrix::crossprod(EX_all, Etheta)                      
    Qb1 <- (XtQX / sig2e1) + diag(1/sig2b1, p)
    Qb1 <- (Qb1 + t(Qb1)) / 2
    Vb1 <- solve(as.matrix(Qb1))
    mb1 <- Vb1 %*% (as.numeric(XtQtheta) / sig2e1)
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, s] <- beta1
    
    ## Update sigma_theta^2
    resid1  <- theta - as.numeric(X %*% beta1)
    Eresid1 <- E_inc %*% resid1                                        
    quad1   <- sum(Eresid1^2)                                          
    shape1 <- a + rankQ/2
    rate1  <- b + 0.5 * quad1
    sig2e1 <- 1 / rgamma(1, shape = shape1, rate = rate1)
    sig2e1Out[s] <- sig2e1
    
    ## Sample RE variance parameter
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    log_dtr <- function(x, a, b, mean = 0, sd = 1){log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))}
    log_target <- function(sig2_val){
      sigma <- sqrt(sig2_val)
      if(!is.finite(sigma) || sigma <= 0) return(-Inf)
      u  <- as.numeric((alpha^(-0.5)/sigma) * (Qsqrt %*% eta2))
      ll <- -N_all * log(sigma) + alpha * sum(u) - alpha * sum(exp(u))
      lp <- log_dtr(sigma, a = 0, b = Inf, mean = 0, sd = 10)
      ll + lp
    }
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    ratio <- (log_target(prop) - log_target(sig2e2)) + (log_q_backward - log_q_forward)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[s] <- sig2e2
    
    ## Fit criterion
    ll[s] <- -2 * sum(dnorm(Y, mean = theta[obs], sd = sqrt(sig2_obs), log = TRUE))  
    if(progress) setTxtProgressBar(pb, s)
  }
  if(progress) close(pb)
  
  keep <- (burn + 1):iter
  theta_bar_obs <- rowMeans(thetaOut[obs, keep, drop = FALSE])         
  sig2_bar_obs  <- rowMeans(sig2Out[obs, keep, drop = FALSE], na.rm = TRUE)
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(Y, mean = theta_bar_obs, sd = sqrt(sig2_bar_obs), log = TRUE))
  
  list(
    Y_hat         = Y,
    obs_idx       = obs,                                               
    Preds         = thetaOut,                                          
    DIC           = DIC,
    sig2          = sig2Out,                                           
    beta1         = beta1Out,
    sig2e1        = sig2e1Out,
    beta2         = beta2Out,
    eta2          = eta2Out,
    sig2e2        = sig2e2Out,
    W             = W,
    E_inc         = E_inc,
    D_var         = D_var,
    Q             = Q,
    Q_var         = Q_var,
    di            = di,
    rankQ         = rankQ,
    ridge_var     = ridge_var,
    N_all         = N_all,
    n_obs         = n_obs,
    state_id      = if(constrained) state_id else NULL,                
    state_levels_obs = if(constrained) state_levels_obs else NULL,     
    state_index_obs  = if(constrained) state_index_obs else NULL,      
    aTot_state    = if(constrained) aTot_lookup else NULL              
  )
}
