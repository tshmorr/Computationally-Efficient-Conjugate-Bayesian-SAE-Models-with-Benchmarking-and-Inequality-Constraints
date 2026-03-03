rCMLG <- function(H=matrix(rnorm(6),3), alpha=c(1,1,1), kappa=c(1,1,1)){
  ## This function will simulate from the cMLG distribution
  notValid <- TRUE
  while(notValid){
    m <- length(kappa)
    w <- log(rgamma(m, shape=alpha, rate=kappa))
    HtH <- t(H) %*% H
    if(Matrix::isDiagonal(HtH)) HtH <- Matrix::Diagonal(x = diag(HtH))
    invHtH <- solve(HtH)
    tHw <- t(H) %*% w
    cand <- as.numeric(invHtH %*% tHw)
    if(!any(is.na(cand))) notValid <- FALSE
  }
  return(cand)
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
  ## This function simulates from a cMLG distribution truncated from below
  repeat {
    emp <- rCMLG(H=H, alpha=alpha, kappa=kappa)
    if(all(emp > cut)) break
  }
  return(emp)
}

## ============================================================
##   Implementing nonspatial model with constraints
## ============================================================
saeFit_con <- function(
    Y, X, Psi, S2, SS,
    a = 0.5, b = 0.5,
    w = 1000, rho = 1000,
    sig2b = 1000,
    iter = 3000, burn = 1000,
    alpha = 1,
    propSD = 10,
    verbose = FALSE, verbose_iter = 1, verbose_every = 1,
    cVec = NULL, aTot = NULL
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
  
  constrained <- !(is.null(cVec) || is.null(aTot))
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, n)
    if(length(cVec) != n) stop("cVec must be length n (or scalar).")
    if(!is.finite(aTot)) stop("aTot must be a finite scalar.")
    if(sum(cVec) > aTot) stop("Infeasible constraints: sum(cVec) > aTot.")
  }
  
  ## initial values ----
  beta1 <- rep(0, p)
  sig2e1 <- median(S2)  
  beta2 <- rep(0, p)
  eta2  <- rep(0, r)
  sig2e2 <- 10
  sig2 <- as.numeric(S2)
  
  if(constrained){
    theta <- cVec + (aTot - sum(cVec)) / n
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
    do_print <- isTRUE(verbose) && (i <= verbose_iter || i %% verbose_every == 0)
    if (do_print) {
      cat("\n====================\n")
      cat("iter:", i, " time:", format(Sys.time(), "%H:%M:%S"), "\n")
      flush.console()
    }
    
    ## ============================================================
    ## sample var function coefs
    ## ============================================================
    H <- rbind(XP, XP, Matrix::bdiag(alpha^(-0.5)/sqrt(sig2b) * Matrix::Diagonal(p), alpha^(-0.5)/sqrt(sig2e2) * Matrix::Diagonal(r)))
    aVec <- c(rep(1/2, n), rep(df/2, n), rep(alpha, p + r))
    kVec <- c((Y - theta)^2 / 2, 0.5 * (SS - 1) * S2, rep(alpha, p + r))
    
    if (do_print) {
      cat("  B1: start rCMLG\n")
      cat("      nrow(H)=", nrow(H),
          " len(a)=", length(aVec),
          " len(k)=", length(kVec), "\n")
      cat("      k summary:",
          "min=", signif(min(kVec), 3),
          "med=", signif(median(kVec), 3),
          "max=", signif(max(kVec), 3), "\n")
      flush.console()
    }
    
    tmp <- rCMLG(H, aVec, kVec)
    beta2 <- beta2Out[, i] <- tmp[1:p]
    eta2  <- eta2Out[, i]  <- tmp[-c(1:p)]
    
    g <- as.numeric(XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <- log(S2) + g
    sig2 <- sig2Out[, i] <- exp(log_sig2)
    
    if (do_print) {
      cat("  B2: end rCMLG\n"); flush.console()
    }
    
    ## ============================================================
    ## Update theta 
    ## ============================================================
    if (do_print) {
      cat("  A1: start theta-update\n"); flush.console()
    }
    
    mu_theta <- as.numeric(X %*% beta1)
    
    if(!constrained){
      v <- 1 / (1/sig2 + 1/sig2e1)
      m <- v * (Y/sig2 + mu_theta/sig2e1)
      theta <- rnorm(n, mean = m, sd = sqrt(v))
    } else {
      dth <- n 
      
      v_dth <- 1 / (1/sig2[dth] + 1/sig2e1)
      m_dth <- v_dth * (Y[dth]/sig2[dth] + mu_theta[dth]/sig2e1)
      
      max_reject <- 2000
      accept <- FALSE
      
      for(rr in 1:max_reject){
        total_first <- sum(theta[1:(dth-1)])
        bad <- FALSE
        
        for(ii in 1:(dth-1)){
          v_i <- 1 / (1/sig2[ii] + 1/sig2e1)
          m_i <- v_i * (Y[ii]/sig2[ii] + mu_theta[ii]/sig2e1)
          
          S_i <- total_first - theta[ii]
          
          v_tilde <- 1 / (1/v_i + 1/v_dth)
          m_tilde <- v_tilde * (m_i / v_i + (aTot - S_i - m_dth) / v_dth)
          
          U_i <- aTot - cVec[dth] - S_i
          if(!is.finite(U_i) || U_i <= cVec[ii]){
            bad <- TRUE
            break
          }
          
          theta_new <- truncnorm::rtruncnorm(
            1, a = cVec[ii], b = U_i,
            mean = m_tilde, sd = sqrt(v_tilde)
          )
          
          total_first <- S_i + theta_new
          theta[ii] <- theta_new
        }
        
        if(bad) next
        
        theta[dth] <- aTot - sum(theta[1:(dth-1)])
        if(theta[dth] < cVec[dth]) next
        
        accept <- TRUE
        break
      }
      
      if(!accept) stop("Constrained theta-update failed to find a feasible draw. Check cVec/aTot.")
    }
    
    # Optional safety check (kept minimal but helpful)
    if(constrained){
      if(abs(sum(theta) - aTot) > 1e-8) stop("sum(theta) != aTot")
      if(any(theta < cVec - 1e-10)) stop("theta violates lower bounds")
    }
    
    thetaOut[, i] <- theta
    
    if (do_print) {
      cat("  A2: end theta-update\n"); flush.console()
    }
    
    ## ============================================================
    ##  Update beta_1 | theta, sig2e1  
    ## ============================================================
    Qb1 <- (crossprod(X) / sig2e1) + diag(1/sig2b, p)
    Qb1 <- (Qb1 + t(Qb1)) / 2
    Vb1 <- solve(Qb1)
    mb1 <- Vb1 %*% (crossprod(X, theta) / sig2e1)
    
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, i] <- beta1
    
    ## ============================================================
    ##  Update sigma_theta^2 
    ## ============================================================
    shape <- a + n/2
    rate  <- b + 0.5 * sum((theta - as.numeric(X %*% beta1))^2)
    sig2e1 <- 1 / rgamma(1, shape = shape, rate = rate)
    sig2e1Out[i] <- sig2e1
    
    ## ============================================================
    ## sample RE variances
    ## ============================================================
    if (do_print) {
      cat("  C1: start MH sig2e2\n"); flush.console()
    }
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    log_dtr <- function(x, a, b, mean = 0, sd = 1){
      log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))
    }
    
    log_target_prop <- dlMLG(q = eta2,
                             V = sqrt(alpha) * sqrt(prop) * diag(r),
                             a = rep(alpha, r),
                             k = rep(alpha, r))
    
    log_target_curr <- dlMLG(q = eta2,
                             V = sqrt(alpha) * sqrt(sig2e2) * diag(r),
                             a = rep(alpha, r),
                             k = rep(alpha, r))
    log_prior_prop <- log_dtr(sqrt(prop), a = 0, b = Inf, mean = 0, sd = 10)
    log_prior_curr <- log_dtr(sqrt(sig2e2), a = 0, b = Inf, mean = 0, sd = 10)
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    
    ratio <- (log_target_prop - log_target_curr) +
      (log_prior_prop  - log_prior_curr) +
      (log_q_backward  - log_q_forward)
    
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    
    if (do_print) {
      cat("  C2: end MH sig2e2\n"); flush.console()
    }
    
    ll[i] <- -2 * sum(dnorm(Y, mean = theta, sd = sqrt(sig2), log = TRUE))
    setTxtProgressBar(pb, i)
  }
  
  keep <- (burn + 1):iter
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(
    Y,
    mean = rowMeans(thetaOut[, keep, drop = FALSE]),
    sd   = sqrt(rowMeans(sig2Out[, keep, drop = FALSE])),
    log  = TRUE
  ))
  
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




## ==================================================================
##   Implementing spatial dependency on the mean with constraints
## ==================================================================
rCMLG <- function(H, alpha, kappa, tol = 1e-12){
  ## This function will simulate from the cMLG distribution
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
  ## This function simulates from a cMLG distribution truncated from below
  repeat {
    emp <- rCMLG(H=H, alpha=alpha, kappa=kappa)
    if(all(emp > cut)) break
  }
  return(emp)
}

saeFit_con_icar <- function(
    Y, X, Psi, S2, SS,
    W = NULL,                 
    a = 0.5, b = 0.5,
    w = 1000, rho = 1000,
    sig2b = 1000,
    iter = 3000, burn = 1000,
    alpha = 1,
    propSD = 10,
    cVec = NULL, aTot = NULL
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
  
  
  if (is.null(W)) stop("Supply adjacency matrix W (n x n) for spatial model.")
  W <- Matrix::Matrix(W, sparse = TRUE)
  W <- (W + Matrix::t(W)) / 2
  Matrix::diag(W) <- 0
  
  di <- Matrix::rowSums(W)                    # d_i
  Q  <- Matrix::Diagonal(x = di) - W          # Q = D - W
  nb <- lapply(1:n, function(i) which(W[i, ] != 0))
  
  
  count_components <- function(nb){
    n0 <- length(nb)
    visited <- rep(FALSE, n0)
    comp <- 0
    for (i in 1:n0){
      if (!visited[i]){
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
  rankQ <- n - nComp   
  
  # ============================================================
  # Constraints 
  # ============================================================
  constrained <- !(is.null(cVec) || is.null(aTot))
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, n)
    if(length(cVec) != n) stop("cVec must be length n (or scalar).")
    if(!is.finite(aTot)) stop("aTot must be a finite scalar.")
    if(sum(cVec) > aTot) stop("Infeasible constraints: sum(cVec) > aTot.")
  }
  
  ## initial values ----
  beta1 <- rep(0, p)
  sig2e1 <- median(S2)
  beta2 <- rep(0, p)
  eta2  <- rep(0, r)
  sig2e2 <- 10
  sig2 <- as.numeric(S2)
  
  
  if(constrained){
    theta <- cVec + (aTot - sum(cVec)) / n
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
    
    ## ============================================================
    ## sample var function coefs
    ## ============================================================
    H <- rbind(
      XP, XP,
      Matrix::bdiag(
        alpha^(-0.5)/sqrt(sig2b)  * Matrix::Diagonal(p),
        alpha^(-0.5)/sqrt(sig2e2) * Matrix::Diagonal(r)
      )
    )
    
    aVec <- c(rep(1/2, n), rep(df/2, n), rep(alpha, p + r))
    kVec <- c((Y - theta)^2 / 2, 0.5 * (SS - 1) * S2, rep(alpha, p + r))
    tmp <- rCMLG(H, aVec, kVec)
    beta2 <- beta2Out[, i] <- tmp[1:p]
    eta2  <- eta2Out[, i]  <- tmp[-c(1:p)]
    
    ## Update sigma^2 (your variance model)
    g <- as.numeric(XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <- log(S2) + g
    sig2 <- sig2Out[, i] <- exp(log_sig2)
    
    ## ============================================================
    ##  Update theta 
    ## ============================================================
    xb <- as.numeric(X %*% beta1)      # x_i^T beta1
    resid <- theta - xb                # eta1 = theta - X beta1
    
    if(!constrained){
      
      for(ii in 1:n){
        nbs <- nb[[ii]]
        sum_nb <- if(length(nbs)) sum(resid[nbs]) else 0
        
        prec <- (1/sig2[ii]) + (di[ii]/sig2e1)
        v_i  <- 1/prec
        
        m_i <- v_i * ( Y[ii]/sig2[ii] + (di[ii]*xb[ii] + sum_nb)/sig2e1 )
        
        theta[ii] <- rnorm(1, mean = m_i, sd = sqrt(v_i))
        resid[ii] <- theta[ii] - xb[ii]
      }
      
    } else {
      
      dth <- n  
      max_reject <- 2000
      accept <- FALSE
      
      for(rr in 1:max_reject){
        
        
        total_first <- sum(theta[1:(dth-1)])
        resid <- theta - xb
        bad <- FALSE
        
        for(ii in 1:(dth-1)){
          
          
          theta_d_tmp <- aTot - total_first
          resid_d_tmp <- theta_d_tmp - xb[dth]
          
          
          nbs_d <- nb[[dth]]
          sum_nb_d <- if(length(nbs_d)) sum(resid[nbs_d]) else 0
          v_dth <- 1 / (1/sig2[dth] + di[dth]/sig2e1)
          m_dth <- v_dth * ( Y[dth]/sig2[dth] + (di[dth]*xb[dth] + sum_nb_d)/sig2e1 )
          
          
          nbs <- nb[[ii]]
          sum_nb <- if(length(nbs)) sum(resid[nbs]) else 0
          v_i <- 1 / (1/sig2[ii] + di[ii]/sig2e1)
          m_i <- v_i * ( Y[ii]/sig2[ii] + (di[ii]*xb[ii] + sum_nb)/sig2e1 )
          
          
          S_i <- total_first - theta[ii]
          
          
          v_tilde <- 1 / (1/v_i + 1/v_dth)
          m_tilde <- v_tilde * (m_i / v_i + (aTot - S_i - m_dth) / v_dth)
          
          
          U_i <- aTot - cVec[dth] - S_i
          if(!is.finite(U_i) || U_i <= cVec[ii]){
            bad <- TRUE
            break
          }
          
          theta_new <- truncnorm::rtruncnorm(
            1, a = cVec[ii], b = U_i,
            mean = m_tilde, sd = sqrt(v_tilde)
          )
          
          theta[ii] <- theta_new
          resid[ii] <- theta_new - xb[ii]
          total_first <- S_i + theta_new
        }
        
        if(bad) next
        
        
        theta[dth] <- aTot - sum(theta[1:(dth-1)])
        resid[dth] <- theta[dth] - xb[dth]
        
        if(theta[dth] < cVec[dth]) next
        
        accept <- TRUE
        break
      }
      
      if(!accept) stop("Constrained theta-update failed. Check cVec/aTot or adjacency W.")
    }
    
    
    if(constrained){
      if(abs(sum(theta) - aTot) > 1e-8) stop("sum(theta) != aTot")
      if(any(theta < cVec - 1e-10)) stop("theta violates lower bounds")
    }
    
    thetaOut[, i] <- theta
    
    
    ## ============================================================
    ##  Update beta_1 | theta, sig2e1  
    ## ============================================================
    XtQX  <- Matrix::crossprod(X, Q %*% X)        # X^T Q X
    XtQth <- Matrix::crossprod(X, Q %*% theta)    # X^T Q theta
    
    Qb1 <- (XtQX / sig2e1) + diag(1/sig2b, p)
    Qb1 <- (Qb1 + t(Qb1)) / 2
    Vb1 <- solve(Qb1)
    mb1 <- Vb1 %*% (XtQth / sig2e1)
    
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, i] <- beta1
    
    ## ============================================================
    ##  Update sigma^2_{e1} 
    ## ============================================================
    xb <- as.numeric(X %*% beta1)
    resid <- theta - xb
    quad <- as.numeric(Matrix::crossprod(resid, Q %*% resid))  # (theta - Xb)^T Q (theta - Xb)
    
    shape <- a + rankQ/2
    rate  <- b + 0.5 * quad
    sig2e1 <- 1 / rgamma(1, shape = shape, rate = rate)
    sig2e1Out[i] <- sig2e1
    
    ## ============================================================
    ## sample RE variances
    ## ============================================================
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    log_dtr <- function(x, a, b, mean = 0, sd = 1){log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))}
    log_target_prop <- dlMLG(q = eta2,
                             V = sqrt(alpha) * sqrt(prop) * diag(r),
                             a = rep(alpha, r),
                             k = rep(alpha, r))
    
    log_target_curr <- dlMLG(q = eta2,
                             V = sqrt(alpha) * sqrt(sig2e2) * diag(r),
                             a = rep(alpha, r),
                             k = rep(alpha, r))
    
    log_prior_prop <- log_dtr(sqrt(prop), a = 0, b = Inf, mean = 0, sd = 10)
    log_prior_curr <- log_dtr(sqrt(sig2e2), a = 0, b = Inf, mean = 0, sd = 10)
    
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    
    ratio <- (log_target_prop - log_target_curr) +
      (log_prior_prop  - log_prior_curr) +
      (log_q_backward  - log_q_forward)
    
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    
    ll[i] <- -2 * sum(dnorm(Y, mean = theta, sd = sqrt(sig2), log = TRUE))
    setTxtProgressBar(pb, i)
  }
  
  keep <- (burn + 1):iter
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(
    Y,
    mean = rowMeans(thetaOut[, keep, drop = FALSE]),
    sd   = sqrt(rowMeans(sig2Out[, keep, drop = FALSE])),
    log  = TRUE
  ))
  
  list(
    Y_hat = Y,
    Preds = thetaOut,
    DIC = DIC,
    sig2 = sig2Out,
    beta1 = beta1Out,
    sig2e1 = sig2e1Out,
    beta2 = beta2Out,
    eta2 = eta2Out,
    sig2e2 = sig2e2Out,
    W = W,
    Q = Q,
    di = di
  )
}

## ==============================================================================
##   Implementing spatial dependency on the  mean and variance with constraints
## ==============================================================================
rCMLG <- function(H, alpha, kappa, tol = 1e-12){
  ## This function will simulate from the cMLG distribution
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
  ## This function checks whether a matrix is (within a tolerance) the identity matrix
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
  ## This function generates one draw from a multivariate normal distribution with mean mu and covariance Sigma
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

## ============================================================
##   Update theta 
## ============================================================
theta_update_spatial_con_plain <- function(theta, A, B, cVec, aTot){
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
    if(!is.finite(U_i) || U_i <= cVec[i]) stop("Infeasible bounds in constrained theta-update.")
    s <- sum(z)
    t12 <- sum(a12 * z)
    Az_i <- as.numeric(A11[i, ] %*% z)
    
    Qdiag_i <- as.numeric(A11[i, i]) - 2 * a12[i] + a22
    if(!is.finite(Qdiag_i) || Qdiag_i <= 0) stop("Nonpositive conditional precision in theta-update.")
    Qzz_i <- Az_i - t12 - a12[i] * s + a22 * s
    
    sum_off <- Qzz_i - Qdiag_i * z[i]
    
    mu_i <- (qz[i] - sum_off) / Qdiag_i
    sd_i <- sqrt(1 / Qdiag_i)
    
    z[i] <- truncnorm::rtruncnorm(
      1, a = cVec[i], b = U_i,
      mean = mu_i, sd = sd_i
    )
  }
  
  theta_new <- numeric(d)
  theta_new[idx] <- z
  theta_new[d] <- aTot - sum(z)
  
  if(theta_new[d] < cVec[d] - 1e-10) stop("theta_d < c_d after constrained update.")
  if(any(theta_new < cVec - 1e-10)) stop("theta violates lower bounds after constrained update.")
  if(abs(sum(theta_new) - aTot) > 1e-8) stop("sum(theta) != aTot after constrained update.")
  
  theta_new
}

saeFit_con_icar2 <- function(
    Y, X, Psi, S2, SS,
    W,                       
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
  n <- length(Y)
  p <- ncol(X)
  XP <- cbind(X, Psi)
  
  ## constraints
  constrained <- !(is.null(cVec) || is.null(aTot))
  if(constrained){
    cVec <- as.numeric(cVec)
    if(length(cVec) == 1) cVec <- rep(cVec, n)
    if(length(cVec) != n) stop("cVec must be length n (or scalar).")
    if(!is.finite(aTot)) stop("aTot must be a finite scalar.")
    if(sum(cVec) > aTot) stop("Infeasible constraints: sum(cVec) > aTot.")
  }
  
  if(!.is_identity(Psi)) stop("This spatial mean+variance ICAR implementation assumes Psi = I_n.")
  Psi <- Matrix::Diagonal(n)
  
  df_vec <- as.numeric(SS - 1)
  if(length(df_vec) == 1) df_vec <- rep(df_vec, n)
  if(length(df_vec) != n) stop("SS must be scalar or length n.")
  
  S2 <- as.numeric(S2)
  if(length(S2) == 1) S2 <- rep(S2, n)
  if(length(S2) != n) stop("S2 must be scalar or length n.")
  
  W <- Matrix::Matrix(W, sparse = TRUE)
  W <- (W + t(W)) / 2
  Matrix::diag(W) <- 0
  di <- Matrix::rowSums(W)
  Q  <- Matrix::Diagonal(x = di) - W
  Q  <- (Q + t(Q)) / 2
  
  nb <- lapply(1:n, function(i) which(W[i, ] != 0))
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
  rankQ <- n - nComp
  
  Q_var <- Q + Matrix::Diagonal(n, x = ridge_var)
  Qsqrt <- psd_sqrt(Q_var)              
  Qsqrt <- Matrix::Matrix(Qsqrt)       
  
  ## initial values
  beta1  <- rep(0, p)
  sig2e1 <- median(S2)                 
  
  beta2  <- rep(0, p)
  eta2   <- rep(0, n)
  sig2e2 <- 10                          
  
  if(constrained){
    theta <- cVec + (aTot - sum(cVec)) / n
    theta[n] <- aTot - sum(theta[1:(n-1)])
    if(theta[n] < cVec[n]) stop("Initial theta not feasible.")
  } else {
    theta <- as.numeric(Y)
  }
  
  ## storage
  thetaOut  <- matrix(NA_real_, n, iter)
  beta1Out  <- matrix(NA_real_, p, iter)
  sig2e1Out <- numeric(iter)
  
  beta2Out  <- matrix(NA_real_, p, iter)
  eta2Out   <- matrix(NA_real_, n, iter)
  sig2Out   <- matrix(NA_real_, n, iter)
  sig2e2Out <- numeric(iter)
  
  ll <- numeric(iter)
  
  pb <- NULL
  if(isTRUE(progress)){
    pb <- txtProgressBar(min = 0, max = iter, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  eps_k <- 1e-12
  
  QX <- Q %*% X
  XtQX <- Matrix::crossprod(X, QX)
  
  for(s in 1:iter){
    
    ## ============================================================
    ## sample RE variances
    ## ============================================================
    expPsiEta <- exp(as.numeric(Psi %*% eta2))  # = exp(eta2)
    
    H_beta2 <- rbind(
      X,
      X,
      alpha^(-0.5)/sqrt(sig2b2) * Matrix::Diagonal(p)
    )
    
    a_beta2 <- c(rep(1/2, n), df_vec/2, rep(alpha, p))
    k_beta2 <- c(
      0.5 * (Y - theta)^2 * expPsiEta,
      (df_vec/2) * S2 * expPsiEta,
      rep(alpha, p)
    )
    k_beta2 <- pmax(k_beta2, eps_k)
    
    beta2 <- rCMLG(H_beta2, a_beta2, k_beta2)
    beta2Out[, s] <- beta2
    
    
    expXBeta2 <- exp(as.numeric(X %*% beta2))
    
    H_eta2 <- rbind(
      Psi,
      Psi,
      alpha^(-0.5)/sqrt(sig2e2) * Qsqrt
    )
    
    a_eta2 <- c(rep(1/2, n), df_vec/2, rep(alpha, n))
    k_eta2 <- c(
      0.5 * (Y - theta)^2 * expXBeta2,
      (df_vec/2) * S2 * expXBeta2,
      rep(alpha, n)
    )
    k_eta2 <- pmax(k_eta2, eps_k)
    
    eta2 <- rCMLG(H_eta2, a_eta2, k_eta2)
    eta2Out[, s] <- eta2
    
    g <- as.numeric(XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <- log(S2) + g
    sig2 <- sig2Out[, s] <- exp(log_sig2)
    
    ## ============================================================
    ## Update theta 
    ## ============================================================
    A <- Matrix::Diagonal(x = 1/sig2) + (1/sig2e1) * Q
    B <- (Y / sig2) + (1/sig2e1) * as.numeric(Q %*% (X %*% beta1))
    
    if(!constrained){
      
      A_dense <- as.matrix(A)
      R <- chol((A_dense + t(A_dense))/2)
      mu <- backsolve(R, forwardsolve(t(R), as.numeric(B)))
      z <- rnorm(n)
      v <- backsolve(R, z)
      theta <- as.numeric(mu + v)
    } else {
      theta <- theta_update_spatial_con_plain(theta, A, B, cVec, aTot)
    }
    thetaOut[, s] <- theta
    
    ## ============================================================
    ## Update beta_1 | theta, sig2e1 
    ## ============================================================
    XtQtheta <- Matrix::crossprod(QX, theta)
    
    Qb1 <- (XtQX / sig2e1) + diag(1/sig2b1, p)
    Qb1 <- (Qb1 + t(Qb1))/2
    Vb1 <- solve(as.matrix(Qb1))
    mb1 <- Vb1 %*% (as.numeric(XtQtheta) / sig2e1)
    
    beta1 <- rmvn1(mb1, Vb1)
    beta1Out[, s] <- beta1
    
    ## ============================================================
    ## Update sigma_theta^2 
    ## ============================================================
    resid1 <- theta - as.numeric(X %*% beta1)
    quad1 <- as.numeric(Matrix::crossprod(resid1, Q %*% resid1))
    
    shape1 <- a + rankQ/2
    rate1  <- b + 0.5 * quad1
    sig2e1 <- 1 / rgamma(1, shape = shape1, rate = rate1)
    sig2e1Out[s] <- sig2e1
    
    ## ============================================================
    ## sample RE variances
    ## ============================================================
    prop <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = sig2e2, sd = propSD)
    
    log_dtr <- function(x, a, b, mean = 0, sd = 1){
      log(truncnorm::dtruncnorm(x, a = a, b = b, mean = mean, sd = sd))
    }
    
    log_target <- function(sig2_val){
      sigma <- sqrt(sig2_val)
      if(!is.finite(sigma) || sigma <= 0) return(-Inf)
      
      u <- as.numeric((alpha^(-0.5)/sigma) * (Qsqrt %*% eta2))
      ll <- -n * log(sigma) + alpha * sum(u) - alpha * sum(exp(u))
      lp <- log_dtr(sigma, a = 0, b = Inf, mean = 0, sd = 10)
      ll + lp
    }
    
    log_q_backward <- log_dtr(sig2e2, a = 0, b = Inf, mean = prop,  sd = propSD)
    log_q_forward  <- log_dtr(prop,  a = 0, b = Inf, mean = sig2e2, sd = propSD)
    
    ratio <- (log_target(prop) - log_target(sig2e2)) + (log_q_backward - log_q_forward)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[s] <- sig2e2
    
    ll[s] <- -2 * sum(dnorm(Y, mean = theta, sd = sqrt(sig2), log = TRUE))
    if(!is.null(pb)) setTxtProgressBar(pb, s)
  }
  
  keep <- (burn + 1):iter
  DIC <- 2 * mean(ll[keep]) + 2 * sum(dnorm(
    Y,
    mean = rowMeans(thetaOut[, keep, drop = FALSE]),
    sd   = sqrt(rowMeans(sig2Out[, keep, drop = FALSE])),
    log  = TRUE
  ))
  
  list(
    Y_hat  = Y,
    Preds  = thetaOut,
    DIC    = DIC,
    sig2   = sig2Out,
    beta1  = beta1Out,
    sig2e1 = sig2e1Out,
    beta2  = beta2Out,
    eta2   = eta2Out,
    sig2e2 = sig2e2Out,
    W = W,
    Q = Q,
    di = di,
    rankQ = rankQ,
    ridge_var = ridge_var
  )
}
