###############################################################################
# Unconstrained Models
###############################################################################

## Unconstrained non-spatial model
saeFit_uncon <- function(Y, X, Psi, S2, SS, a=0.5, b=0.5, w=1000, rho=1000, sig2b=1000, iter=3000, burn=1000, alpha=1, propSD=10){
  df <- (SS-1)
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- 1
  sig2e2 <- 10
  XP <- cbind(X, Psi)
  ICAR <- Diagonal(n)
  cholICAR <- solve(t(chol(ICAR)))
  aE <- c(rep(alpha, n), w)
  kE <- c(rep(alpha, n), rho)
  
  sig2Out <- matrix(NA, nrow=n, ncol=iter)
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  beta1Out <- beta2Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- eta2Out <- matrix(NA, nrow=r, ncol=iter)
  sig2e2Out <- rep(NA, iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    ## sample mean function coefs
    Q <- t(XP) %*% Diagonal(x = 1/sig2) %*% XP + bdiag((1/sig2b) * Diagonal(p), (1/sig2e1) * Diagonal(r))
    Q  <- (Q + t(Q)) / 2
    eps <- 1e-8 * mean(diag(Q))
    Q <- Q + eps * Diagonal(nrow(Q))
    var <- solve(Q)
    var <- as.matrix(var)
    var <- (var + t(var)) / 2
    mu <- var%*%t(XP)%*%Diagonal(x=1/sig2)%*%Y
    tmp <- as.numeric(rmvnorm(1, mean=mu, sigma=var))
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i]  <- XP%*%c(beta1,eta1)

    ## sample var function coefs
    H <- rbind(XP, XP, bdiag(alpha^(-0.5)/sqrt(sig2b)*Diagonal(p),  alpha^(-0.5)/sqrt(sig2e2)*Diagonal(r)))
    a <- c(rep(1/2,n), df/2, rep(alpha,p+r))
    k <- c((Y-theta)^2/2, 1/2 *(SS-1)*S2, rep(alpha,p+r))
    tmp <- rCMLG(H, a, k)
    beta2 <- beta2Out[,i] <- tmp[1:p]
    eta2 <- eta2Out[,i] <- tmp[-c(1:p)]
    g <- (XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <-  log(S2) + g
    sig2     <- sig2Out[, i] <- exp(log_sig2)

    ## sample RE variances
    sig2e1 <- 1/rgamma(1, a+r/2, b+sum(eta1^2)/2)
    prop <- rtruncnorm(1, 0, mean=sig2e2, sd=propSD)
    ratio <- dlMLG(q=eta2, V=sqrt(alpha)*sqrt(prop)*diag(r), a=rep(alpha,r), k=rep(alpha,r)) -
      dlMLG(eta2, sqrt(alpha)*sqrt(sig2e2)*diag(r), rep(alpha,r), rep(alpha,r)) +
      dtruncnorm(sqrt(prop), 0, sd=10) - dtruncnorm(sqrt(sig2e2), 0, sd=10) +
      dtruncnorm(sig2e2, 0, mean=prop, sd=propSD) - dtruncnorm(prop, 0, mean=sig2e2, sd=propSD)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(sig2), log=T))
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(rowMeans(sig2Out[,-c(1:burn)])), log=T))
  return(list(Preds=thetaOut, DIC=DIC, sig2=sig2Out, beta1=beta1Out, beta2=beta2Out, eta2=eta2Out, sig2e2=sig2e2Out))
}



## Unconstrained spatial ICAR model
saeFit_uncon_icar2 <- function(Y, X, Psi, S2, SS, W, a=0.5, b=0.5, w=1, rho=1, sig2b=1000, iter=3000, burn=1000, alpha=1, propSD=10){
  df <- (SS-1)
  n <- length(Y)
  p <- ncol(X)
  r <- ncol(Psi)
  sig2 <- S2
  sig2e1 <- 1
  sig2e2 <- 10
  XP <- cbind(X, Psi)
  dep <- 0.99
  ICAR <- (diag(rowSums(W)) - dep*W)
  cholICAR <- t(chol(ICAR))
  Vinv <- solve(cholICAR)
  aE <- c(rep(alpha, n), w)
  kE <- c(rep(alpha, n), rho)
  
  depOut <- rep(NA, iter)
  beta1 <- rep(1,p)
  sig2Out <- matrix(NA, nrow=n, ncol=iter)
  thetaOut <- matrix(NA, nrow=n, ncol=iter)
  beta1Out <- beta2Out <- matrix(NA, nrow=p, ncol=iter)
  eta1Out <- eta2Out <- matrix(NA, nrow=r, ncol=iter)
  sig2e2Out <- rep(NA,iter)
  ll <- rep(NA, iter)
  pb <- txtProgressBar(min=0, max=iter, style=3)
  for(i in 1:iter){
    Q <- forceSymmetric(t(XP)%*%Diagonal(x=1/sig2)%*%XP + bdiag(1/sig2b*Diagonal(p), 1/sig2e1*ICAR))
    Q  <- (Q + t(Q)) / 2
    eps <- 1e-8 * mean(diag(Q))
    Q <- Q + eps * Diagonal(nrow(Q))
    U <- chol(Q)
    a <- t(XP)%*%Diagonal(x=1/sig2)%*%Y
    bb <- rnorm(nrow(Q))
    tmp <- backsolve(U, backsolve(U, a, transpose = TRUE) + bb)
    beta1 <- beta1Out[,i] <- tmp[1:p]
    eta1 <- eta1Out[,i] <- tmp[-c(1:p)]
    theta <- thetaOut[,i]  <- XP%*%c(beta1,eta1)
   
    prop <- runif(1)
    ratio <- dmvnorm(eta1, sigma=sig2e1*solve(diag(rowSums(W)) - prop*W), log=T) -
      dmvnorm(eta1, sigma=sig2e1*solve(diag(rowSums(W)) - dep*W), log=T)
    if(runif(1) < exp(ratio)) dep <- prop
    depOut[i] <- dep
    ICAR <- diag(rowSums(W)) - dep*W
    cholICAR <- t(chol(ICAR))
    Vinv <- solve(cholICAR)
    
    ## sample var function coefs
    H <- rbind(XP, XP, bdiag(alpha^(-0.5)/sqrt(sig2b)*Diagonal(p),  alpha^(-0.5)/sqrt(sig2e2)*t(cholICAR)))
    a <- c(rep(1/2,n), df/2, rep(alpha,p+r))
    k <- c((Y-theta)^2/2, 1/2 *(SS-1)*S2, rep(alpha,p+r))
    tmp <- rCMLG(H, a, k)
    beta2 <- beta2Out[,i] <- tmp[1:p]
    eta2 <- eta2Out[,i] <- tmp[-c(1:p)]
    g <- (XP %*% c(beta2, eta2))
    g <- g - mean(g)
    log_sig2 <-  log(S2) + g
    sig2     <- sig2Out[, i] <- exp(log_sig2)
    
    ## sample RE variances
    sig2e1 <- 1/rgamma(1, a+r/2, b+t(eta1)%*%ICAR%*%eta1/2)
    prop <- rtruncnorm(1, 0, mean=sig2e2, sd=propSD)
    ratio <- dlMLG(eta2, sqrt(alpha)*sqrt(prop)*Vinv, rep(alpha,r), rep(alpha,r)) -
      dlMLG(eta2, sqrt(alpha)*sqrt(sig2e2)*Vinv, rep(alpha,r), rep(alpha,r)) +
      dtruncnorm(sqrt(prop), 0, sd=10) - dtruncnorm(sqrt(sig2e2), 0, sd=10) +
      dtruncnorm(sig2e2, 0, mean=prop, sd=propSD) - dtruncnorm(prop, 0, mean=sig2e2, sd=propSD)
    if(runif(1) < exp(ratio)) sig2e2 <- prop
    sig2e2Out[i] <- sig2e2
    
    ll[i] <- -2*sum(dnorm(Y, mean = theta, sd=sqrt(sig2), log=T))
    setTxtProgressBar(pb, i)
  }
  DIC <- 2*mean(ll[-c(1:burn)]) + 2*sum(dnorm(Y, mean=rowMeans(thetaOut[,-c(1:burn)]), sd=sqrt(rowMeans(sig2Out[,-c(1:burn)])), log=T))
  return(list(Preds=thetaOut, DIC=DIC, sig2=sig2Out, beta1=beta1Out, beta2=beta2Out, eta1=eta1Out, eta2=eta2Out, dep=depOut, sig2e2=sig2e2Out))
}
