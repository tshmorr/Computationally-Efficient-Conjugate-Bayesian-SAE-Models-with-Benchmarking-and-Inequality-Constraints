
library(tidycensus)
library(tigris)
library(spdep)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(Matrix)
library(mvtnorm)
library(readr)
library(sampling)
library(mase)
library(truncnorm)
library(statmod)
library(matrixStats)
library(VGAM)
library(rstan)
library(expm)
library(parallel)
library(MASS)
library(spatialreg)
library(stringr)
library(purrr)
library(Rcpp)
library(spdep)
source("MCMC.R")

emp <-read_rds("empdata.rds")
ia <- tigris::counties(state = "Iowa", cb = TRUE, class = "sf")
ia$NAME <- tolower(trimws(gsub("\\s*county\\b", "", ia$NAME, ignore.case = TRUE)))
ia <- ia %>% dplyr::filter(tolower(NAME) %in% tolower(emp$NAME)) %>% dplyr::slice(match(tolower(emp$NAME), tolower(NAME)))
W <-as.matrix(as_dgRMatrix_listw(nb2listw(poly2nb(ia), style='B')))
XmatA <- cbind(1, emp$amtoffertilizer)
aTot <- sum(emp$theta_true)

fit_all_models <- function(sim_i){
  fits  <- list()
  times <- numeric(2)
  names(times) <- c("con","icar2")
  times["con"] <- system.time(fits$con <- saeFit_con(Y=sim_i$Y_hat, X=XmatA, Psi=diag(length(sim_i$theta_true)), S2=sim_i$Var, SS=sim_i$N, iter=3000, burn=1000, cVec=sim_i$ci, aTot=aTot))["elapsed"]
  times["icar2"] <- system.time(fits$icar2 <- saeFit_con_icar2(Y=sim_i$Y_hat, X=XmatA, Psi=diag(length(sim_i$theta_true)), S2=sim_i$Var, SS=sim_i$N, W=W, iter=3000, burn=1000, cVec=sim_i$ci, aTot=aTot))["elapsed"]
  list(fits = fits, time = times)
}

sim_ynoisy <- function(empdata) {
  empdata$Y_hat <- empdata$theta_true + rnorm(n = nrow(empdata), mean = 0, sd = sqrt(empdata$Var))
  empdata
}

run_models <- function(base_data) {
  sim_i <- sim_ynoisy(base_data)
  out <- fit_all_models(sim_i)
  list(fits = out$fits, time = out$time)
}

empirical_sim <- function(base_data, R = 100) {
  pb <- utils::txtProgressBar(min = 0, max = R, style = 3); on.exit(close(pb), add = TRUE)
  res <- vector("list", R)
  for (i in seq_len(R)) {
    message(sprintf("Rep %d / %d", i, R))
    res[[i]] <- run_models(base_data)
    utils::setTxtProgressBar(pb, i)}
  names(res) <- paste0("rep", seq_len(R))
  cat("\n")
  list(results = res)
}

set.seed(3)
emp_sim_results<- empirical_sim(emp, R = 100)

# RMSE
rmse <- function(est, true) {sqrt(mean((est - true)^2))}
n_reps <- length(emp_sim_results$results)
halm_rmse  <- numeric(n_reps)
icar2_rmse <- numeric(n_reps)
direct_rmse <- numeric(n_reps)
preds_con  <- vector("list", n_reps)
preds_icar2 <- vector("list", n_reps)

for (t in seq_len(n_reps)) {
  preds_con[[t]]  <- rowMeans(emp_sim_results$results[[t]]$fits$con$Preds[, -c(1:1000)])
  halm_rmse[t]    <- rmse(preds_con[[t]], emp$theta_true)
  preds_icar2[[t]] <- rowMeans(emp_sim_results$results[[t]]$fits$icar2$Preds[, -c(1:1000)])
  icar2_rmse[t]    <- rmse(preds_icar2[[t]], emp$theta_true)
  direct_rmse[t]  <- rmse(emp_sim_results$results[[t]]$fits$con$Y_hat, emp$theta_true)
}
c(halm  = mean(halm_rmse), icar2 = mean(icar2_rmse), direct = mean(direct_rmse))


# COVERAGE & INTERVAL SCORE
coverage_rate <- function(l, u, theta) {mean(theta >= l & theta <= u)}
interval_score <- function(l, u, theta, alpha = 0.05) {
  width <- u - l; under <- (2/alpha) * (l - theta) * (theta < l); over <- (2/alpha) * (theta - u) * (theta > u)
  width + under + over
}

n_reps <- length(emp_sim_results$results)
cover_con  <- numeric(n_reps); is_con  <- numeric(n_reps)
cover_icar2 <- numeric(n_reps); is_icar2 <- numeric(n_reps)
cover_direct <- numeric(n_reps); is_direct <- numeric(n_reps)
sd_direct <- sqrt(as.numeric(unlist(emp$Var)))

for (t in seq_len(n_reps)) {
  draws_con <- emp_sim_results$results[[t]]$fits$con$Preds[, -c(1:1000)]
  low_con   <- apply(draws_con, 1, quantile, probs = 0.025)
  high_con  <- apply(draws_con, 1, quantile, probs = 0.975)
  cover_con[t] <- coverage_rate(low_con, high_con, emp$theta_true)
  is_con[t]    <- mean(interval_score(low_con, high_con, emp$theta_true))
  draws_icar2 <- emp_sim_results$results[[t]]$fits$icar2$Preds[, -c(1:1000)]
  low_icar2   <- apply(draws_icar2, 1, quantile, probs = 0.025)
  high_icar2  <- apply(draws_icar2, 1, quantile, probs = 0.975)
  cover_icar2[t] <- coverage_rate(low_icar2, high_icar2, emp$theta_true)
  is_icar2[t]    <- mean(interval_score(low_icar2, high_icar2, emp$theta_true))
  yhat <- emp_sim_results$results[[t]]$fits$con$Y_hat
  low_direct  <- yhat - 1.96 * sd_direct
  high_direct <- yhat + 1.96 * sd_direct
  cover_direct[t] <- coverage_rate(low_direct, high_direct, emp$theta_true)
  is_direct[t]    <- mean(interval_score(low_direct, high_direct, emp$theta_true))
}

coverage_means <- c(con = mean(cover_con), icar2 = mean(cover_icar2), direct = mean(cover_direct))
intervalscore_means <- c(con = mean(is_con), icar2 = mean(is_icar2), direct = mean(is_direct))
coverage_means
intervalscore_means

# RUNTIME
con_time <- sapply(emp_sim_results$results, function(x) x$time[1])
icar2_time <- sapply(emp_sim_results$results, function(x) x$time[2])
runtimes <- c(halm_time = mean(con_time), hicar2_time = mean(icar2_time))
runtimes


