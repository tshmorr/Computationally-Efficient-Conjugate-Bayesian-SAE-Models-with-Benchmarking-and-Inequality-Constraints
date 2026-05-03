library(readr)
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
library(data.table)
library(readxl)
library(janitor)
library(scales)
library(sf)
library(cowplot)
library(FNN)
source('mcmc_chalm_nonspatial.r')

## -------------------------------------------------
## Read analysis df
## -------------------------------------------------
df <- readRDS("analysis_data.RDS")
ggplot(df) + geom_sf() + theme_minimal()

## Addressing missing data issues
## 1. Lower bound
df$cvec[is.na(df$cvec)] <- 0

## 2. CV
cv_defaults <- c(IL = 1.3, WI = 4.2, OH = 2.2, MO = 1.2, MN = 1.3, MI = 4.3, IA = 3.2, IN = 1.0)
i <- is.na(df$`CV (%)`)
df$`CV (%)`[i] <- cv_defaults[df$STUSPS[i]]

## 3. Production Bushels for other counties
other_counties <- list(IL = 108146367, IN = 987654, IA = 14193900, MI = 16292256, MN = 127868076, MO = 114921840, OH = 9839148, WI = 410337)
temp <- knn_constrained_allocation_by_state(counties_sf = df, value_col = "production_bushels", state_col = "STUSPS", all_other_totals = other_counties, k = 5)
temp <- temp %>% st_drop_geometry() %>% dplyr::select(STATEFP, COUNTYFP, value_imputed, allocation_status)
idx_na <- which(is.na(df$production_bushels))
match_idx <- match(paste(df$STATEFP[idx_na], df$COUNTYFP[idx_na]), paste(temp$STATEFP, temp$COUNTYFP))
df$production_bushels[idx_na] <- temp$value_imputed[match_idx]
df <- df %>% dplyr::select(-se, -Var)

## 4. SE and Variance
df$se <- (df$`CV (%)`)/100 * df$production_bushels
df$Var <- df$se^2

## -------------------------------------------------
## Model Fitting
## -------------------------------------------------
#aTot_state <- c("IL" = 1716079078, "IA"= 2419901737, "MN" = 1324138682, "IN" = 879918936, "MI" = 302808059, "OH" = 584933770, "WI" = 505376067, "MO" = 301078583)
aTot_state <- c("IL" = 2143284486, "IA"= 2511484113, "MN" = 1478451268, "IN" = 1039725085, "MI" = 467281181, "OH" = 608637574, "WI" = 533043125, "MO" = 301078583)
XmatA <- cbind(1, df$Amtoffert)
n <- length(df$production_bushels)
state_id <- df$STUSPS
cVec <- df$cvec

fit <- saeFit_con_multistate(Y = df$production_bushels, X = XmatA, Psi = diag(n), S2 = df$Var, SS = df$SS, a = 0.5, b = 0.5, w = 1000, rho = 1000, sig2b = 1000, iter = 3000, burn = 1000, alpha = 1, propSD = 10, cVec = cVec, aTot_state = aTot_state, state_id = state_id)
boundaries <- tigris::counties(state = c("Iowa","Illinois","Indiana","Michigan","Minnesota","Missouri","Ohio","Wisconsin"), cb = TRUE, class = "sf")
boundaries <- boundaries[, -c(3,4,5,7,9,10,11,12)]
ggplot(boundaries) + geom_sf() + theme_minimal()

boundaries$NAME <- tolower(trimws(gsub("\\s*county\\b", "", boundaries$NAME, ignore.case = TRUE)))
### Plot results
burn <- 1000
keep <- (burn + 1):ncol(fit$Preds)
post_theta <- fit$Preds[, keep, drop = FALSE]

analysis_df <- df %>%mutate(direct_est = production_bushels, direct_se = sqrt(Var), model_est = rowMeans(post_theta), model_se = apply(post_theta, 1, sd), model_lwr = apply(post_theta, 1, quantile, probs = 0.025),model_upr  = apply(post_theta, 1, quantile, probs = 0.975),ratio_model_direct = ifelse(direct_est > 0, model_est / direct_est, NA_real_),diff_model_direct  = model_est - direct_est,benchmark_total = unname(aTot_state[STUSPS]))
plot_df <- analysis_df %>% st_drop_geometry()
state_outline <- boundaries %>%group_by(STUSPS) %>%summarise(.groups = "drop")

## -------------------------------------------------
## Graphs
## -------------------------------------------------
# map of the model-based estimates and the posterior standard errors
p_map_est_one <- ggplot() + geom_sf(data = boundaries, fill = "grey85", color = "black", linewidth = 0.08) + geom_sf(data = analysis_df,aes(fill = model_est),color = "white",linewidth = 0.08) +geom_sf(data = state_outline,fill = NA,color = "black",linewidth = 0.45) +
  scale_fill_viridis_c(limits = c(1000, 75000000), labels = label_number(scale_cut = cut_short_scale()),na.value = "grey85") +labs(title = "Model-based Estimates",fill = "Bushels") +
  theme(legend.position = "right",plot.title = element_text(face = "bold"))

p_map_se_one <- ggplot() + geom_sf(data = boundaries, fill = "grey85", color = "black", linewidth = 0.08) +  geom_sf(data = analysis_df,aes(fill = model_se),color = "white",linewidth = 0.08) +geom_sf(data = state_outline,fill = NA,color = "black",linewidth = 0.45) +
  scale_fill_viridis_c(labels = label_number(scale_cut = cut_short_scale()),na.value = "grey85") +labs(title = "Posterior Standard Errors of Model-based Estimates",fill = "SE") +
 theme(legend.position = "right",plot.title = element_text(face = "bold"))

# map of the direct estimates and the standard errors of the direct estimates
p_map_direct_est <- ggplot() + geom_sf(data = boundaries, fill = "grey85", color = "black", linewidth = 0.08) +  geom_sf(data = analysis_df, aes(fill = direct_est), color = "white", linewidth = 0.08) +geom_sf(data = state_outline, fill = NA, color = "black", linewidth = 0.45) +
  scale_fill_viridis_c(limits = c(1000, 75000000), labels = label_number(scale_cut = cut_short_scale()),na.value = "grey85") +labs(title = "Direct Estimates",fill = "Bushels") +
  theme(legend.position = "right",plot.title = element_text(face = "bold"))

p_map_direct_se <- ggplot() + geom_sf(data = boundaries, fill = "grey85", color = "black", linewidth = 0.08) + geom_sf(data = analysis_df, aes(fill = direct_se), color = "white", linewidth = 0.08) +geom_sf(data = state_outline, fill = NA, color = "black", linewidth = 0.45) +
  scale_fill_viridis_c(labels = label_number(scale_cut = cut_short_scale()),na.value = "grey85") +labs(title = "Standard Error of Direct Estimates",fill = "SE") +
 theme(legend.position = "right",plot.title = element_text(face = "bold"))

p_two_maps2 <- cowplot::plot_grid(p_map_direct_est, p_map_est_one, ncol = 2)#, labels = c("A", "B"), label_size = 12)
#notes_text <- paste("Notes: The model based estimates are based on the non spatial constrained HALM.\n", " Gray areas indicate counties with unpublished corn production estimates", sep = "")

final_plot2 <- ggdraw() + draw_label("Comparing Corn Production Estimates: Direct Estimates Vs. Estimates from the CHALM Model", x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, fontface = "bold", size = 24) +
  draw_plot(p_two_maps2, x = 0.02, y = 0.18, width = 0.96, height = 0.72) #+ draw_label(notes_text, x = 0.02, y = 0.04, hjust = 0, vjust = 0, size = 9)
final_plot2

# Uncertainty plots
p_two_maps3 <- cowplot::plot_grid(p_map_direct_se, p_map_se_one, ncol = 2)#, labels = c("C", "D"), label_size = 12)
#notes_text <- paste("Notes: The posterior standard errors of moodel based estimates are based on the non spatial constrained HALM.\n", " Gray areas indicate counties with unpublished corn production estimates", sep = "")

final_plot3 <- ggdraw() + draw_label("Comparing Corn Production Uncertainty Estimates: Direct Vs. Posterior (from the CHALM Model)", x = 0.5, y = 0.98, hjust = 0.5, vjust = 1, fontface = "bold", size = 20) +
  draw_plot(p_two_maps3, x = 0.02, y = 0.18, width = 0.96, height = 0.72) #+ draw_label(notes_text, x = 0.02, y = 0.04, hjust = 0, vjust = 0, size = 9)
final_plot3


