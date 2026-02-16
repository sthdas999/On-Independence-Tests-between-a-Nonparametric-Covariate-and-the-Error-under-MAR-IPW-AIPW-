
############################################################
# Independence Testing in Partially Linear Model under MAR
# Author: Sthitadhi Das
# Description: CC, IPW, AIPW estimation + residual-based tests
############################################################

# Required libraries
library(np)        # Nonparametric regression
library(splines)
library(MASS)
library(boot)

set.seed(123)

############################################################
# 1. Load Data
############################################################
# Replace with your actual dataset path
# data <- read.csv("your_data.csv")

# Example structure assumption:
# Y  : response (glucose)
# V1 : age (nonparametric covariate X)
# V2,V3,V4 : parametric covariates

############################################################
# 2. Generate MAR Missingness
############################################################
generate_mar <- function(data, miss_rate = 0.2) {
  X <- data$V1
  V3 <- data$V3
  
  linpred <- -1 + 0.02*X + 0.01*V3
  pi_hat <- 1 / (1 + exp(-linpred))
  
  delta <- rbinom(nrow(data), 1, pi_hat)
  data$delta <- delta
  data$Y_obs <- ifelse(delta == 1, data$Y, NA)
  
  return(data)
}

############################################################
# 3. Propensity Score Estimation
############################################################
estimate_propensity <- function(data) {
  fit <- glm(delta ~ V1 + V2 + V3 + V4,
             family = binomial,
             data = data)
  data$pi_est <- predict(fit, type = "response")
  return(data)
}

############################################################
# 4. Partially Linear Estimation
############################################################
fit_pl_model <- function(data, weights = NULL) {
  
  # Nonparametric part
  np_fit <- npreg(Y_obs ~ V1, data = data, regtype = "ll")
  m_hat <- fitted(np_fit)
  
  # Parametric part (weighted if needed)
  Z <- cbind(1, data$V2, data$V3, data$V4)
  Y_adj <- data$Y_obs - m_hat
  
  if (is.null(weights)) {
    beta_hat <- solve(t(Z) %*% Z) %*% t(Z) %*% Y_adj
  } else {
    W <- diag(weights)
    beta_hat <- solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W %*% Y_adj
  }
  
  residuals <- Y_adj - Z %*% beta_hat
  return(list(beta = beta_hat, residuals = residuals))
}

############################################################
# 5. CC, IPW, AIPW Estimation
############################################################
run_estimators <- function(data) {
  
  # Complete Case
  cc_data <- subset(data, delta == 1)
  cc_fit <- fit_pl_model(cc_data)
  
  # IPW
  ipw_weights <- data$delta / data$pi_est
  ipw_fit <- fit_pl_model(data, weights = ipw_weights)
  
  # AIPW (simple augmentation illustration)
  outcome_model <- lm(Y_obs ~ V1 + V2 + V3 + V4, data = data)
  mu_hat <- predict(outcome_model, newdata = data)
  
  aipw_weights <- data$delta / data$pi_est
  Y_aug <- mu_hat + aipw_weights * (data$Y_obs - mu_hat)
  data$Y_obs <- Y_aug
  
  aipw_fit <- fit_pl_model(data)
  
  return(list(CC = cc_fit,
              IPW = ipw_fit,
              AIPW = aipw_fit))
}

############################################################
# 6. Simple Residual-Based Independence Test (Permutation)
############################################################
perm_test <- function(residuals, X, B = 1000) {
  obs_stat <- cor(residuals, X, method = "spearman")
  perm_stats <- replicate(B, {
    cor(sample(residuals), X, method = "spearman")
  })
  p_value <- mean(abs(perm_stats) >= abs(obs_stat))
  return(p_value)
}

############################################################
# End of Script
############################################################
cat("R script for MAR semiparametric independence testing loaded successfully.\n")
