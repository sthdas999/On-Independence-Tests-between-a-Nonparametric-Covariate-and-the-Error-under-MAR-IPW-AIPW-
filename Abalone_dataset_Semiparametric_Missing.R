############################################################
## Real Data Study: Abalone Dataset (MAR, CC, IPW, AIPW)
## Author: Sthitadhi Das
############################################################

rm(list = ls())
set.seed(123)

############################################################
## 1. Load required libraries
############################################################
library(tidyverse)
library(np)           # nonparametric regression
library(MASS)
library(energy)       # distance covariance
library(Hmisc)

############################################################
## 2. Load Abalone dataset
############################################################
# Download from Kaggle and place abalone.csv in working dir
abalone <- read.csv("abalone.csv")

# Keep only continuous variables
abalone <- abalone %>%
  dplyr::select(Rings, Length, Diameter, Height, Whole.weight)

# Rename for notation consistency
colnames(abalone) <- c("Y", "V1", "V2", "V3", "V4")

n <- nrow(abalone)

############################################################
## 3. Construct MAR missingness in response
############################################################
alpha0 <- -0.5
alpha1 <- 2.0
alpha2 <- -1.5

linpred <- alpha0 + alpha1 * scale(abalone$V1) +
  alpha2 * scale(abalone$V4)

pi_true <- plogis(linpred)
delta <- rbinom(n, 1, pi_true)

abalone$delta <- delta
abalone$Y_obs <- ifelse(delta == 1, abalone$Y, NA)

mean(delta == 0)   # missing rate ~ 20%

############################################################
## 4. Semiparametric regression function
############################################################
np_fit <- function(Y, V){
  bw <- npregbw(
    xdat = V,
    ydat = Y,
    regtype = "lc",
    bwmethod = "cv.ls"
  )
  fit <- npreg(bw)
  list(fit = fit, fitted = fitted(fit))
}

############################################################
## 5. Complete Case (CC) estimation
############################################################
cc_data <- abalone %>% filter(delta == 1)

fit_cc <- np_fit(cc_data$Y, cc_data[, c("V1","V2","V3","V4")])

beta0_cc <- mean(cc_data$Y - fit_cc$fitted)
res_cc <- cc_data$Y - beta0_cc - fit_cc$fitted

############################################################
## 6. IPW estimation
############################################################
# Estimate propensity score
ps_model <- glm(delta ~ V1 + V4, family = binomial, data = abalone)
pi_hat <- fitted(ps_model)

weights_ipw <- abalone$delta / pi_hat

fit_ipw <- np_fit(
  Y = abalone$Y_obs,
  V = abalone[, c("V1","V2","V3","V4")]
)

beta0_ipw <- weighted.mean(
  abalone$Y_obs - fit_ipw$fitted,
  weights = weights_ipw,
  na.rm = TRUE
)

res_ipw <- abalone$Y_obs - beta0_ipw - fit_ipw$fitted

############################################################
## 7. AIPW estimation
############################################################
# Outcome regression
fit_or <- np_fit(
  Y = cc_data$Y,
  V = cc_data[, c("V1","V2","V3","V4")]
)

m_hat <- predict(
  fit_or$fit,
  newdata = abalone[, c("V1","V2","V3","V4")]
)

aipw_term <- abalone$delta / pi_hat * (abalone$Y_obs - m_hat) + m_hat

beta0_aipw <- mean(aipw_term, na.rm = TRUE)

res_aipw <- aipw_term - beta0_aipw

############################################################
## 8. Dependence test statistics
############################################################

# Rank-based statistic
rank_test <- function(res, V){
  r <- rank(res, na.last = "keep")
  s <- apply(V, 2, rank)
  cor(r, rowMeans(s), use = "complete.obs")
}

# Distance covariance
dcov_test <- function(res, V){
  energy::dcov(res, V)
}

# Spacing-based statistic (1D projection)
spacing_test <- function(res, V){
  z <- scale(V) %*% rep(1/ncol(V), ncol(V))
  o <- order(z)
  mean(diff(res[o])^2, na.rm = TRUE)
}

############################################################
## 9. Permutation test
############################################################
perm_test <- function(stat_fun, res, V, B = 1000){
  obs <- stat_fun(res, V)
  perm <- replicate(B, {
    stat_fun(sample(res), V)
  })
  mean(abs(perm) >= abs(obs))
}

############################################################
## 10. Compute p-values
############################################################
V_cc <- as.matrix(cc_data[, c("V1","V2","V3","V4")])
V_all <- as.matrix(abalone[, c("V1","V2","V3","V4")])

results <- tibble(
  Method = c("CC","IPW","AIPW"),
  Rank_p = c(
    perm_test(rank_test, res_cc, V_cc),
    perm_test(rank_test, res_ipw, V_all),
    perm_test(rank_test, res_aipw, V_all)
  ),
  DCOV_p = c(
    perm_test(dcov_test, res_cc, V_cc),
    perm_test(dcov_test, res_ipw, V_all),
    perm_test(dcov_test, res_aipw, V_all)
  ),
  Spacing_p = c(
    perm_test(spacing_test, res_cc, V_cc),
    perm_test(spacing_test, res_ipw, V_all),
    perm_test(spacing_test, res_aipw, V_all)
  )
)

############################################################
## 11. Final outputs
############################################################

cat("\n--- Intercept Estimates ---\n")
cat("CC   :", round(beta0_cc, 3), "\n")
cat("IPW  :", round(beta0_ipw, 3), "\n")
cat("AIPW :", round(beta0_aipw, 3), "\n")

cat("\n--- Permutation p-values ---\n")
print(results)
