############################################################
## Real Data Analysis with MAR Missingness
## CC / IPW / AIPW + Permutation Tests
## Fully editable template
############################################################

set.seed(123)

############################################################
## 1. LOAD LIBRARIES
############################################################

library(stats)
library(MASS)
library(dplyr)
library(knitr)

############################################################
## 2. LOAD OR DEFINE DATA
############################################################
## Replace this block with your real blood glucose data

# Example placeholders
n <- 520
V <- cbind(
  age = rnorm(n, 45, 12),
  bmi = rnorm(n, 27, 4),
  bp  = rnorm(n, 80, 10),
  ins = rnorm(n, 15, 6)
)

y <- 120 + 0.4 * V[,2] + 0.3 * V[,4] + rnorm(n, 0, 8)

############################################################
## 3. MAR MISSINGNESS GENERATOR
############################################################

generate_MAR <- function(y, V, miss_rate = 0.2) {
  
  lin_pred <- -1 + 0.02 * V[,1] + 0.01 * V[,3]
  prob_obs <- plogis(lin_pred)
  
  # Scale to target missingness
  prob_obs <- prob_obs * (1 - miss_rate) / mean(prob_obs)
  prob_obs[prob_obs > 1] <- 1
  
  delta <- rbinom(length(y), 1, prob_obs)
  y_obs <- y
  y_obs[delta == 0] <- NA
  
  list(y = y_obs, delta = delta)
}

############################################################
## 4. RESIDUAL CONSTRUCTION
############################################################

get_residuals <- function(y, V, delta, method = c("CC","IPW","AIPW")) {
  
  method <- match.arg(method)
  
  # Propensity score
  ps_fit <- glm(delta ~ V[,1] + V[,3], family = binomial)
  pi_hat <- fitted(ps_fit)
  
  # Outcome regression (replace by kernel estimator if needed)
  mu_hat <- lm(y ~ V)$fitted.values
  
  if (method == "CC") {
    idx <- which(delta == 1)
    res <- y[idx] - mu_hat[idx]
    V_use <- V[idx, ]
  }
  
  if (method == "IPW") {
    res <- delta * (y - mu_hat) / pi_hat
    V_use <- V
  }
  
  if (method == "AIPW") {
    res <- delta * (y - mu_hat) / pi_hat + mu_hat - mean(mu_hat, na.rm = TRUE)
    V_use <- V
  }
  
  list(res = res, V = V_use)
}

############################################################
## 5. TEST STATISTICS (EDIT THESE)
############################################################

## ----- Classical versions -----

Tn1 <- function(res, V) {
  cor(rank(res), rank(V[,1]), use = "complete.obs")
}

Tn2 <- function(res, V) {
  cor(rank(res), rank(V[,2]), use = "complete.obs")
}

Tn3 <- function(res, V) {
  cor(rank(res), rank(rowMeans(V)), use = "complete.obs")
}

## ----- Spacing-based versions (PLACEHOLDERS) -----

Tn1_spacing <- function(res, V) {
  mean(diff(sort(res)) * diff(sort(V[,1])))
}

Tn2_spacing <- function(res, V) {
  mean(diff(sort(res)) * diff(sort(V[,2])))
}

Tn3_spacing <- function(res, V) {
  mean(diff(sort(res)) * diff(sort(rowMeans(V))))
}

############################################################
## 6. PERMUTATION P-VALUE FUNCTION
############################################################

perm_pvalue <- function(stat_fun, res, V, B = 1000) {
  
  T_obs <- stat_fun(res, V)
  
  T_perm <- replicate(B, {
    perm_res <- sample(res)
    stat_fun(perm_res, V)
  })
  
  mean(abs(T_perm) >= abs(T_obs))
}

############################################################
## 7. MAIN ANALYSIS LOOP
############################################################

missing_levels <- c(0.10, 0.20, 0.30)
methods <- c("CC","IPW","AIPW")

results <- data.frame()

for (miss in missing_levels) {
  
  mar <- generate_MAR(y, V, miss_rate = miss)
  
  for (m in methods) {
    
    obj <- get_residuals(mar$y, V, mar$delta, m)
    
    pvals <- c(
      Tn1  = perm_pvalue(Tn1,  obj$res, obj$V),
      Tn2  = perm_pvalue(Tn2,  obj$res, obj$V),
      Tn3  = perm_pvalue(Tn3,  obj$res, obj$V),
      Tn1S = perm_pvalue(Tn1_spacing, obj$res, obj$V),
      Tn2S = perm_pvalue(Tn2_spacing, obj$res, obj$V),
      Tn3S = perm_pvalue(Tn3_spacing, obj$res, obj$V)
    )
    
    results <- rbind(
      results,
      data.frame(
        Missingness = paste0(100*miss,"%"),
        Method = m,
        t(pvals)
      )
    )
  }
}

############################################################
## 8. FINAL TABLE (LaTeX-READY)
############################################################

results_rounded <- results %>%
  mutate(across(3:8, ~ round(., 3)))

kable(results_rounded,
      caption = "Permutation-based $p$-values under MAR missingness")

############################################################
## END OF SCRIPT
############################################################
