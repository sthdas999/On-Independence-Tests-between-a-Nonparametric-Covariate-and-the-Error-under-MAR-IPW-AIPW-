
## =====================================================
## Real Data Analysis under MAR
## IPW, AIPW and Spacing-based Independence Test
## Blood Glucose Dataset (Kaggle)
## =====================================================

library(readxl)
library(dplyr)
library(np)

set.seed(123)

## 1. Load data
data <- read_excel("Dataset for People for their Blood Glucose Level with their Superficial body feature readings..xlsx")

## Rename variables (adjust if column names differ)
data <- data %>%
  rename(
    Y  = `Blood Glucose Level`,
    Age = Age,
    BMI = BMI,
    BP  = `Blood Pressure`,
    PA  = `Physical Activity`,
    FH  = `Family History`
  )

## 2. MAR missingness generation (if Y is complete)
pi_true <- plogis(-2 + 0.03 * Age + 0.02 * BP)
delta <- rbinom(nrow(data), 1, pi_true)

data$Y_obs <- ifelse(delta == 1, data$Y, NA)
data$delta <- delta

## 3. Propensity score model
ps_model <- glm(delta ~ Age + BP + PA + FH,
                family = binomial, data = data)

data$ps_hat <- predict(ps_model, type = "response")

## 4. Nonparametric component m(Z)
cc_data <- data %>% filter(!is.na(Y_obs))

np_fit <- npreg(Y_obs ~ factor(PA) + factor(FH),
                regtype = "lc", data = cc_data)

m_hat <- predict(np_fit, newdata = data)

## 5. IPW estimator
X <- model.matrix(~ Age + BMI + BP, data = data)
W <- diag(data$delta / data$ps_hat)

beta_ipw <- solve(t(X) %*% W %*% X) %*%
            t(X) %*% W %*% (data$Y_obs - m_hat)

## 6. AIPW estimator
Y_tilde <- m_hat + data$delta / data$ps_hat * (data$Y_obs - m_hat)

beta_aipw <- solve(t(X) %*% X) %*% t(X) %*% Y_tilde

## 7. Residuals
res_ipw  <- data$Y_obs - X %*% beta_ipw - m_hat
res_aipw <- Y_tilde     - X %*% beta_aipw - m_hat

## 8. Spacing-based Kendall tau
spacing_tau <- function(u, v) {
  u_rank <- rank(u) / (length(u) + 1)
  v_rank <- rank(v) / (length(v) + 1)
  du <- diff(sort(u_rank))
  dv <- diff(sort(v_rank))
  cor(du, dv, method = "kendall")
}

tau_ipw  <- spacing_tau(data$PA[!is.na(res_ipw)],
                        res_ipw[!is.na(res_ipw)])

tau_aipw <- spacing_tau(data$PA,
                        res_aipw)

## 9. Results
cat("IPW estimates:\n")
print(beta_ipw)

cat("\nAIPW estimates:\n")
print(beta_aipw)

cat("\nSpacing-based Kendall tau (IPW residuals):", tau_ipw, "\n")
cat("Spacing-based Kendall tau (AIPW residuals):", tau_aipw, "\n")
