############################################################
# Real Data Driven Study: Blood Glucose Level Analysis
# Semiparametric IPW / AIPW under MAR
# Missingness only in Y
# Continuous nonparametric covariates
############################################################

# -----------------------------
# 1. Load required libraries
# -----------------------------
library(tidyverse)
library(np)
library(survey)

# -----------------------------
# 2. Read and preprocess data
# -----------------------------
# Update file path if necessary
data <- read.csv("blood_glucose.csv")

# Rename variables (typical Kaggle diabetes naming)
data <- data %>%
  rename(
    Y        = Glucose,
    Age      = Age,
    BMI      = BMI,
    BP       = BloodPressure,
    Insulin  = Insulin
  )

# Keep only continuous covariates
data <- data %>%
  select(Y, Age, BMI, BP, Insulin)

# Remove biologically invalid zero values
data <- data %>%
  filter(BMI > 0, BP > 0, Insulin > 0)

cat("Sample size after cleaning:", nrow(data), "\n")

# -----------------------------
# 3. Induce MAR missingness in Y
# -----------------------------
set.seed(123)

# MAR mechanism depending only on covariates
linpred <- -1.5 + 0.02 * data$Age + 0.015 * data$BP
pi_obs  <- plogis(linpred)

# Missingness indicator
data$delta <- rbinom(nrow(data), size = 1, prob = pi_obs)

# Observed response
data$Y_obs <- ifelse(data$delta == 1, data$Y, NA)

cat("Observed proportion of Y:", mean(data$delta), "\n")

# -----------------------------
# 4. Propensity score estimation
# -----------------------------
ps_model <- glm(delta ~ Age + BMI + BP + Insulin,
                family = binomial,
                data = data)

data$pi_hat <- predict(ps_model, type = "response")

# -----------------------------
# 5. IPW nonparametric regression
# -----------------------------
np_ipw <- npreg(
  Y_obs ~ Age + BMI + BP + Insulin,
  data = data,
  weights = data$delta / data$pi_hat
)

data$m_hat_ipw <- fitted(np_ipw)

# -----------------------------
# 6. AIPW estimation
# -----------------------------
# Outcome regression using complete cases
np_or <- npreg(
  Y_obs ~ Age + BMI + BP + Insulin,
  data = data[data$delta == 1, ]
)

data$m_hat_or <- predict(np_or, newdata = data)

# AIPW pseudo-outcome
data$Y_aipw <- with(data,
                    delta * Y_obs / pi_hat +
                      (1 - delta / pi_hat) * m_hat_or
)

# Final AIPW nonparametric fit
np_aipw <- npreg(
  Y_aipw ~ Age + BMI + BP + Insulin,
  data = data
)

data$m_hat_aipw <- fitted(np_aipw)

# -----------------------------
# 7. Spacing-based residuals
# -----------------------------
data$resid_aipw <- data$Y_aipw - data$m_hat_aipw

# Order by a continuous covariate (Age)
ord <- order(data$Age)

V_ord <- data$Age[ord]
R_ord <- data$resid_aipw[ord]

# Third-order spacing residuals
D3 <- R_ord[4:length(R_ord)] -
  3 * R_ord[3:(length(R_ord)-1)] +
  3 * R_ord[2:(length(R_ord)-2)] -
  R_ord[1:(length(R_ord)-3)]

# -----------------------------
# 8. Diagnostic plot
# -----------------------------
plot(V_ord[-(1:3)], D3,
     xlab = "Age (ordered)",
     ylab = "Third-order spacing residual",
     main = "Spacing-based residual diagnostic")
abline(h = 0, col = "red")

cat("Analysis completed successfully.\n")
############################################################
