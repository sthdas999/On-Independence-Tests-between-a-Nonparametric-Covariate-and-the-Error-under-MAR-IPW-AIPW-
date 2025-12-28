
# sim_study.R
# Simulation framework for finite-sample size and power study
# Requires: tidyverse, foreach, doParallel, broom, xtable

library(tidyverse)
library(foreach)
library(doParallel)
library(broom)
library(xtable)

# ---------- User choices ----------
ns <- c(100, 200, 400)
R <- 500           # Monte Carlo replications
B <- 500           # bootstrap samples
alpha_vals <- c(0, 1)    # spacing exponent options
methods <- c("NW","ILLS","IPW","AIPW","IPW-S","AIPW-S")
tests <- c("Tn1","Tn2","Tn3")
alternatives <- c("null","symmetric","asymmetric")
alt_strengths <- c(0, 0.5, 1)
seed0 <- 20251212

# ---------- Parallel setup ----------
ncores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncores)
registerDoParallel(cl)

set.seed(seed0)

# ---------- Data generating mechanism ----------
generate_data <- function(n, alt = "null", alt_strength = 0) {
  X <- runif(n, -1, 1)
  m0 <- function(x) 0.5 * x

  if (alt == "symmetric") {
    eps <- rnorm(n, 0, sqrt(1 + alt_strength * X^2))
  } else if (alt == "asymmetric") {
    eps_raw <- rchisq(n, df = 1 + alt_strength) - (1 + alt_strength)
    eps <- scale(eps_raw)[, 1]
  } else {
    eps <- rnorm(n, 0, 1)
  }

  Y <- m0(X) + eps
  data.frame(X = X, Y = Y)
}

# ---------- Estimation / imputation methods ----------
estimate_method <- function(dat, method_label) {
  n <- nrow(dat)

  if (method_label == "NW") {
    fit <- loess(Y ~ X, data = dat, span = 0.5)
    pred <- predict(fit, dat$X)

  } else if (method_label == "ILLS") {
    fit <- loess(Y ~ X, data = dat, span = 0.3)
    pred <- predict(fit, dat$X)

  } else if (method_label %in% c("IPW","IPW-S")) {
    prop <- rep(0.8, n)
    pred <- rep(mean(dat$Y / prop), n)

  } else if (method_label %in% c("AIPW","AIPW-S")) {
    fit <- lm(Y ~ X, data = dat)
    reg_pred <- predict(fit, dat)
    prop <- rep(0.8, n)
    pred <- reg_pred + (dat$Y - reg_pred) / prop

  } else {
    pred <- rep(mean(dat$Y), n)
  }

  list(pred = pred)
}

# ---------- Test statistics ----------
compute_Tn1 <- function(dat, est) {
  res <- dat$Y - est$pred
  var(res)
}

compute_Tn2 <- function(dat, est, alpha = 0) {
  res <- dat$Y - est$pred
  s <- sort(abs(res))
  mean(s^alpha)
}

compute_Tn3 <- function(dat, est) {
  res <- dat$Y - est$pred
  mean(res^2)
}

# ---------- Bootstrap p-value ----------
bootstrap_pvalue <- function(dat, est_fun, T_fun, B = 500,
                              alpha = NULL, method_label = NULL) {
  est0 <- est_fun(dat, method_label)
  T0 <- if (is.null(alpha)) T_fun(dat, est0) else T_fun(dat, est0, alpha)

  n <- nrow(dat)
  Tb <- numeric(B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    datb <- dat[idx, , drop = FALSE]
    estb <- est_fun(datb, method_label)
    Tb[b] <- if (is.null(alpha)) T_fun(datb, estb) else T_fun(datb, estb, alpha)
  }

  mean(abs(Tb - mean(Tb)) >= abs(T0 - mean(Tb)))
}

# ---------- Single replication ----------
simulate_one <- function(n, alt_type, alt_strength,
                         method_label, alpha, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  dat <- generate_data(n, alt = alt_type, alt_strength = alt_strength)

  tibble(
    n = n,
    alt = alt_type,
    alt_strength = alt_strength,
    method = method_label,
    alpha = alpha,
    p_Tn1 = bootstrap_pvalue(dat, estimate_method, compute_Tn1,
                             B = B, method_label = method_label),
    p_Tn2 = bootstrap_pvalue(dat, estimate_method, compute_Tn2,
                             B = B, alpha = alpha, method_label = method_label),
    p_Tn3 = bootstrap_pvalue(dat, estimate_method, compute_Tn3,
                             B = B, method_label = method_label)
  )
}

# ---------- Main simulation loop ----------
results <- foreach(n = ns, .combine = bind_rows,
                   .packages = c("tidyverse")) %:%
  foreach(alt = alternatives, .combine = bind_rows) %:%
  foreach(delta = alt_strengths, .combine = bind_rows) %:%
  foreach(method = methods, .combine = bind_rows) %:%
  foreach(alpha = alpha_vals, .combine = bind_rows) %dopar% {
    map_dfr(seq_len(R), function(r) {
      simulate_one(n, alt, delta, method, alpha,
                   seed = seed0 + r)
    })
  }

# ---------- Shutdown ----------
stopCluster(cl)

# Optional save
# saveRDS(results, file = "simulation_results.rds")
