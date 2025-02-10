# plottting : customize theme -------------------------------------
library(ggplot2)

title_size <- 13.2  # base size 11 * 1.2 (default for theme_bw())

my_theme <- list(
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(5, "mm"),
        plot.title = element_text(size = title_size),
        strip.background = element_blank())
)

# Scoring / Evaluation functions ------------------------------------

# handle zero forecast specifically
s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(score)
}

s_pois_da <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(sum(score) / n_days)       # daily averages
}

s_quad <- function(X, Y) {
  return((X - Y)^2)
}

s_quad_da <- function(X, Y) {
  return(sum((X - Y)^2) / n_days)   # daily averages
}

s_theta <- function(x, y, theta) {
  # Elementary scoring function for the mean following Ehm et al. (2016)
  if (sum(y) == 0) {
    s <- theta * sum(x > theta)
  } else {
    s <- sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x) )
  }
  return(s)
}

s_theta_da <- function(x, y, theta) {
  # Elementary scoring function for the mean following Ehm et al. (2016)
  if (sum(y) == 0) {
    s <- theta * sum(x > theta)
  } else {
    s <- sum(pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x))
  }
  return(s / n_days) # daily averages
}

inf_gain <- function(X1, X2, Y) {
  X1_t <- rowSums(X1)
  X2_t <- rowSums(X2)
  return(
    rowSums(Y * (log(X1) - log(X2))) - (X1_t - X2_t)
  )
}

cum_inf_gain <- function(X1, X2, Y) {
  return(cumsum(inf_gain(X1, X2, Y)))
}

inf_gain_per_event <- function(X1, X2, Y) {
  n_t <- rowSums(Y)
  X1_t <- rowSums(X1)
  X2_t <- rowSums(X2)
  return(
    rowSums(Y * (log(X1) - log(X2))) / n_t - (X1_t - X2_t) / n_t
  )
}

cum_igpe <- function(X1, X2, Y) {
  n_t <- cumsum(rowSums(Y))
  X1_t <- cumsum(rowSums(X1))
  X2_t <- cumsum(rowSums(X2))
  return(
    cumsum(rowSums(Y * (log(X1) - log(X2)))) / n_t - (X1_t - X2_t) / n_t
  )
}

igpe <- function(X1, X2, Y) {
  N <- sum(Y)
  return(sum(Y * (log(X1) - log(X2))) / N - sum(X1 - X2) / N)   # "bin-centric" perspective
}

igpe2 <- function(X1, X2, Y) {
  N <- sum(Y)
  y_vec <- as.vector(Y)
  # "event-centric" perspective: each earthquake forms separate summand
  x1_vec <- rep(as.vector(X1), y_vec)
  x2_vec <- rep(as.vector(X2), y_vec)
  return(sum(log(x1_vec) - log(x2_vec)) / N - sum(X1 - X2) / N)
}

# Tests ------------------------------------------------------------------

csep_test <- function(fcst1, fcst2, y, scf = NULL) {
  N <- sum(y)
  diff_logs <- log(fcst1) - log(fcst2)
  I_N_ij <- sum(y * diff_logs) / N - (sum(fcst1) - sum(fcst2)) / N
  test_var <- 1 / (N - 1) * sum(y * diff_logs^2) - 1 / (N^2 - N) * sum(y * diff_logs)^2

  test_stat <- I_N_ij / sqrt(test_var) * sqrt(N)
  pval <- 1 - pt(test_stat, N - 1)
  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(test_var), mean_diff = I_N_ij))
}

dm_test <- function(fcst1, fcst2, y, scf, max_lag = 6) {
  if (!is.null(dim(fcst1)) && dim(fcst1) == 2) {     # forecasts are already aggregated on daily basis
    diff_scores <- scf(fcst2, y) - scf(fcst1, y)
  } else {
    diff_scores <- rowSums(scf(fcst2, y) - scf(fcst1, y))
  }
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = max_lag, type = "covariance", plot = F)$acf
  dm_var <- sum(c(auto_covs[1], 2 * auto_covs[-1]))

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(nrow(fcst1))
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var), mean_diff = mean_diff_score))
}
