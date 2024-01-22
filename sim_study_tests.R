library(ggplot2)
library(dplyr)

s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(score)
}


mix_forecasts <- function(fcst1, fcst2) {
  pick <- rbinom(length(fcst1), 1, 0.5)

  return(list(
    a = ifelse(pick == 1, fcst1, fcst2),
    b = ifelse(pick == 0, fcst1, fcst2)
  ))
}


geo_test <- function(fcst1, fcst2, y) {
  N <- sum(y)
  diff_logs <- log(fcst1) - log(fcst2)
  I_N_ij <- sum(y * diff_logs) / N - (sum(fcst1) - sum(fcst2)) / N
  test_var <- 1 / (N - 1) * sum((diff_logs - mean(diff_logs))^2)

  test_stat <- I_N_ij / sqrt(test_var) * sqrt(N)
  pval <- 1 - pt(test_stat, N - 1)
  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(test_var), mean_diff = I_N_ij))
}


dm_test <- function(fcst1, fcst2, y) {
  diff_scores <- rowSums(s_pois(fcst1, y) - s_pois(fcst2, y))
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = 7, type = "covariance", plot = F)$acf
  dm_var <- auto_covs[1] + 2 * sum(auto_covs[-1])

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(n_days)
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var),
                    mean_diff = mean_diff_score))
}


sim_tests <- function(B = 100) {
  tests <- list("geo_test" = geo_test, "dm_test" = dm_test)
  collect_results <- data.frame()

  y <- obs
  lm <- models[[which(model_names == "LM")]]
  mix1 <- as.vector(models[[which(model_names == "LG")]])
  mix2 <- as.vector(models[[which(model_names == "FCM")]])
  fcst_names <- c("LM", "MixA", "MixB")

  for (i in 1:B) {
    cat("*")
    mixed_fcsts <- mix_forecasts(mix1, mix2)
    fcst_data <- setNames(
      list(lm, matrix(mixed_fcsts[[1]], nrow = nrow(y)),
           matrix(mixed_fcsts[[2]], nrow = nrow(y))),
      fcst_names
    )

    for (stat_test in names(tests)) {
      for (fcst1 in 1:(length(fcst_names) - 1)) {
        for (fcst2 in (fcst1 + 1):length(fcst_names)) {

          test_vals <- tests[[stat_test]](fcst_data[[fcst_names[fcst1]]],
                                          fcst_data[[fcst_names[fcst2]]],
                                          y)

          collect_results <- rbind(collect_results,
                                   cbind(test_vals, I = i, Type = stat_test,
                                         cmp = paste(fcst_names[fcst1], "vs.", fcst_names[fcst2])))
        }
      }
    }
  }
  return(collect_results)
}


analyze_autocorrelation_fcn <- function(B = 10) {
  acf_data <- data.frame()

  for (i in 1:B) {
    fcst_data <- get_forecasts()
    score_diffs <- s_pois(fcst_data$a, fcst_data$y) - s_pois(fcst_data$b, fcst_data$y)

    acf_obj <- acf(rowSums(score_diffs), type = "covariance", plot = F)

    acf_data <- rbind(
      acf_data, data.frame(lag = acf_obj$lag, value = acf_obj$acf, I = i)
    )
  }
  return(acf_data)
}

B <- 1
acf_data <- analyze_autocorrelation_fcn(B = B)
ggplot(acf_data) +
    facet_wrap()
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(x = lag, y = value, xend = lag, yend = 0)) +
    theme_bw()


plot_results <- function(df) {
  df %>%
    mutate(Type = ifelse(Type == "geo_test", "CSEP T-test",
                         "Diebold-Mariano test"),
           Type = factor(Type, ordered = T,
                         levels = c("CSEP T-test", "Diebold-Mariano test")),
           cmp = factor(cmp, ordered = T,
                        levels = c("LM vs. MixA", "LM vs. MixB", "MixA vs. MixB"))) %>%
    ggplot() +
    facet_grid(Type~cmp) +
    geom_histogram(aes(x = pval, y = ..density..)) +
    xlab("p value") +
    ylab("count") +
    ggtitle("Simulation Study") +
    theme_bw() +
    theme(strip.background = element_blank())
}

cmp_run_sim <- compiler::cmpfun(sim_tests)

collect_results <- cmp_run_sim(B = 100)

pl <- plot_results(collect_results)
pl

ggsave("./../test/sim_study3/ss_tests.pdf", width = 150, height = 120, unit = "mm")
