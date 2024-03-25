library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(lubridate)

s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(score)
}

fpath <- "./figures9"
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

################################################################################

get_mix_forecasts <- function(fcst1, fcst2) {
  pick <- rbinom(n_days, 1, 0.5)
  a <- matrix(NA, nrow = nrow(fcst1), ncol = ncol(fcst1))
  a[pick == 1, ] <- fcst1[pick == 1, ]
  a[pick == 0, ] <- fcst2[pick == 0, ]
  b <- matrix(NA, nrow = nrow(fcst1), ncol = ncol(fcst1))
  b[pick == 1, ] <- fcst2[pick == 1, ]
  b[pick == 0, ] <- fcst1[pick == 0, ]
  return(list(a = a, b = b))
}

csep_test <- function(fcst1, fcst2, y) {
  N <- sum(y)
  diff_logs <- log(fcst1) - log(fcst2)
  I_N_ij <- sum(y * diff_logs) / N - (sum(fcst1) - sum(fcst2)) / N
  test_var <- 1 / (N - 1) * sum(y * diff_logs^2) - 1 / (N^2 - N) * sum(y * diff_logs)^2

  test_stat <- I_N_ij / sqrt(test_var) * sqrt(N)
  pval <- 1 - pt(test_stat, N - 1)
  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(test_var), mean_diff = I_N_ij))
}

dm_test <- function(fcst1, fcst2, y, scf) {
  diff_scores <- rowSums(scf(fcst2, y) - scf(fcst1, y))
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = 6, type = "covariance", plot = F)$acf
  dm_var <- auto_covs[1] + 2 * sum(auto_covs[-1])

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(n_days)
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var),
                    mean_diff = mean_diff_score))
}

sim_tests <- function(B = 400, weekday = "Mo") {
  collect_results <- data.frame()

  filter_days <- (lubridate::wday(times, label = T) == weekday)

  y <- obs
  lm <- models[[which(model_names == "LM")]]
  mix1 <- models[[which(model_names == "LG")]]
  mix2 <- models[[which(model_names == "FCM")]]
  fcst_names <- c("LM", "MixA", "MixB")

  for (i in 1:B) {
    cat("*")
    mixed_fcsts <- get_mix_forecasts(mix1, mix2)

    fcst_data <- setNames(list(lm, mixed_fcsts[[1]], mixed_fcsts[[2]]), fcst_names)

    for (fcst1 in 1:(length(fcst_names) - 1)) {
      for (fcst2 in (fcst1 + 1):length(fcst_names)) {
        cmp_name <- paste(fcst_names[fcst1], "vs.", fcst_names[fcst2])
        test_vals1 <- csep_test(fcst_data[[fcst_names[fcst1]]], fcst_data[[fcst_names[fcst2]]], y)
        test_vals2 <- dm_test(fcst_data[[fcst_names[fcst1]]], fcst_data[[fcst_names[fcst2]]], y, s_pois)
        test_vals3 <- csep_test(fcst_data[[fcst_names[fcst1]]][filter_days,],
                                fcst_data[[fcst_names[fcst2]]][filter_days,], y[filter_days,])
        collect_results <- rbind(
          collect_results,
          cbind(test_vals1, I = i, cmp = cmp_name, t = "CSEP", m = "all"),
          cbind(test_vals2, I = i, cmp = cmp_name, t = "DM", m = "all"),
          cbind(test_vals3, I = i, cmp = cmp_name, t = "CSEP", m = "Mo")
        )
      }
    }
  }
  return(collect_results)
}

# Plots ------------------------------------------------------------------------

analyze_autocorrelation_fcn <- function(B = 10, s_func = s_pois) {
  acf_data <- data.frame()

  for (i in 1:B) {
    fcst_data <- mix_forecasts(as.vector(models[[which(model_names == "LG")]]),
                               as.vector(models[[which(model_names == "FCM")]]))
    fcst_a <- matrix(fcst_data$a, nrow = nrow(obs))
    fcst_b <- matrix(fcst_data$b, nrow = nrow(obs))
    score_diffs <- s_func(fcst_a, obs) - s_func(fcst_b, obs)

    acf_obj <- acf(rowSums(score_diffs), type = "covariance", plot = F)

    acf_data <- rbind(
      acf_data, data.frame(lag = acf_obj$lag, value = acf_obj$acf, I = i)
    )
  }

  ggplot(acf_data) +
    facet_wrap(~I) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(aes(x = lag, y = value, xend = lag, yend = 0)) +
    theme_bw()
}

plot_results <- function(df) {
  bins <- seq(0, 1, length.out = 20 + 1)

  new_cmp <- c("MixA vs. MixB" = "MixA vs. MixB", "LM vs. MixA" = "MixA vs. LM",
               "LM vs. MixB" = "MixB vs. LM")
  ord <- c("MixA vs. MixB", "MixA vs. LM", "MixB vs. LM")

  df %>%
    mutate(cmp = factor(new_cmp[cmp], ordered = T, levels = ord)) %>%
    ggplot() +
    #  facet_grid(Type~cmp) +
    facet_wrap(~cmp, scale="free_y") +
    geom_histogram(aes(x = pval, y = ..density..), breaks = bins) +
    geom_vline(xintercept = c(0.05, 0.95), color = "#f8766d", linetype = "dashed",
               size = 0.3) +
    scale_x_continuous(breaks = 0:4 / 4, labels = c("0", "0.25", "0.5", "0.75", "1")) +
    #scale_y_continuous(breaks = NULL, labels = NULL, minor_breaks = (0:8) * 2.5) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(NULL) +
    my_theme +
    theme(strip.background = element_blank())
}

cmp_run_sim <- compiler::cmpfun(sim_tests)

set.seed(999)
results <- cmp_run_sim(B = 200)
write.csv(results, "./figures9/simstudy_tests_400.csv")

my_plot <- plot_results(filter(results, t == "DM", m == "all"))
file_path <- file.path(fpath, "Fig5_sim-DieboldMariano.pdf")
ggsave(file_path, width = 140, height = 50, unit = "mm")

my_plot <- plot_results(filter(results, t == "CSEP", m == "all"))
file_path <- file.path(fpath, "Fig_CSEP-t-Test.pdf")
ggsave(file_path, width = 140, height = 50, unit = "mm")

my_plot <- plot_results(filter(results, t == "CSEP", m == "Mo"))
file_path <- file.path(fpath, "Fig_CSEP-t-Test-Mo.pdf")
ggsave(file_path, width = 140, height = 50, unit = "mm")

results %>%
  filter(cmp == "MixA vs. MixB", t == "CSEP") %>%
  summarise(A = sum(pval < 0.05), B = sum(pval > 0.95))

rm(cmp_run_sim, cmp_cons_fcsts, results, my_plot)
