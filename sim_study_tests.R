library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(lubridate)

fpath <- "./figures"
tpath <- "./save_results"

source("data_prep.R")   # make sure to load data of LM, LG and FMC model
source("utils.R")

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

# Plots functions ------------------------------------------------------------------------

plot_results <- function(df) {
  bins <- seq(0, 1, length.out = 20 + 1)

  new_cmp <- c("MixA vs. MixB" = "MixA vs. MixB", "LM vs. MixA" = "MixA vs. LM",
               "LM vs. MixB" = "MixB vs. LM")
  ord <- c("MixA vs. MixB", "MixA vs. LM", "MixB vs. LM")

  df %>%
    mutate(cmp = factor(new_cmp[cmp], ordered = T, levels = ord)) %>%
    ggplot() +
    facet_wrap(~cmp) +
    geom_histogram(aes(x = pval, y = ..density..), breaks = bins) +
    geom_vline(xintercept = c(0.05, 0.95), color = "#f8766d", linetype = "dashed",
               size = 0.3) +
    scale_x_continuous(breaks = 0:4 / 4, labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(breaks = NULL, labels = NULL, minor_breaks = (0:8) * 2.5) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(NULL) +
    my_theme +
    theme(strip.background = element_blank())
}

cmp_run_sim <- compiler::cmpfun(sim_tests)

set.seed(1)
results <- cmp_run_sim(B = 400)
write.csv(results, file.path(tpath, "simstudy-tests_400.csv"))

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
  filter(cmp == "MixA vs. MixB") %>%
  group_by(t, m) %>%
  summarise(A = sum(pval < 0.05, na.rm = T), B = sum(pval > 0.95, na.rm = T))


# investigate CSEP : all vs. Mondays -----------------------------------------------------------------------------------
new_names <- c("mean_diff" = "IGPE", "sd" = "standard dev.", "zval" = "z-statistic", "pval" = "p-value")
filt <- "MixA vs. MixB"
x1 <- c(-7, qt(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 259 - 1), 7)
x2 <- c(-7, qt(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 1813 - 1), 7)
n <- length(x1)
pval_band <- data.frame(metric = factor("z-statistic", ordered = T, levels = unname(new_names)),
                        ymin = -7.5, ymax = -7.25, min1 = x1[-n], max1 = x1[-1], min2 = x2[-n], max2 = x2[-1])
results %>%
  filter(t == "CSEP", cmp == filt) %>%
  pivot_longer(cols = c(mean_diff, zval, pval, sd), names_to = "metric") %>%
  pivot_wider(id_cols = c(I, metric), names_from = m, values_from = value) %>%
  mutate(metric = factor(new_names[metric], ordered = T, levels = unname(new_names))) %>%
  ggplot() +
  facet_wrap(~metric, scales = "free") +
  geom_point(aes(x = Mo, y = all), alpha = 0.5) +
  geom_rect(data = pval_band, aes(xmin = min1, xmax = max1, ymin = ymin, ymax = ymax), color = "black", fill = NA) +
  geom_rect(data = pval_band, aes(xmin = ymin, xmax = ymax, ymin = min2, ymax = max2), color = "black", fill = NA) +
  ggtitle("CSEP Tests") +
  labs(subtitle = filt, color = "p-value", fill = "p-value") +
  xlab("Filterd for Monday") +
  ylab("All values") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("./../test/filter_weekday/csep_mixAmixB.pdf", width = 240, height = 200, unit = "mm")

results %>%
  filter(t == "CSEP") %>%
  pivot_longer(cols = c(mean_diff, zval, pval, sd), names_to = "metric") %>%
  pivot_wider(id_cols = c(I, metric, cmp), names_from = m, values_from = value) %>%
  mutate(metric = factor(new_names[metric], ordered = T, levels = unname(new_names))) %>%
  ggplot() +
  facet_wrap(~metric, scales = "free") +
  geom_point(aes(x = Mo, y = all, color = cmp), alpha = 0.5) +
  ggtitle("CSEP Tests") +
  xlab("Filterd for Monday") +
  ylab("All values") +
  labs(color = "Comparison") +
  theme_bw() +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave("./../test/filter_weekday/csep_all.pdf", width = 240, height = 200, unit = "mm")

results %>%
  filter(t == "CSEP", cmp == filt) %>%
  pivot_longer(cols = c(mean_diff, zval, pval, sd), names_to = "metric") %>%
  mutate(metric = factor(new_names[metric], ordered = T, levels = unname(new_names))) %>%
  ggplot() +
  facet_wrap(~metric, scales = "free") +
  geom_histogram(aes(x = value, fill = m), position = "dodge", bins = 20) +
  labs(fill = "Version", subtitle = filt) +
  xlab(NULL) +
  ggtitle("Marginal Distributions") +
  theme_bw() +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave("./../test/filter_weekday/csep_marginals_lmmixB.pdf", width = 240, height = 200, unit = "mm")

rm(cmp_run_sim, cmp_cons_fcsts, results, my_plot)
