library(monotone)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

set.seed(111)

MEAN_OBS <- 3.656171e-05

s_quad <- function(x, y) return(mean((x - y)^2))
s_pois <- function(X, y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (y != 0)
  score <- -y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(mean(score))
}

scf_list <- list(quad = s_quad, pois = s_pois)

generate_scenario <- function(n) {
  lambda_0 <- MEAN_OBS

  z <- rnorm(n)
  mu <- 2 * lambda_0 * pnorm(z)
  lambda <- rexp(n, rate = 1 / mu)
  y <- rpois(n, lambda)

  df <- data.frame(
    y = y,
    mean_fcst = lambda,
    under_fcst = 0.5 * lambda,
    over_fcst = 1.5 * lambda,
    less_disc = mu
  )
}

generate_scenario <- function(n) {
  lambda_0 <- MEAN_OBS

  z <- rnorm(n)
  mu <- 2 * lambda_0 * pnorm(z)
  lambda <- rgamma(n, shape = mu, rate = 1)
  y <- rpois(n, lambda)

  df <- data.frame(
    y = y,
    mean_fcst = lambda,
    under_fcst = 0.99 * lambda,
    over_fcst = 1.01 * lambda,
    less_disc = mu
  )
}

filter_jumps <- function(v) {
  n <- length(v)
  return(c(T, v[c(-1, -n)] - v[c(-n+1, -n)] > 0, T))
}

S_theta <- function(x, y, theta) {
  # Elementary scoring function for the mean following Ehm et al. (2016)
  if ( sum(y) == 0 ) {
    s <- theta * sum(x > theta)
  } else {
    s <- sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x) )
  }
  return(s)
}

run_simulation_study <- function(n_vals, B = 100) {
  scores <- data.frame()
  reliability <- data.frame()
  murphy <- data.frame()

  # grid for murphy diagram
  n_theta <- 100
  log_grid <- seq(-24, 4, len = n_theta)
  grd <- exp(log_grid)

  for (n in n_vals) {
    cat("\n:", n, ": ")
    for (b in 1:B) {
      cat("*")
      data <- generate_scenario(n)
      y <- data$y
      fcsts <- select(data, -y)

      for(forecaster in names(fcsts)) {
        x <- fcsts[, forecaster]
        ord <- order(x, y, decreasing = c(FALSE, TRUE))
        x <- x[ord]
        y <- y[ord]
        x_rc <- monotone(y)

        for (scf_name in names(scf_list)) {
          scf <- scf_list[[scf_name]]
          s <- scf(x, y)
          s_rc <- scf(x_rc, y)
          s_mg <- scf(mean(y), y)

          scores <- rbind(
            scores,
            data.frame(n = n, fcst = forecaster, Scoring = scf_name,
                       Type = c("Score", "MCB", "DSC"), I = b,
                       value = c(s, s - s_rc, s_mg - s_rc))
          )
        }
        filt_jumps <- filter_jumps(x_rc)

        reliability <- rbind(
          reliability,
          data.frame(n = n, fcst = forecaster, I = b, x = x[filt_jumps],
                     x_rc = x_rc[filt_jumps])
        )

        elem_scores <- sapply(grd, function(theta) S_theta(x, y, theta))

        murphy <- rbind(
          murphy,
          data.frame(n = n, fcst = forecaster, I = b, x = log_grid, y = elem_scores)
        )
      }
    }
  }
  return(list(scores = scores, reliability = reliability, murphy = murphy))
}


plot_forecaster <- function(data) {
  data %>%
    arrange(mean_fcst) %>%
    mutate(y = factor(y), X = 1:nrow(.)) %>%
    pivot_longer(cols = c(-X, -y), names_to="Forecaster") %>%
    ggplot() +
    facet_wrap(~y, scales = "free") +
    geom_point(aes(x = X, y = value, color = Forecaster)) +
    theme_bw() +
    theme(strip.background = element_blank())
  ggsave("./../test/sim_study/forecasts.pdf", width = 240, height = 150, unit = "mm")
}

plot_results <- function(results) {
  pois <- ggplot(filter(results, Scoring == "pois")) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_ribbon(aes(x = n, ymin = lower_q, ymax = upper_q, fill = fcst),
                alpha = 0.25, show.legend = F) +
    geom_line(aes(x = n, y = mean_val, color = fcst, group = fcst),
              show.legend = F, size = 0.7) +
    ggtitle("Poisson Score") +
    xlab("n") +
    ylab("score") +
    scale_x_log10() +
    theme_bw() +
    theme(strip.background = element_blank())

  quad <- ggplot(filter(results, Scoring == "quad")) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_ribbon(aes(x = n, ymin = lower_q, ymax = upper_q, fill = fcst),
                alpha = 0.25) +
    geom_line(aes(x = n, y = mean_val, color = fcst, group = fcst), size = 0.7) +
    ggtitle("Quadratic Score") +
    xlab("n") +
    ylab("score") +
    labs(color = "Forecaster", fill = "Forecaster") +
    scale_x_log10() +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position = "bottom")

  my_plot <- grid.arrange(pois, quad, nrow = 2, heights = c(0.45, 0.55),
                          top = "Simulation Study")
  return(my_plot)
}

plot_one_n <- function(results, n_pick) {
  pois <- ggplot(filter(results, Scoring == "pois", n == n_pick)) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_pointrange(aes(x = fcst, y = mean_val, ymin = lower_q, ymax = upper_q,
                        color = fcst), size = 0.75, show.legend = F) +
    ggtitle("Poisson Score") +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_log10() +
    theme_bw() +
    theme(strip.background = element_blank())

  quad <- ggplot(filter(results, Scoring == "quad", n == n_pick)) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_pointrange(aes(x = fcst, y = mean_val, ymin = lower_q, ymax = upper_q,
                        color = fcst), size = 0.75, show.legend = F) +
    ggtitle("Quadratic Score") +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_log10() +
    theme_bw() +
    theme(strip.background = element_blank())

  my_plot <- grid.arrange(pois, quad, nrow = 2, top = "Simulation Study")
  return(my_plot)
}

cmp_run_sim <- compiler::cmpfun(run_simulation_study)

l_results <- cmp_run_sim(10^7, B = 10)


results_summary <- scores %>%
    group_by(n, fcst, Scoring, Type) %>%
    summarise(mean_val = mean(value),
              lower_q = quantile(value, 0.05),
              upper_q = quantile(value, 0.95),
              .groups = "drop")
my_plot <- plot_results(results)

my_plot <- plot_one_n(results, 10^7)

my_plot

write.csv(results, "./../test/sim_study/results_4.csv")
results <- read.csv("./../test/sim_study/results_4.csv")
ggsave("./../test/sim_study/simulation_study_4_log.pdf", width = 240, height = 180,
       unit = "mm", plot = my_plot)
