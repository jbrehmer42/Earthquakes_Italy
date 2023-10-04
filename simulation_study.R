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
    under_fcst = 0.5 * lambda,
    over_fcst = 1.5 * lambda,
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

# Simulation study run ---------------------------------------------------------

run_simulation_study <- function(n_vals, B = 100) {
  scores <- data.frame()
  reliability <- data.frame()
  murphy <- data.frame()

  # grid for murphy diagram
  n_theta <- 100
  log_grid <- seq(-15, 4, len = n_theta - 1)
  grd <- c(0, 10^log_grid)

  for (n in n_vals) {
    cat("\n:", n, ": ")
    for (b in 1:B) {
      cat("*")
      data <- generate_scenario(n)
      y_save <- data$y          # we keep sorting y, so store a reference version
      fcsts <- select(data, -y)
      mean_y <- mean(y_save)

      for(forecaster in names(fcsts)) {
        x <- fcsts[, forecaster]
        ord <- order(x, y_save, decreasing = c(FALSE, TRUE))
        x <- x[ord]
        y <- y_save[ord]
        x_rc <- monotone(y)

        for (scf_name in names(scf_list)) {
          scf <- scf_list[[scf_name]]
          s <- scf(x, y)
          s_rc <- scf(x_rc, y)
          s_mg <- scf(mean_y, y)

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
        elem_scores_rc <- sapply(grd, function(theta) S_theta(x_rc, y, theta))
        elem_scores_mg <- sapply(grd, function(theta) S_theta(mean_y, y, theta))

        murphy <- rbind(
          murphy,
          data.frame(n = n, fcst = forecaster, I = b, x = rep(grd, 3),
                     Type = rep(c("Score", "MCB", "DSC"), each = length(grd)),
                     y = c(elem_scores, elem_scores - elem_scores_rc,
                           elem_scores_mg - elem_scores_rc))
        )
      }
    }
  }
  return(list(scores = scores, reliability = reliability, murphy = murphy))
}

# Plots ------------------------------------------------------------------------


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
  # expects summary statistics
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

plot_scores <- function(results) {
  pois <- ggplot(filter(results, Scoring == "pois")) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_point(aes(x = fcst, y = value, color = fcst, shape = factor(I))) +
    ggtitle("Poisson Score") +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_log10() +
    labs(shape = "Simmulation run") +
    guides(color = "none", shape = "none") +
    theme_bw() +
    theme(strip.background = element_blank())

  quad <- ggplot(filter(results, Scoring == "quad")) +
    facet_wrap(~Type, scales = "free_y", nrow = 1) +
    geom_point(aes(x = fcst, y = value, color = fcst, shape = factor(I))) +
    ggtitle("Quadratic Score") +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_log10() +
    labs(shape = "Simmulation run") +
    guides(color = "none") +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position = "bottom")

  my_plot <- grid.arrange(pois, quad, nrow = 2, top = "Simulation Study",
                          heights = c(0.45, 0.55))
  return(my_plot)
}

plot_rel <- function(reliability) {
  d <- 10^6
  my_trans <- function(x) sign(x) * log(abs(x) * d  + 1, base = 10)

  my_breaks <- c(0, 10^c(-6, -4, -2, 0))
  my_breaks <- my_trans(my_breaks)
  my_labels <- c("0", paste0("1e", c(-6, -4, -2)), "1")

  my_plot <- ggplot(mapping = aes(x = my_trans(x), y = my_trans(x_rc),
                                  color = fcst, linetype = factor(I))) +
    facet_wrap(~fcst) +
    geom_step(data = reliability, mapping = aes(group = paste(fcst, I)), alpha = 0.5) +
    geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.7,
                  linetype = "dashed") +
    xlab("Forecasted mean") +
    scale_x_continuous(breaks = my_breaks, labels = my_labels) +
    scale_y_continuous(breaks = my_breaks, labels = my_labels) +
    ylab("Conditional mean") +
    labs(color = "Forecaster", linetype = "Simulation run") +
    guides(color = "none") +
    ggtitle("Reliability Diagram") +
    theme_bw() +
    theme(strip.background = element_blank(), legend.position = "bottom")

  return(my_plot)
}

plot_murphy <- function(murphy) {
  d <- 10^10
  my_trans <- function(x) sign(x) * log(abs(x) * d  + 1, base = 10)

  my_breaks <- c(0, 10^c(-8, -6, -4, -2, 0, 2))
  my_breaks <- my_trans(my_breaks)
  my_labels <- c("0", paste0("1e", c(-8, -6, -4, -2)), "1", "1e2")

  my_plot <- ggplot(murphy,
                    aes(x = my_trans(x), y = y, color = fcst)) +
    facet_grid(paste("Run", I)~factor(Type, ordered = T, levels = c("Score", "MCB", "DSC")),
               scales = "free_y") +
    # geom_line(data = murphy, mapping = aes(group = paste(fcst, I)), alpha = 0.1) +
    geom_line() +
    scale_x_continuous(breaks = my_breaks, labels = my_labels) +
    xlab(expression(theta)) +
    ylab("Elementary score") +
    labs(color = "Forecaster", linetype = "Simulation run") +
    ggtitle("Murphy Diagram") +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_blank())

  return(my_plot)
}

# Analysis ---------------------------------------------------------------------

cmp_run_sim <- compiler::cmpfun(run_simulation_study)

l_results <- cmp_run_sim(10^7, B = 5)

df_scores <- l_results$scores %>%
    group_by(n, fcst, Scoring, Type) %>%
    summarise(mean_val = mean(value),
              lower_q = quantile(value, 0.05),
              upper_q = quantile(value, 0.95),
              .groups = "drop")

# my_plot <- plot_results(df_scores)

add <- "exp_n7_B5"
my_plot <- plot_scores(l_results$scores)
ggsave(paste0("./../test/sim_study/sim_", add, "_scores.pdf"),
       width = 240, height = 180, unit = "mm", plot = my_plot)

my_plot <- plot_murphy(l_results$murphy)
ggsave(paste0("./../test/sim_study/sim_", add, "_murphy.pdf"),
       width = 240, height = 180, unit = "mm", plot = my_plot)

my_plot <- plot_rel(l_results$reliability)
ggsave(paste0("./../test/sim_study/sim_", add, "_reliability.pdf"),
       width = 240, height = 180, unit = "mm", plot = my_plot)

data <- generate_scenario(10^7)
table(data$y)
print(mean(data$mean_fcst > data$less_disc))


write.csv(l_results$scores, paste0("./../test/sim_study/results_", add, "_scores.csv"))
write.csv(l_results$murphy, paste0("./../test/sim_study/results_", add, "_murphy.csv"))
write.csv(l_results$reliability, paste0("./../test/sim_study/results_", add, "_rel.csv"))
# results <- read.csv("./../test/sim_study/results_4.csv")
