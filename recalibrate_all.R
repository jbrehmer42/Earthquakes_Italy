library(dplyr)
library(data.table)
library(tidyr)
library(monotone)
library(scales)
library(ggplot2)
library(gridExtra)
library(grid)

# write zero score as it will not change sums
my_s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(sum(score) / n_days)
}

fpath <- "./figures_rec"

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")
title_size <- 13.2  # base size 11 * 1.2 (default for theme_bw())
my_theme <- list(
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = title_size),
        strip.background = element_blank())
)

# recalibration routines =======================================================

# adapted from
# https://github.com/dwolffram/replication-ARSIA2023/blob/main/R/reliability_functions.R

# we can compress step function by only storing first and last value and jumps!
filter_jumps <- function(v) {
  n <- length(v)
  return(c(T, v[c(-1, -n)] - v[c(-n+1, -n)] > 0, T))
}

fit_isotonic_all <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)

  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)
  return(data.table(x = x, y = y, x_rc = x_rc))
}

fit_isotonic_by_cell <- function(x, y) {
  # fit isotonic regression for each cell separately
  return(do.call(
    rbind, lapply(1:ncol(x), function(c) fit_isotonic_all(x[, c], y[, c]))
  ))
}

resample_residuals <- function(x, y, B, get_iso_fit) {
  collect_vals <- list()

  res <- y - x
  mean_res <- mean(res)
  for (i in 1:B) {
    y <- x + matrix(sample(as.vector(res), prod(dim(y)), replace = TRUE),
                    nrow = nrow(y), ncol = ncol(y))

    collect_vals[[i]] <- get_iso_fit(x, y) %>%
      filter(filter_jumps(x_rc)) %>%
      transmute(x = x, y = pmax(0, x_rc - mean_res), I = paste0("R", i))
  }
  return(collect_vals)
}

resample_daily_residuals <- function(x, y, B, get_iso_fit) {
  # resample rows, i.e. each acknowledge spatial distribution of errors
  collect_vals <- list()

  mean_res <- mean(res)
  for (i in 1:B) {
    y <- x + res[sample(1:nrow(res), nrow(res), replace = TRUE), ]
    collect_vals[[i]] <- get_iso_fit(x, y) %>%
      filter(filter_jumps(x_rc)) %>%
      transmute(x = x, y = pmax(0, x_rc - mean_res), I = paste0("R", i))
  }
  return(collect_vals)
}

resample_trams_residuals <- function(x, y, B, get_iso_fit) {
  # resample transformed resiudals, transformation has to be defined in global scope
  collect_vals <- list()

  res <- trans$transform(y) - trans$transform(x)
  for (i in 1:B) {
    y <- trans$inverse(trans$transform(x) + matrix(sample(as.vector(res), prod(dim(y)), replace = TRUE),
                                                   nrow = nrow(y), ncol = ncol(y)))

    collect_vals[[i]] <- get_iso_fit(x, y) %>%
      filter(filter_jumps(x_rc)) %>%
      transmute(x = x, y = pmax(0, x_rc - mean_res), I = paste0("R", i))
  }
  return(collect_vals)
}

resample_adapt_cond_distr <- function(x, y, B, get_iso_fit) {
  x <- as.vector(x)
  y <- as.vector(y)
  collect_vals <- list()

  t <- table(y) / length(y)
  vals <- as.numeric(names(t))
  f_y <- as.numeric(t)
  eps <- sum(f_y[-1] * vals[-1]) / x - 1
  n_split <- 5      # as not everything fits into memory, split
  split <- ceiling(seq(1, length(y), length.out = n_split))
  cmp <- matrix(rep(cumsum(f_y[-length(f_y)]), split[2] - split[1] + 1),
                ncol = length(f_y) - 1, byrow = T)
  for (i in 1:B) {
    for (j in 2:n_split) {
      i_vec <- split[j - 1]:split[j]
      u <- runif(length(i_vec), 0, 1 + eps[i_vec])
      y[i_vec] <- rowSums(u - eps[i_vec] >= cmp[1:pmin(length(i_vec), nrow(cmp)), ])
    }
    collect_vals[[i]] <- get_iso_fit(x, y) %>%
      filter(filter_jumps(x_rc)) %>%
      transmute(x = x, y = pmax(0, x_rc), I = paste0("R", i))
  }
  return(collect_vals)
}

reldiag <- function(x, y, n_resamples = 99, region_level = 0.9, score = my_s_pois,
                    get_iso_fit = fit_isotonic_all, res_for_cons = resample_residuals) {
  my_iso_fit <- get_iso_fit(x, y)
  s <- score(my_iso_fit$x, my_iso_fit$y)
  s_rc <- score(my_iso_fit$x_rc, my_iso_fit$y)
  s_mg <- score(mean(my_iso_fit$y), my_iso_fit$y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  low <- floor(n_resamples * (1 - region_level) / 2)
  up <- n_resamples - low
  low <- pmax(low, 1)

  collect_vals <- list(my_iso_fit %>%
    filter(filter_jumps(x_rc)) %>% transmute(x = x, y = x_rc, I = "Fit"))
  # bootstrap fits for consistency bands
  bootrstrap_fits <- res_for_cons(x, y, n_resamples, get_iso_fit)
  collect_vals <- append(collect_vals, bootrstrap_fits)

  # build data table with the collected values on a joint grid
  results <- do.call(rbind, collect_vals) %>%
    pivot_wider(id_cols = x, values_from = y, names_from = I, values_fill = NA) %>%
    arrange(x) %>%                  # sort x-value
    setnafill(type = "locf") %>%    # use last available observation to fill NA (step function!)
    apply(1, function(row) {
      sorted <- sort(row[c(-1, -2)])
      c(row[1], row[2], sorted[low], sorted[up])
    }) %>% t() %>%  # apply writes results in columns
    as.data.table()

  colnames(results) <- c("x", "x_rc", "lower", "upper")
  stats <- data.frame(Score = s, MCB = mcb, DSC = dsc, UNC = unc)
  return(list(results = results, stats = stats))
}

# recalibrate values ===========================================================

# Define different settings: pick the one you want to execute

# Setting 1: Analyze all values ------------------------------------------------
my_obs <- obs
my_models <- models
n_resamples <- 20

# only use every 7-th value, as forecasts span 7-day periods------------
mod <- 7      # use 1 to 7
seventh <- 0:floor(nrow(obs) / 7)
# for a random control: mod <- sample(1:7, length(seventh), replace = T)
pick <- (7 * seventh + mod)[7 * seventh + mod <= nrow(obs)]
weekday <- lubridate::wday(times[mod], label = T, abbr = F)
my_obs <- obs[pick, ]
my_models <- lapply(models, function(x) x[pick, ])
n_resamples <- 20

# in any case define breaks and labels
breaks_x <- list(c(0, 10^c(-5, 0)), c(0, 10^c(-8, -6, -5, 0)),
                 c(0, 10^c(-6, -5, 0)), c(0, 10^c(-6, -5, 0)))
breaks_y <- list(c(0, 10^c(-6, -5, -4, 0)), c(0, 10^c(-8, -7, -6, -5, -4, 0)),
                 c(0, 10^c(-6, -5, 0)), c(0, 10^c(-6, -5, -4, 0)))
minor_breaks <- c(0, 10^(-10:0))
my_labeller <- function(l) {
  labels <- character(length(l))
  labels[l > 0 & l < 1] <- paste0("1e", log(l[l > 0 & l < 1], base = 10))
  labels[l == 1] <- "1"
  labels[l == 0] <- "0"
  return(labels)
}

# daily values -----------------------------------------------------------------
my_obs <- matrix(rowSums(obs), ncol = 1)
my_models <- lapply(models, function(x) matrix(rowSums(x), ncol = 1))
n_resamples <- 5000

breaks_x <- list(c(0, 0.2, 0.3, 4), c(0, 0.2, 0.3, 4),
                 c(0, 0.2, 0.3, 4), c(0, 0.2, 0.3, 4))
breaks_y <- list(c(0, 0.2, 0.3, 4), c(0 ,0.3, 4),
                 c(0, 0.3, 4), c(0, 0.3, 4))
minor_breaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
my_labeller <- function(l) paste(l)

# ------------------------------------------------------------------------------


pick_iso <- fit_isotonic_all

pick_res <- resample_adapt_cond_distr
pick_res <- resample_residuals

reldiag_cmp <- compiler::cmpfun(function(x, y) {
  reldiag(x, y, n_resamples, get_iso_fit = pick_iso, res_for_cons = pick_res)
})

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  res <- reldiag_cmp(my_models[[i]], my_obs)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.2e", res$stats[1, ]),
                        collapse = "\n"))
  )
}

# or just load already recalibrated values
recal_models <- read.csv("./../tmp_results/recal_models_100_bag-Tilmann.csv")
collect_stats <- read.csv("./../tmp_results/collect_stats.csv")

# plot calibration curve under certain data transformations ====================
# use empirical CDF transform
col_ecdfs <- list()
for (i in 1:length(my_models)) {
  t <- table(my_models[[i]])
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / prod(dim(my_models[[i]])),
                               M = model_names[i])
}

# use individual empirical CDFs for plotting -----------------------------------

collect_recal_plots <- list()

for (i in 1:length(my_models)) {
  df_ecdf <- col_ecdfs[[i]] %>%
    slice(floor(seq(1, nrow(.), length.out = min(nrow(.), 10^5)))) %>%   # sample on grid for plotting
    select(x, y)

  my_ecdf <- function(x) {
    # first known values, then values we want to evaluate --> fill with last observation
    rbind(cbind(df_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>%
      arrange(x) %>%
      setnafill(type = "locf") %>%
      filter(a > 0) %>%
      arrange(a) %>%
      pull(y)
  }

  xmin <- 0.7
  xmax <- 1.0
  ymin <- 0.0
  ymax <- 0.25

  my_breaks <- seq(0, 1, length.out = 9)

  my_hist <- ggplot(data.table(x = as.vector(my_models[[i]]))) +
    geom_histogram(aes(x = my_ecdf(x)), fill = "gray", col = "black", size = 0.2,
                   breaks = my_breaks) +
    theme_classic(base_size = 5.5) +
    theme(axis.line.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.background = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  inset_histograms <- layer(
    data = data.frame(Model = model_names[i], x = 0), stat = StatIdentity,
    position = PositionIdentity, geom = GeomCustomAnn, inherit.aes = TRUE,
    params = list(grob = ggplotGrob(my_hist), xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax))

  t_breaks_y <- my_ecdf(breaks_y[[i]])
  labels_y <- my_labeller(breaks_y[[i]])
  t_breaks_x <- my_ecdf(breaks_x[[i]])
  labels_x <- my_labeller(breaks_x[[i]])
  minor_breaks <- my_ecdf(minor_breaks)

  main_plot <- ggplot(filter(recal_models, Model == model_names[i]),
                      aes(x = my_ecdf(x))) +
    facet_wrap(~Model, nrow = 1) +
    geom_ribbon(aes(ymin = my_ecdf(lower), ymax = my_ecdf(upper), fill = Model),
                alpha = 0.33, show.legend = FALSE) +
    geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
                linetype = "dashed") +
    geom_step(aes(y = my_ecdf(x_rc), color = Model), size = 0.3, show.legend = FALSE) +
    scale_color_manual(values = model_colors) +
    scale_fill_manual(values = model_colors) +
    scale_x_continuous(breaks = t_breaks_x, labels = labels_x, minor_breaks = minor_breaks) +
    scale_y_continuous(breaks = t_breaks_y, labels = labels_y, minor_breaks = minor_breaks) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(NULL) +
    geom_text(data = filter(collect_stats, Model == model_names[i]),
              mapping = aes(x = 0.05, y = 0.75, label = label),
              size = 8 * 0.36, hjust = 0, vjust = 0) +
    my_theme +
    theme(aspect.ratio = 1)

  collect_recal_plots[[i]] <- main_plot + inset_histograms
}

combine <- grid.arrange(collect_recal_plots[[1]],
                        collect_recal_plots[[2]],
                        collect_recal_plots[[3]],
                        collect_recal_plots[[4]],
                        nrow = 2,
                        top = textGrob("Reliability Diagram",
                                       gp = gpar(fontsize = title_size)),
                        bottom = textGrob("Forecasted mean",
                                       gp = gpar(fontsize = 11)),
                        left = textGrob("Conditional mean", rot = 90,
                                       gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "all_Tilmann.pdf")
ggsave(file_path, width = 145, height = 160, unit = "mm", plot = combine)

# log log scale ----------------------------------------------------------------
main_plot <- ggplot(recal_models, aes(x = x)) +
  facet_wrap(~Model, nrow = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model),
              alpha = 0.33, show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_step(aes(y = x_rc, color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_log10(limits = c(0.05, 7.0)) + # for daily forecasts
  scale_y_log10(limits = c(0.05, 7.0)) +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram (loglog)") +
  my_theme +
  theme(aspect.ratio = 1)

file_path <- file.path(fpath, "daily_resres_log.pdf")
ggsave(file_path, width = 145, height = 160, unit = "mm", plot = main_plot)
