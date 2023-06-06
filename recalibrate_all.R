library(monotone)     # for fast isotonic regression
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(grid)
# for fast handling of large data.frames, therefore don't only use the data type
# data.table but also its native functions
library(data.table)

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
  # create data table -> sort -> add isotonic fit
  dt <- data.table(x = as.vector(x), y = as.vector(y))[
    order(x, -y)][, x_rc := monotone(y)]
  return(dt)
}

fit_isotonic_by_cell <- function(x, y) {
  # fit isotonic regression for each cell separately

  # check whether all observations are 0 or less (in resampling we can have that)
  zero_cols <- apply(y, 2, function(col) all(col <= 0))

  recal_by_cell <- do.call(rbind, lapply(1:ncol(x), function(col) {
    if (zero_cols[col]) {
      return(data.table(x = x[, col], y = y[, col], x_rc = 0))   # recalibration is just 0's
    } else {
      return(fit_isotonic_all(x[, col], y[, col]))
    }
  }))
  # we have for the same forecast value, possibly #cells different recalibrated values
  # --> estimate conditional averages with isotonic regression (i.e., we are estimating
  # conditional averages of conditional averages ~iterated conditional expectation)
  # ! normal averaging destroys isonoticity as we could have ((1,2), (2,3)) and
  #  ((2, 0), (3,1)) -pointwise-avg-> ((1,2), (2,1.5), (3,1))
  recal_by_cell <- recal_by_cell[order(x, -x_rc),][, x_rc := monotone(x_rc)]
  return(recal_by_cell)
}

resample_residuals <- function(x, y, B, get_iso_fit) {
  collect_vals <- list()

  res <- y - x
  mean_res <- mean(res)
  for (i in 1:B) {
    y <- x + matrix(sample(as.vector(res), prod(dim(y)), replace = TRUE),
                    nrow = nrow(y), ncol = ncol(y))

    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = pmax(0, x_rc - mean_res), I = paste0("R", i))]
  }
  return(collect_vals)
}

resample_daily_residuals <- function(x, y, B, get_iso_fit) {
  # resample rows, i.e. each acknowledge spatial distribution of errors
  collect_vals <- list()

  res <- y - x
  # by subtracting mean residual, accorind forecast is unconditionally and hence
  # conditionally (we assume forecast and residuals are independent) calibrated
  res <- res - matrix(rep(colMeans(res), nrow(res)), nrow(res), byrow = T)
  for (i in 1:B) {
    y <- x + res[sample(1:nrow(res), nrow(res), replace = TRUE), ]
    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = pmax(0, x_rc), I = paste0("R", i))]
  }
  return(collect_vals)
}

resample_cell_residuals <- function(x, y, B, get_iso_fit) {
  # resample, but only use cell-specific residuals
  collect_vals <- list()

  res <- y - x
  # by subtracting mean residual, accorind forecast is unconditionally and hence
  # conditionally (we assume forecast and residuals are independent) calibrated
  res <- res - matrix(rep(colMeans(res), nrow(res)), nrow(res), byrow = T)
  for (i in 1:B) {
    y <- x + apply(res, 2, function(col) sample(col, length(col), replace = T))
    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = pmax(0, x_rc), I = paste0("R", i))]
  }
  return(collect_vals)
}

resample_trans_residuals <- function(x, y, B, get_iso_fit) {
  # resample transformed resiudals, transformation has to be defined in global scope
  collect_vals <- list()

  res <- trans$transform(y) - trans$transform(x)
  for (i in 1:B) {
    y <- trans$inverse(
      trans$transform(x) + matrix(sample(as.vector(res), prod(dim(y)), replace = TRUE),
                                  nrow = nrow(y), ncol = ncol(y))
    )

    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = pmax(0, x_rc - mean_res), I = paste0("R", i))]
  }
  return(collect_vals)
}

resample_adapt_ucond_distr <- function(x, y, B, get_iso_fit) {
  # Adapt unconditional distribution of observed earthquakes y to construct a
  # mean calibrated probabilitstic forecast for each mean forecast x.
  # We achieve this by adapting the weight of the 0 of the uncondtional distribution,
  # i.e., if x > mean(ucond_y), we make the weight of 0 smaller and if
  # x < mean(ucond_y), we make the weight larger, until we have equality, so that
  # the mean of our adapted distribution just corresponds to the mean forecast.
  # Having this mean calibrated forecast, we can resample y by drawing from
  # the adapted distribution.
  # This method can not be used with fit_isotonic_by_cell
  # One Problem of this method is that for small forecasts x the according ditribution
  # puts virtally all mass in 0, i.e., resampling returns only 0's. If we would put
  # higher earthquake numbers firts to zero, then maybe the 1 would have some
  # non-negligible mass, so that resampling does not only return 1's
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
    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = x_rc, I = paste0("R", i))]
  }
  return(collect_vals)
}

resample_adapt_ucond_distr_by_cell <- function(x, y, B, get_iso_fit) {
  # Adaption of resample_adapt_ucond_distr for fit_isotonic_by_cell

  collect_vals <- list()

  n_rows <- nrow(y)
  n_cols <- ncol(y)
  zero_cols <- colSums(y) == 0
  t <- lapply(1: n_cols, function(col) {
    if (zero_cols[col]) {                                 # we observed always 0 earthquakes
      return(c("0" = 1 - 1 / n_rows, "1" = 1 / n_rows))   # add one 1, so that we can raise
    } else {
      return(table(y[, col]) / n_rows)
    }
  })
  vals <- lapply(t, function(tab) as.numeric(names(tab)))   # values that are assumed
  f_y <- lapply(t, function(tab) as.numeric(tab))           # with according probabilities
  rm(t, zero_cols)
  eps <- lapply(1:n_cols, function(col) sum(f_y[[col]][-1] * vals[[col]][-1]) / x[, col] - 1)

  for (i in 1:B) {
    y <- do.call(cbind, lapply(1: n_cols, function(col) {
      u <- runif(n_rows, 0, 1 + eps[[col]])
      my_f_y <- f_y[[col]]
      return(
        rowSums(u - eps[[col]] >= matrix(rep(cumsum(my_f_y[-length(my_f_y)]), n_rows),
                                         nrow = n_rows, byrow = T)))
    }))
    collect_vals[[i]] <- get_iso_fit(x, y)[
      filter_jumps(x_rc), .(x = x, x_rc = x_rc, I = paste0("R", i))]
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

  collect_vals <- list(my_iso_fit[filter_jumps(x_rc), .(x = x, x_rc = x_rc, I = "Fit")])
  # bootstrap fits for consistency bands
  bootstrap_fits <- res_for_cons(x, y, n_resamples, get_iso_fit)
  collect_vals <- append(collect_vals, bootstrap_fits)

  sum_stats <- function(row) {
    sorted <- sort(as.numeric(row[1, -"Fit"]))
    return(list(x_rc = row[1, Fit], lower = sorted[low], upper = sorted[up]))
  }
  # build data table with the collected values on a joint grid
  results <- do.call(rbind, collect_vals)
  # pivot wider and sort (to use step function logic to fill missing vals)
  results <- dcast(results, x ~ I, value.var = "x_rc", fill = NA)[order(x)]
  # use last available observation to fill NA (step fcn!), and calc statistics on each row
  results <- setnafill(results, type = "locf")[, sum_stats(.SD), keyby=x]

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

# runtimes per fit: reldiag fits n_resamples + 1 isotonic regressions, and
# resamples therefore n_resamples times
#
# rough runtime durations with n_resamples = 5 averaged over the models
#       (1)   (2)   (3)   (4)   (5)
# (A)   100s  60s   45s   60s   2min
# (B)   -     -     -     -     5min

# * for FMC it took 40min, probably since it issues the most different forecasts
# (2.4 * 10^6) and in (B) we group_by forecast to reconcile forecasts from
# different cells

# (B1) is not possible as resampling discards spatial structure,
# (B2, B3, B4) is computationally not possible, as we are fitting #cells isotonic
# fits, and residual resampling does not produce any trivial (only 0's) cells
# unlike (B5)

pick_iso <- fit_isotonic_all        # (A)
pick_iso <- fit_isotonic_by_cell    # (B)

pick_res <- resample_adapt_ucond_distr  # (1)
pick_res <- resample_residuals          # (2)
pick_res <- resample_daily_residuals    # (3)
pick_res <- resample_cell_residuals     # (4)
pick_res <- resample_adapt_ucond_distr_by_cell  # (5)

reldiag_cmp <- compiler::cmpfun(function(x, y) {
  reldiag(x, y, n_resamples, get_iso_fit = pick_iso, res_for_cons = pick_res)
})

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  print(Sys.time())
  cat("-", model_names[i], "\n")
  res <- reldiag_cmp(my_models[[i]], my_obs)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.3f", res$stats[1, ]),
                        collapse = "\n"))
  )
  print(Sys.time())
}

# or just load already recalibrated values
recal_models <- read.csv("./../tmp_results/recal_models_all-Tilmann-100-new.csv", row.names = 1)
collect_stats <- read.csv("./../tmp_results/collect_stats_new.csv", row.names = 1)

# plot calibration curve under certain data transformations ====================
# use empirical CDF transform
col_ecdfs <- list()
for (i in 1:length(my_models)) {
  t <- table(my_models[[i]])
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / prod(dim(my_models[[i]])),
                               M = model_names[i])
}

# or just loaded already computed values
for (i in 1:length(models)) {
  col_ecdfs[[i]] <- read.csv(paste0("./../tmp_results/ecdf_", model_names[i], ".csv"), row.names = 1)
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
              mapping = aes(x = 0.02, y = 0.72, label = label),
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

file_path <- file.path(fpath, "all_Til_meanTrans_100.pdf")
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
