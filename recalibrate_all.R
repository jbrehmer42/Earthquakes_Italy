# write zero score as it will not change sums
my_s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(mean(score))
}

# recalibration routines =======================================================
# (1) resample residuals -------------------------------------------------------

# adapted from
# https://github.com/dwolffram/replication-ARSIA2023/blob/main/R/reliability_functions.R
reldiag <- function(x, y, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # we can compress step function by only storing first value and jumps!
  filter_jumps <- function(v) return(c(T, v[-1] - v[-length(v)] > 0))

  score <- my_s_pois
  # score <- function(x, y) mean(S_quad(x, y))

  x_rc <- monotone(y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  n_samples <- n_resamples + 1 # total number of samples including observed sample
  low <- floor(n_samples * (1 - region_level) / 2)
  up <- n_samples - low
  low <- pmax(low, 1)

  jumps <- filter_jumps(x_rc)
  collect_vals <- list(data.table(x = x[jumps], y = x_rc[jumps], I = "Fit"))

  res <- y - x
  mean_res <- mean(res)
  for (i in 2:n_samples) {
    # y <- round(pmax(0, x + sample(res, length(y), replace = TRUE)))
    y <- x + sample(res, length(y), replace = TRUE)
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    jumps <- filter_jumps(x_rc)
    collect_vals[[i]] <- data.table(x = x[jumps], y = pmax(0, x_rc[jumps] - mean_res),
                                    I = paste0("R", i))
  }
  # build joint data table with the collected values
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

# (2) adapt unconditional distribution as suggested by Tilmann -----------------
reldiag <- function(x, y, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # we can compress step function by only storing first value and jumps!
  filter_jumps <- function(v) return(c(T, v[-1] - v[-length(v)] > 0))

  score <- my_s_pois
  # score <- function(x, y) mean(S_quad(x, y))

  x_rc <- monotone(y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  n_samples <- n_resamples + 1 # total number of samples including observed sample
  low <- floor(n_samples * (1 - region_level) / 2)
  up <- n_samples - low
  low <- pmax(low, 1)

  jumps <- filter_jumps(x_rc)
  collect_vals <- list(data.table(x = x[jumps], y = x_rc[jumps], I = "Fit"))

  t <- table(y) / length(y)
  vals <- as.numeric(names(t))
  f_y <- as.numeric(t)
  eps <- sum(f_y[-1] * vals[-1]) / x - 1
  n_split <- 5      # as not everything fits into memory, split
  split <- ceiling(seq(1, length(y), length.out = n_split))
  cmp <- matrix(rep(cumsum(f_y[-length(f_y)]), split[2] - split[1] + 1),
                ncol = length(f_y) - 1, byrow = T)
  for (i in 2:n_samples) {
    for (j in 2:n_split) {
      i_vec <- split[j - 1]:split[j]
      u <- runif(length(i_vec), 0, 1 + eps[i_vec])
      y[i_vec] <- rowSums(u - eps[i_vec] >= cmp[1:pmin(length(i_vec), nrow(cmp)), ])
    }
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    jumps <- filter_jumps(x_rc)
    collect_vals[[i]] <- data.table(x = x[jumps], y = x_rc[jumps],
                                    I = paste0("R", i))
  }
  # build joint data table with the collected values
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

# all values
my_obs <- as.vector(obs)
my_models <- lapply(models, as.vector)
n_resamples <- 20
# daily values
my_obs <- rowSums(obs)
my_models <- lapply(models, rowSums)
n_resamples <- 5000
# only use every 7-th value, as forecasts span 7-day periods
mod <- 7      # use 1 to 7
seventh <- 0:floor(nrow(obs) / 7)
pick <- (7 * seventh + mod)[7 * seventh + mod <= nrow(obs)]
weekday <- lubridate::wday(times[mod], label = T, abbr = F)

my_obs <- as.vector(obs[pick, ])
my_models <- lapply(models, function(x) as.vector(x[pick, ]))
n_resamples <- 20


reldiag_cmp <- compiler::cmpfun(reldiag)

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  res <- reldiag_cmp(my_models[[i]], my_obs, n_resamples = n_resamples)
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
for (i in 1:length(models)) {
  t <- table(my_models[[i]])
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / length(my_models[[i]]),
                               M = model_names[i])
}

# (1) use mean empirical CDF ---------------------------------------------------

mean_ecdf <- do.call(rbind, col_ecdfs) %>%
  pivot_wider(id_cols = x, names_from = M, values_from = y, values_fill = NA) %>%
  arrange(x) %>%
  setnafill(type = "locf") %>%
  slice(floor(seq(1, nrow(.), length.out = 10^5))) %>%   # sample on grid for plotting
  transmute(x = x, y = rowMeans(cbind(LM, FMC, LG, SMA)))

my_ecdf <- function(x) {
  # first known values, then values we want to evaluate --> fill with last observation
  rbind(cbind(mean_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>%
    arrange(x) %>%
    setnafill(type = "locf") %>%
    filter(a > 0) %>%
    arrange(a) %>%
    pull(y)
}

# create inset-histograms
inset_histograms <- list()
for (i in 1:length(models)) {
  xmin <- 0.7
  xmax <- 1.0
  ymin <- 0.0
  ymax <- 0.25

  my_breaks <- seq(0, 1, length.out = 9)

  my_hist <- ggplot(data.table(x = my_models[[i]])) +
    geom_histogram(aes(x = my_ecdf(x)), fill = "gray", col = "black", size = 0.2,
                   breaks = my_breaks) +
    theme_classic(base_size = 5.5) +
    theme(axis.line.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.background = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  inset_histograms[[i]] <- layer(
    data = data.frame(Model = model_names[i], x = 0), stat = StatIdentity,
    position = PositionIdentity, geom = GeomCustomAnn, inherit.aes = TRUE,
    params = list(grob = ggplotGrob(my_hist), xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax))
}

breaks_y <- c(0, 10^c(-7, -6, -5, -4, 0))
t_breaks_y <- my_ecdf(breaks_y)
labels_y <- c("0", paste0("1e-", c(7, 6, 5, 4)), "1")
breaks_x <- c(0, 10^c(-6, -5, 0))
t_breaks_x <- my_ecdf(breaks_x)
labels_x <- c("0", paste0("1e-", c(6, 5)), "1")
minor_breaks <- my_ecdf(c(0, 10^(-10:0)))

main_plot <- ggplot(recal_models, aes(x = my_ecdf(x))) +
  facet_wrap(~Model, nrow = 1) +
  geom_ribbon(aes(ymin = my_ecdf(lower), ymax = my_ecdf(upper), fill = Model), alpha = 0.33,
              show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_step(aes(y = my_ecdf(x_rc), color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(breaks = t_breaks_x, labels = labels_x, minor_breaks = minor_breaks) +
  scale_y_continuous(breaks = t_breaks_y, labels = labels_y, minor_breaks = minor_breaks) +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram") +
  geom_text(data = collect_stats, mapping = aes(x = 0.05, y = 0.75, label = label),
            size = 8 * 0.36, hjust = 0, vjust = 0) +
  theme_bw() +
  theme(strip.background = element_blank(), aspect.ratio = 1)

combine_plots <- main_plot + inset_histograms

file_path <- file.path(fpath, "rc_probTrans_avg.pdf")
ggsave(file_path, width = 310, height = 90, unit = "mm", plot = combine_plots)

# use individual empirical CDFs ------------------------------------------------

collect_recal_plots <- list()

# for non-aggregated data
breaks_x <- list(c(0, 10^c(-5, 0)), c(0, 10^c(-8, -6, -5, 0)),
                 c(0, 10^c(-6, -5, 0)), c(0, 10^c(-6, -5, 0)))
breaks_y <- list(c(0, 10^c(-6, -5, -4, 0)), c(0, 10^c(-8, -7, -6, -5, -4, 0)),
                 c(0, 10^c(-6, -5, 0)), c(0, 10^c(-6, -5, -4, 0)))
my_labeller <- function(l) {
  labels <- character(length(l))
  labels[l > 0 & l < 1] <- paste0("1e", log(l[l > 0 & l < 1], base = 10))
  labels[l == 1] <- "1"
  labels[l == 0] <- "0"
  return(labels)
}


for (i in 1:length(models)) {
  mean_ecdf <- col_ecdfs[[i]] %>%
    slice(floor(seq(1, nrow(.), length.out = min(nrow(.), 10^5)))) %>%   # sample on grid for plotting
    select(x, y)

  my_ecdf <- function(x) {
    # first known values, then values we want to evaluate --> fill with last observation
    rbind(cbind(mean_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>%
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

  my_hist <- ggplot(data.table(x = my_models[[i]])) +
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
  minor_breaks <- my_ecdf(c(0, 10^(-10:0)))

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
    theme_bw() +
    theme(strip.background = element_blank(), aspect.ratio = 1)

  collect_recal_plots[[i]] <- main_plot + inset_histograms
}

combine <- grid.arrange(collect_recal_plots[[2]] + ylab("Conditional mean"),
                        collect_recal_plots[[3]],
                        collect_recal_plots[[1]], collect_recal_plots[[4]],
                        nrow = 1,
                        top = textGrob("Reliability diagram",
                                       gp = gpar(fontsize=15)),
                        bottom = textGrob("Forecasted mean",
                                       gp = gpar(fontsize=11)))

file_path <- file.path(fpath, "rc_probTrans_ind.pdf")
ggsave(file_path, width = 310, height = 90, unit = "mm", plot = combine)

# log log scale ----------------------------------------------------------------
main_plot <- ggplot(recal_models, aes(x = x)) +
  facet_wrap(~Model, nrow = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model),
              alpha = 0.33, show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_step(aes(y = x_rc, color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_log10(limits = c(0.05, 5.0)) +
  scale_y_log10(limits = c(0.05, 5.0)) +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram") +
  theme_bw() +
  theme(strip.background = element_blank(), aspect.ratio = 1)

file_path <- file.path(fpath, "rc_loglog_daily_resres_round_pmax.pdf")
ggsave(file_path, width = 310, height = 90, unit = "mm",
       plot = main_plot + ggtitle("Reliability Diagram (loglog)"))


# ==============================================================================
# resample from transformed residuals with empirical transform

# resample from transformed residuals
reldiag <- function(x, y, trans, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # we can compress step function by only storing first value and jumps!
  filter_jumps <- function(v) return(c(T, v[-1] - v[-length(v)] > 0))

  score <- my_s_pois
  # score <- function(x, y) mean(S_quad(x, y))

  x_rc <- monotone(y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  n_samples <- n_resamples + 1 # total number of samples including observed sample
  low <- floor(n_samples * (1 - region_level) / 2)
  up <- n_samples - low
  low <- pmax(low, 1)

  jumps <- filter_jumps(x_rc)
  collect_vals <- list(data.table(x = x[jumps], y = x_rc[jumps], I = "Fit"))

  res <- trans$transform(y) - trans$transform(x)
  mean_res <- mean(res)
  for (i in 2:n_samples) {
    y <- trans$inverse(trans$transform(x) + sample(res, length(y), replace = TRUE))
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    jumps <- filter_jumps(x_rc)
    collect_vals[[i]] <- data.table(x = x[jumps], y = pmax(0, x_rc[jumps]),
                                    I = paste0("R", i))
  }
  # build joint data table with the collected values
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

reldiag_cmp <- compiler::cmpfun(reldiag)

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  mean_ecdf <- col_ecdfs[[i]] %>%
    slice(floor(seq(1, nrow(.), length.out = min(nrow(.), 10^5)))) %>%   # sample on grid for plotting
    select(x, y)

  my_ecdf <- function(x) {
    # first known values, then values we want to evaluate --> fill with last observation
    rbind(cbind(mean_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>%
      arrange(x) %>%
      setnafill(type = "locf") %>%
      filter(a > 0) %>%
      arrange(a) %>%
      pull(y)
  }
  my_equa <- function(y) {
  # first values for evaluation, then known values --> fill with next obersvation
  rbind(data.table(x = NA, y = y, a = 1:length(y)), cbind(mean_ecdf, a = -1.0)) %>%
    arrange(y) %>%
    setnafill(type = "nocb") %>%
    filter(a > 0) %>%
    arrange(y) %>%
    pull(y)
  }

  trans <- trans_new("ecdf", my_ecdf, my_equa)

  res <- reldiag_cmp(as.vector(models[[i]]), as.vector(obs), trans, n_resamples = 1000)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.2e", res$stats[1, ]),
                        collapse = "\n"))
  )
}
