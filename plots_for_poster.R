library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

library(gridExtra)        # grid.arrange : put ggplots next to each other
library(monotone)         # for fast isotonic regression
library(sf)               # st_as_sf : convert data.frame to geographic format sf
library(rnaturalearth)    # ne_countries - to load country data
library(scales)           # trans_new - custom data transformation

source("data_prep.R")
source("functions_eval.R")

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")

fpath <- "./figures_pos"

base_size <- 20 / 1.2
label_size <- base_size * 0.8

my_theme <- list(
  theme_bw(base_size = base_size) +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05),
        plot.title = element_text(size = 20))
)

################################################################################
# Visualize spatial distribution of earthquakes in training data
################################################################################

lon_lim <- range(cells$LON)
lat_lim <- range(cells$LAT)

# RAM intensive: PyCharms takes up to 3 GB
europe <- ne_countries(scale = "medium", returnclass = "sf")

cells_and_events <- events %>%
  group_by(N) %>%
  summarise(Count = n(), .groups = "drop") %>%
  right_join(cells, by = "N")

spat_distr <- ggplot() +
  geom_sf(data = europe, fill = NA, color = "gray", alpha = 0.4, size = 0.2) +
  geom_tile(data = cells_and_events, aes(x = LON, y = LAT, fill = Count), color = "#f5a536") +
  geom_sf(data = filter(europe, name == "Italy"), fill = NA, color = "black",
          alpha = 0.4, size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_viridis_c(name = "Count", na.value = NA) +
  ggtitle("M4+ Earthquakes 2005-2020") +
  my_theme +
  theme(legend.title = element_text(size = 8))

mag_distr <- events %>%
  mutate(bin_mag = cut_width(MAG, 0.5, boundary = 4.0)) %>%
  group_by(bin_mag) %>%
  summarise(Count = n(), .groups = "drop") %>%
  ggplot() +
  geom_col(aes(x = Count, y = bin_mag), fill = "#100030") +
  scale_x_log10(breaks = c(1, 10, 100), minor_breaks = c(1:10, (2:10) * 10)) +
  ylab("Magnitude") +
  ggtitle("Magnitude Distribution") +
  my_theme +
  theme(aspect.ratio = 1.2)

# now look at one forecasts spatially
i_time <- 1449  # 3 days later Mag 6.1 earthquake
pred_one_day <- data.frame(LM = models[[which(model_names == "LM")]][i_time, ]) %>%
  mutate(N = 1:nrow(.)) %>%
  right_join(cells, by = "N") %>%
  select(LON, LAT, LM)

events_by_cell <- events %>%
  filter(TS >= times[i_time], TS < times[i_time] + days(7)) %>%
  group_by(N) %>%
  summarise(Count = n(), .groups = "drop") %>%
  left_join(cells, by = "N") %>%
  select(LON, LAT, Count)

one_pred <- ggplot() +
  geom_sf(data = europe, color = "gray", alpha = 0.4,
          size = 0.2, fill = NA) +
  geom_tile(data = pred_one_day, aes(x = LON, y = LAT, fill = LM), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", alpha = 0.4,
          size = 0.2, fill = NA) +
  geom_tile(data = events_by_cell, aes(x = LON, y = LAT, color = "Obs. earthquakes"), fill = NA) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_viridis_c(name = "Pred.\nmean",
                       breaks = 10^(c(-6, -4, -2)), labels = paste0("e", c(-6, -4, -2)),
                       trans = "log10", option = "magma") +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  ggtitle(paste0("LM Model (",  ymd(times[i_time]), ")")) +
  my_theme +
  theme(legend.position = "right", legend.title = element_text(size = label_size))

file_path <- file.path(fpath, "Poster_Fig1.pdf")
ggsave(file_path, width = 110, height = 110, unit = "mm", plot = one_pred)

rm(cells_and_events, events_by_cell, pred_one_day, spat_distr, mag_distr,
   one_pred)

################################################################################
# Visualize Poisson Scores over time and score differences spatially
################################################################################

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

scores <- do.call(cbind, lapply(models, function(X) colMeans(S_pois(X, obs))))
scores <- data.frame(scores)
colnames(scores) <- model_names

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  mutate(LON = cells$LON, LAT = cells$LAT) %>%
  select(LON, LAT, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")
diff_scores$Model <- paste(cmp_model, "vs.", diff_scores$Model)

my_colors <- c("#800303", "#f51818", "#ffffff", "#057ffa")
limits <- range(diff_scores$value)
my_breaks <- c(-0.01, -0.001, -0.0001, 0, 0.0001, 0.001)
my_labels <- c("-1e-2", "-1e-3", "-1e-4", " 0", " 1e-4", " 1e-3")

# need log transform for positive and negative values (see ?modulus_trans)
# but need to scale with d to get sufficient resolution
d1 <- 10^5
my_trans <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d1  + 1),
  function(y) sign(y) / d1 * (exp(abs(y)) - 1)
)
col_breaks <- my_trans$transform(c(limits[1], -limits[2], 0, limits[2]))

eq_loc <- events %>%
  group_by(N) %>%
  summarise(Count = n(), .groups = "drop")  %>%
  left_join(cells, by = "N") %>%
  select(LON, LAT, Count)

spat_plot <- ggplot() +
  facet_grid(~Model) +
  geom_tile(data = diff_scores,
            aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_tile(data = eq_loc, aes(x = LON, y = LAT, color = "Obs. earthquakes"),
            fill = NA) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", fill = NA,
          size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE, clip = 'off') +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "Score\ndifference",
                       trans = my_trans, colors = my_colors, values = rescale(col_breaks),
                       breaks = my_breaks, labels = my_labels,
                       guide = guide_colorbar(barheight = unit(35, "mm"))) +
  ggtitle("Average Poisson Score Differences") +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  annotate(geom = "text", label = "FMC/LG/SMA\ndominates", x = 30, y = 44.5,
           size = label_size / .pt * 0.8, color = "#057ffa") +
  annotate(geom = "text", label = "LM\ndominates", x = 30, y = 42,
           size = label_size / .pt * 0.8, color = "#f51818") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = label_size),
        plot.margin = margin(5.5, 25.5, 0, 5.5))

# and now temporally

scores <- do.call(cbind, lapply(models, function(X) rowMeans(S_pois(X, obs))))
scores <- data.frame(scores) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(scores) <- c(model_names, "X", "earthquake")

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

d2 <- 10^7
my_trans2 <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d2  + 1),
  function(y) sign(y) / d2 * (exp(abs(y)) - 1)
)
my_breaks <- c(-10^(c(-2, -4, -6)), 0, 10^(c(-6, -4)))
minor_breaks <- c(-10^(-2:-7), 0, 10^(-7:-3))
my_labels <- c(paste("-1e-", c(2, 4, 6)), "0", paste("1e-", c(6, 4)))

temp_plot <- ggplot(diff_scores) +
  geom_vline(data = filter(diff_scores, earthquake > 0),
             aes(xintercept = X, linetype = "Obs. earthquakes"),
             alpha = 0.2, color = "gray", size = 0.3) +
  geom_point(aes(x = X, y = value, color = Model), size = 0.75, alpha = 0.4) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores))) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models),
                     guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
  scale_linetype_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = 1),
                        labels = "",
                        guide = guide_legend(override.aes = list(alpha = 1, size = 0.5),
                                             order = 2)) +
  scale_y_continuous(trans = my_trans2, breaks = my_breaks, labels = my_labels,
                     minor_breaks = minor_breaks) +
  annotate(geom = "text", label = "FMC/LG/SMA dominates", x = 50, y = 0.001,
           size = label_size / .pt * 0.8, hjust = 0) +
  annotate(geom = "text", label = "LM dominates", x = 50, y = -0.005,
           size = label_size / .pt * 0.8, hjust = 0) +
  xlab(NULL) +
  ylab("log-transformed score") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "right", legend.title = element_text(size = label_size),
        legend.text = element_text(size = label_size),
        plot.margin = margin(5.5, 42.7, 5.5, 16.5))

combine_plots <- grid.arrange(spat_plot, temp_plot, nrow = 2, heights = c(11, 6) / 17)

file_path <- file.path(fpath, "Poster_Fig5.pdf")
ggsave(file_path, width = 310, height = 170, unit = "mm", plot = combine_plots)

rm(scores, diff_scores, combine_plots, spat_plot, temp_plot)

################################################################################
# Visualize reliability diagram
################################################################################

# write zero score as it will not change sums
my_s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(mean(score))
}

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

reldiag_cmp <- compiler::cmpfun(reldiag)

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  res <- reldiag_cmp(as.vector(models[[i]]), as.vector(obs), n_resamples = 20)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.2e", res$stats[1, ]),
                        collapse = "\n"))
  )
}

recal_models <- read.csv("./../tmp_results/recal_models_all-Til-100-5.csv",
                         row.names = 1) %>%
  filter(Model %in% model_names)
collect_stats <- read.csv("./../tmp_results/collect_stats_5.csv",
                          row.names = 1) %>%
  filter(Model %in% model_names)

# now plot

d <- 10^9
my_trans <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d  + 1),
  function(y) sign(y) / d * (exp(abs(y)) - 1)
)

plot_min <- my_trans$transform(min(c(recal_models$x, recal_models$lower)))
plot_max <- my_trans$transform(max(c(recal_models$x, recal_models$upper)))
hist_breaks <- seq(plot_min, plot_max, length.out = 9)

# create inset-histograms
inset_histograms <- list()
for (i in 1:length(models)) {
  xmin <- my_trans$transform(10^-4)
  xmax <- my_trans$transform(0.5)
  ymin <- my_trans$transform(10^-9)
  ymax <- my_trans$transform(10^-6)

  my_hist <- ggplot(data.table(x = my_trans$transform(as.vector(models[[i]])))) +
    geom_histogram(aes(x = x), fill = "gray", col = "black", size = 0.2,
                   breaks = hist_breaks) +
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

my_breaks <- c(0, 10^c(-8, -6, -4, -2, 0))
my_labels <- c("0", paste0("1e", c(-8, -6, -4, -2)), "1")
minor_breaks <- c(0, 10^(-10:0))

main_plot <- ggplot(recal_models, aes(x = x)) +
  facet_wrap(~Model, nrow = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.33,
              show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_step(aes(y = x_rc, color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_continuous(trans = my_trans, breaks = my_breaks,
                     minor_breaks = minor_breaks, labels = my_labels) +
  scale_y_continuous(trans = my_trans, breaks = my_breaks,
                     minor_breaks = minor_breaks, labels = my_labels) +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram") +
  geom_text(data = collect_stats, mapping = aes(x = 10^(-9), y = 0.002, label = label),
            size = label_size / .pt * 0.8, hjust = 0, vjust = 0) +
  my_theme +
  theme(strip.background = element_blank(), aspect.ratio = 1)

combine_plots <- main_plot + inset_histograms

file_path <- file.path(fpath, "Poster_Fig6.pdf")
ggsave(file_path, width = 310, height = 100, unit = "mm", plot = combine_plots)

# for daily forecasts comparison with result from Jonas, use quadratic scoring fcn!

recal_models <- data.table()
collect_stats <- data.table()

for (i in 1:length(models)) {
  res <- reldiag_cmp(rowSums(models[[i]]), rowSums(obs), n_resamples = 999)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.2e", res$stats[1, ]),
                        collapse = "\n"))
  )
}

ggplot(recal_models, aes(x = x)) +
  facet_wrap(~factor(Model, ordered = TRUE, levels = c("LM", "FMC", "LG", "SMA")),
                     nrow = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.33,
              show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_line(aes(y = x_rc, color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.05, NA)) +
  geom_text(data = collect_stats, mapping = aes(x = 0.2, y = 2.0, label = label),
            size = 8 * 0.36, hjust = 0, vjust = 0) +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram") +
  my_theme +
  theme(aspect.ratio = 1)

ggsave("./../test/send2/reliability_daily.pdf", width = 150, height = 180,
       unit = "mm")

# recalibrate single cells: investigate cells with many earthquakes and cells with no
# earthquakes

y_ord <- order(colSums(obs))

pick <- y_ord[floor(seq(1, length(y_ord), length.out = 5))][5]

recal_models <- data.table()
collect_stats <- data.table()

for (i in 1:length(models)) {
  res <- reldiag_cmp(models[[i]][, pick], obs[, pick], n_resamples = 999)
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i], res$stats,
          label = paste(names(res$stats), c("", " ", " ", " "),
                        sprintf("%.2e", res$stats[1, ]),
                        collapse = "\n"))
  )
}

ggplot(recal_models, aes(x = x)) +
  facet_wrap(~factor(Model, ordered = TRUE, levels = c("LM", "FMC", "LG", "SMA")),
                     nrow = 2, scales = "free") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.33,
              show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_line(aes(y = x_rc, color = Model), size = 0.3, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram") +
  my_theme +
  theme(aspect.ratio = 1)

rm(recal_models, collect_stats)

################################################################################
# Visualize Murphy Diagram
################################################################################

S_elem <- compiler::cmpfun(S_theta)

n_theta <- 100
log_grid <- seq(-24, 4, len = n_theta)
grd <- exp(log_grid)

# use for loops, as otherwise memory demand is too high
murphy_df <- matrix(0, nrow = length(grd), ncol = length(models))
for (m in 1:length(models)) {
  print(m)
  for (t in 1:n_theta) {
    cat("*")
    murphy_df[t, m] <- S_elem(models[[m]], obs, grd[t])
  }
}
colnames(murphy_df) <- model_names

# S_theta sums, but we want averages
murphy_df <- murphy_df / prod(dim(obs))

murphy_df <- read.csv("./../tmp_results/murphy_df.csv", row.names = 1) %>%
  select(all_of(model_names))

data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  scale_color_manual(name = NULL, values = model_colors) +
  xlab(expression(paste("Threshold log", theta))) +
  ylab("Elementary score") +
  ggtitle("Murphy Diagram") +
  my_theme +
  theme(legend.position = c(0.01, 1.03), legend.justification = c(0, 1),
        legend.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig7.pdf")
ggsave(file_path, width = 130, height = 80, unit = "mm")

data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(theta, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models)) +
  xlab(expression(paste("Threshold log", theta))) +
  ylab("Elementary score") +
  ggtitle("Murphy Difference Diagram") +
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
        legend.text = element_text(size = 8),
        legend.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig7_Diff.pdf")
ggsave(file_path, width = 110, height = 75, unit = "mm")

rm(murphy_df)
