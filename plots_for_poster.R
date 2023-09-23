library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

library(gridExtra)        # grid.arrange : put ggplots next to each other
library(monotone)         # for fast isotonic regression
library(sf)               # st_as_sf : convert data.frame to geographic format sf
library(rnaturalearth)    # ne_countries - to load country data
library(scales)           # trans_new - custom data transformation
library(showtext)         # use custom fonts

source("data_prep.R")
source("functions_eval.R")

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")

font_add("roboto", "./../Poster/fonts/Roboto/Roboto-Light.ttf")
font_add("roboto_bold", "./../Poster/fonts/Roboto/Roboto-Medium.ttf")
font_add("MiriamLibre", "./../Poster/fonts/MiriamLibre/MiriamLibre-Regular.ttf")
showtext_auto()

fpath <- "./figures_pos"

base_size <- 28 / 1.2
label_size <- base_size * 0.8

my_theme <- list(
  theme_bw(base_size = base_size) +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05),
        legend.title = element_text(family = "roboto", size = label_size),
        plot.title = element_text(family = "roboto"),
        axis.title = element_text(family = "roboto"),
        axis.text = element_text(family = "MiriamLibre"),
        legend.text = element_text(family = "MiriamLibre"))
)

################################################################################
# Visualize spatial distribution of earthquakes in training data
################################################################################

lon_lim <- range(cells$LON)
lat_lim <- range(cells$LAT)

# RAM intensive: PyCharms takes up to 3 GB
europe <- ne_countries(scale = "medium", returnclass = "sf")

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
                       breaks = 10^(-6:-2), labels = c("1e-6", "", "1e-4", "", "1e-2"),
                       trans = "log10", option = "magma",
                       guide = guide_colorbar(barheight = unit(40, "mm"))) +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  ggtitle(paste0("LM Model (",  ymd(times[i_time]), ")")) +
  my_theme +
  theme(legend.position = "right")

file_path <- file.path(fpath, "Poster_Fig1.pdf")
ggsave(file_path, width = 180, height = 140, unit = "mm", plot = one_pred)

rm(events_by_cell, pred_one_day, one_pred)

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

anno1_text <- data.frame(text = "FMC/\nLG/SMA\npreferred", x = 28.3, y = 45, Model = "SMA")
anno2_text <- data.frame(text = "LM\npreferred", x = 28.3, y = 41.5, Model = "SMA")

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
                       guide = guide_colorbar(barheight = unit(50, "mm"))) +
  ggtitle("Mean Poisson Score Difference") +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  geom_text(data = anno1_text, mapping = aes(x = x, y = y, label = text),
           size = label_size / .pt * 0.8, color = "#057ffa", family = "roboto") +
  geom_text(data = anno2_text, mapping = aes(x = x, y = y, label = text),
           size = label_size / .pt * 0.8, color = "#f51818", family = "roboto") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(5.5, 5.5, 3, 3),
        plot.margin = margin(5.5, 50, 0, 5.5),
        plot.title = element_text(family = "roboto", margin = margin(0, 0, 2, 0)))

file_path <- file.path(fpath, "Poster_Fig2.pdf")
ggsave(file_path, width = 370, height = 145, unit = "mm", plot = spat_plot)

# and now temporally

scores_t <- do.call(cbind, lapply(models, function(X) rowMeans(S_pois(X, obs))))
scores_t <- data.frame(scores_t) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(scores_t) <- c(model_names, "X", "earthquake")

diff_scores_t <- scores_t %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

d2 <- 10^7
my_trans2 <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d2  + 1),
  function(y) sign(y) / d2 * (exp(abs(y)) - 1)
)
my_breaks_t <- c(-10^(c(-2, -4, -6)), 0, 10^(c(-6, -4)))
minor_breaks_t <- c(-10^(-2:-7), 0, 10^(-7:-3))
my_labels_t <- c(paste("-1e-", c(2, 4, 6)), "0", paste("1e-", c(6, 4)))

point_size <- 2.5
point_alpha <- 0.4

temp_plot <- ggplot(diff_scores_t) +
  # plot points going with no earthquakes first!
  # (that is why we separated geom point in two to define order of drawing groups)
  geom_point(data = filter(diff_scores_t, !earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  # add black borders to triangle used for earthquake days
  geom_point(data = filter(diff_scores_t, earthquake), aes(x = X, y = value),
             color = "black", shape = 2, size = point_size, stroke = 0.5,
             alpha = point_alpha) +
  # plot triangles
  geom_point(data = filter(diff_scores_t, earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores_t$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores_t))) +
  scale_y_continuous(trans = my_trans2, breaks = my_breaks_t, labels = my_labels_t,
                     minor_breaks = minor_breaks_t) +
  scale_color_manual(name = "LM vs.", values = model_colors, breaks = ana_models,
                     guide = guide_legend(order = 1, override.aes = list(alpha = 1))) +
  scale_shape_manual(name = "Earthquake\ncount", labels = c("FALSE" = "= 0", "TRUE" = "> 0"),
                     values = c("FALSE" = 16, "TRUE" = 17)) +
  annotate(geom = "text", label = "FMC/LG/SMA preferred", x = 50, y = 0.001,
           size = label_size / .pt * 0.8, hjust = 0) +
  annotate(geom = "text", label = "LM preferred", x = 50, y = -0.005,
           size = label_size / .pt * 0.8, hjust = 0) +
  xlab(NULL) +
  ylab("Score difference") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "right", legend.text = element_text(size = label_size),
        plot.margin = margin(5.5, 35, 5.5, 5.5))

file_path <- file.path(fpath, "Poster_Fig3.pdf")
ggsave(file_path, width = 370, height = 105, unit = "mm", plot = temp_plot)

rm(scores, scores_t, diff_scores, diff_scores_t, spat_plot, temp_plot)

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

murphy_diag <- data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  scale_color_manual(name = NULL, values = model_colors) +
  xlab(expression(paste("log", theta))) +
  ylab("Elementary score") +
  ggtitle("Murphy Diagram") +
  my_theme +
  theme(legend.position = c(0.01, 1.03), legend.justification = c(0, 1),
        legend.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig4.pdf")
ggsave(file_path, width = 150, height = 135, unit = "mm", plot = murphy_diag)

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

murphy_diff_diag <- data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(theta, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models)) +
  xlab(expression(paste("log", theta))) +
  ylab("Elementary score") +
  ggtitle("Murphy Difference Diagram") +
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
        legend.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig4_Diff.pdf")
ggsave(file_path, width = 150, height = 140, unit = "mm", plot = murphy_diff_diag)

rm(murphy_df, murphy_diag, murphy_diff_diag)

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

collect_stats <- mutate(collect_stats,
                         label = paste0("Score ", sprintf("%0.3f", Score),
                                        "\nMCB ", sprintf("%0.3f", MCB),
                                        "\nDSC ", sprintf("%0.3f", DSC),
                                        "\nUNC ", sprintf("%0.3f", UNC)))

main_plot <- ggplot(recal_models, aes(x = x)) +
  facet_wrap(~Model, nrow = 2) +
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
  geom_text(data = collect_stats, mapping = aes(x = 10^(-9), y = 0.001, label = label),
            size = label_size / .pt * 0.9, hjust = 0, vjust = 0, family = "MiriamLibre") +
  my_theme +
  theme(strip.background = element_blank(), aspect.ratio = 1,
        panel.spacing.x = unit(10, "mm"), plot.background = element_blank())

combine_plots <- main_plot + inset_histograms

file_path <- file.path(fpath, "Poster_Fig5.pdf")
ggsave(file_path, width = 270, height = 245, unit = "mm", plot = combine_plots)

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

rm(recal_models, collect_stats, main_plot, inset_histograms)
