library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

library(gridExtra)        # grid.arrange : put ggplots next to each other
library(monotone)         # for fast isotonic regression
library(sf)               # st_as_sf : convert data.frame to geographic format sf
library(rnaturalearth)    # ne_countries - to load country data
library(scales)           # seq_gradient_pal - custom color gradients
                          # trans_new - custom data transformation

source("data_prep.R")
source("functions_eval.R")

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")

fpath <- "./figures"

my_theme <- list(
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05))
)

################################################################################
# Visualize spatial distribution of earthquakes in training data
################################################################################

lon_lim <- range(cells$LON)
lat_lim <- range(cells$LAT)

# RAM intensive: PyCharms takes up to 3 GB
europe <- ne_countries(scale = "medium", returnclass = "sf")

m_low <- "#5e0306"
m_high <- "#fc474d"
mag_colors <- seq_gradient_pal(m_low, m_high, "Lab")(seq(0,1,length.out=6))

# scale size maps 4 to area = 4, but we want to scale it a bit down
size_trans <- trans_new(
  "linearScaleDown", function(x) x - 2, function(x) x + 2
)
new_range <- range(size_trans$transform(events$MAG))

ggplot() +
  # plot has a blue background: first draw countries white again
  geom_sf(data = europe, fill = "white", size = 0.2) +
  # before painting them in light green
  geom_sf(data = europe, fill = "#04b50a", color = "black", alpha = 0.4, size = 0.2) +
  geom_tile(data = cells, aes(x = LON, y = LAT), color = "#f5a536", fill = NA) +
  geom_sf(data = st_as_sf(events, coords = c("LON", "LAT"), crs = 4326),
          aes(size = MAG, color = MAG), alpha = 0.6, shape = 13) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = NULL) +
  scale_color_gradient(low = m_low, high = m_high) +
  scale_size(range = new_range, trans = size_trans) +
  guides(color = "none", size = guide_legend(title = "Magnitude",
                                             override.aes = list(color = mag_colors))) +
  ggtitle("M4+ Earthquakes") +
  theme_bw() +
  my_theme +
  theme(panel.background = element_rect(fill = rgb(0.227, 0.643, 0.941, 0.4)))

file_path <- file.path(fpath, "Poster_Fig1.pdf")
ggsave(file_path, width = 140, height = 140, unit = "mm")

################################################################################
# Visualize predictions over time and over grid cells
################################################################################

# first look at forecasts temporally
pred_by_day <- data.frame(do.call(cbind, lapply(models, rowSums))) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(pred_by_day) <- c(model_names, "X", "earthquake")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0
eq_shape <- 3   # shape of marks for earthquakes

pred_by_day %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  geom_point(data = filter(pred_by_day, earthquake),
             aes(x = X, y = 0.15, color = "Observed earthquakes"),
             size = 0.7, stroke = 0.4, shape = eq_shape) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(
    name = NULL,
    values = c(model_colors, "Observed earthquakes" = "black"),
    guide = guide_legend(override.aes = list(linetype = c(rep(1, 4), 0),
                                             shape = c(rep(NA, 4), eq_shape)), nrow = 4)
  ) +
  xlab("") +
  ylab("Predicted mean") +
  ggtitle("Predicted Mean Number of Events for all of Italy") +
  theme_bw() +
  my_theme +
  theme(legend.justification = c(0, 1), legend.position = c(0.01, 1.02),
        legend.direction = "vertical",
        legend.key.size = unit(0.5, "lines"),
        legend.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig2.pdf")
ggsave(file_path, width = 140, height = 90, unit = "mm")

# now look at forecasts spatially
i_time <- 1449  # 3 days later Mag 6.1 earthquake
pred_one_day <- data.frame(do.call(cbind, lapply(models, function(m) m[i_time,]))) %>%
  mutate(N = 1:nrow(.))
colnames(pred_one_day) <- c(model_names, "N")

pred_by_cell_long <- pred_one_day %>%
  right_join(cells, by = "N") %>%
  select(LON, LAT, all_of(model_names)) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model")

events_by_cell <- events %>%
  filter(TS >= times[i_time], TS < times[i_time] + days(7)) %>%
  group_by(N) %>%
  summarise(Count = n()) %>%
  left_join(cells, by = "N") %>%
  select(LON, LAT, Count)

ggplot() +
  facet_wrap(~Model, ncol = 2) +
  geom_tile(data = pred_by_cell_long,
            aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", alpha = 0.4,
          size = 0.2, fill = NA) +
  geom_tile(data = events_by_cell, aes(x = LON, LAT, color = "Observed earthquakes"), fill = NA) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_viridis_c(name = "Predicted\nmean",
                       breaks = 10^(c(-9, -7, -5, -3)), labels = paste0("e", c(-9, -7, -5, -3)),
                       trans = "log10", option = "magma") +
  scale_color_manual(name = "Observed\nearthquakes", values = c("Observed earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  ggtitle(paste("Predictions for the", times[i_time], "-", times[i_time] + days(7))) +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = 9))

file_path <- file.path(fpath, "Poster_Fig3.pdf")
ggsave(file_path, width = 200, height = 180, unit = "mm")

rm(pred_by_day, pred_by_cell_long, events_by_cell, pred_one_day)

################################################################################
# Visualize Poisson Scores over time and score differences spatially
################################################################################

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

scores <- do.call(cbind, lapply(models, function(X) rowMeans(S_pois(X, obs))))
scores <- data.frame(scores) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(scores) <- c(model_names, "X", "earthquake")

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

## render_strips?

no_eq <- diff_scores %>%
  mutate(value = ifelse(earthquake, NA, value)) %>%
  ggplot() +
  facet_wrap(~"Days with no earthquake") +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year])) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(c("", " ", ""), ana_models, "vs.", cmp_model)) +
  scale_y_continuous(labels = scientific) +
  xlab("") +
  ylab("Score") +
  ggtitle("Daily Average Poisson Score Differences") +
  theme_bw() +
  my_theme +
  theme(legend.justification = c(0, 1), legend.position = c(0.01, 0.99),
        legend.direction = "vertical", strip.background = element_blank(),
        legend.key.size = unit(0.5, "lines"))

yes_eq <- diff_scores %>%
  mutate(value = ifelse(earthquake, value, NA)) %>%
  ggplot() +
  facet_wrap(~"Days with at least one earthquake") +
  geom_point(aes(x = X, y = value, color = Model), show.legend = FALSE, size = 0.4) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year])) +
  scale_color_manual(name = NULL, values = model_colors) +
  scale_y_continuous(labels = scientific) +
  xlab("") +
  ylab(NULL) +
  ggtitle("") +
  theme_bw() +
  my_theme +
  theme(strip.background = element_blank())

stack_plots <- grid.arrange(no_eq, yes_eq, nrow = 1)
file_path <- file.path(fpath, "Poster_Fig4.pdf")
ggsave(file_path, width = 200, height = 110, unit = "mm", plot = stack_plots)

# and now spatially

scores <- do.call(cbind, lapply(models, function(X) colMeans(S_pois(X, obs))))
scores <- data.frame(scores)
colnames(scores) <- model_names

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  mutate(LON = cells$LON, LAT = cells$LAT) %>%
  select(LON, LAT, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")
diff_scores$Model <- paste(diff_scores$Model, "vs.", cmp_model)

my_colors <- c("#ff9603", "#f51818", "#ffffff", "#057ffa")
limits <- range(diff_scores$value)
my_breaks <- c(-0.01, -0.001, 0, 0.001)
my_labels <- c("-1e-2", "-1e-3", " 0", " 1e-3")

# need log transform for positive and negative values (see ?modulus_trans)
# but need to scale with d to get sufficient resolution
d <- 10^5
my_trans <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d  + 1),
  function(y) sign(y) / d * (exp(abs(y)) - 1)
)
col_breaks <- my_trans$transform(c(limits[1], -limits[2], 0, limits[2]))

ggplot() +
  facet_grid(~Model) +
  geom_tile(data = diff_scores,
            aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", fill = NA,
          size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "Score\ndifference",
                       trans = my_trans, colors = my_colors, values = rescale(col_breaks),
                       breaks = my_breaks, labels = my_labels) +
  ggtitle("Average Score Differences") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank())

file_path <- file.path(fpath, "Poster_Fig5.pdf")
ggsave(file_path, width = 180, height = 80, unit = "mm")

rm(scores, scores_acc, obs_acc, m_neigh, diff_scores, diff_scores_acc, combine)

################################################################################
# Visualize reliability diagram
################################################################################

# adapted from
# https://github.com/dwolffram/replication-ARSIA2023/blob/main/R/reliability_functions.R
reldiag <- function(x, y, n_resamples = 999, region_level = 0.9) {

  recal <- function(x, y) {
    ord <- order(x, y, decreasing = c(FALSE,TRUE))
    x_rc <- monotone(y[ord])
    return(x_rc[order(ord)])  # get original order of x vectors back
  }

  score <- function(x, y) mean(S_pois(x, y), na.rm = TRUE)
  # score <- function(x, y) mean(S_quad(x, y), na.rm = TRUE)

  x_rc <- recal(x, y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  n_samples <- n_resamples + 1 # total number of samples including observed sample
  low <- floor(n_samples * (1 - region_level) / 2)
  up <- n_samples - low

  bootstrap <- sapply(1:n_resamples, function(i) x + sample(y - x, length(y), replace = TRUE)) %>%
    as.data.table() %>% mutate_all(function(c) pmax(0, c)) %>%
    mutate_all(function(c) recal(x, c)) %>%   # calibrate boostrap samples
    mutate(old = x_rc) %>%
    apply(1, sort)       # apply will transpose data.table!

  results <- data.frame(x = x, x_rc = x_rc, lower = bootstrap[low,],
                        upper = bootstrap[up,])
  stats <- data.frame(Score = s, MCB = mcb, DSC = dsc, UNC = unc)
  return(list(results = results, stats = stats))
}

recal_models <- data.table()
collect_stats <- data.table()

for (i in 1:length(models)) {
  res <- reldiag(rowSums(models[[i]]), rowSums(obs))
  recal_models <- rbind(recal_models, cbind(Model = model_names[i], res$results))
  collect_stats <- rbind(
    collect_stats,
    cbind(Model = model_names[i],
          label = paste(names(res$stats), format(round(res$stats, 3), nsmall = 3,
                                                 scientific = FALSE),
                        collapse = "\n"))
  )
}

ggplot(recal_models, aes(x = x)) +
  facet_wrap(~Model, nrow = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.33,
              show.legend = FALSE) +
  geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
              linetype = "dashed") +
  geom_line(aes(y = x_rc, color = Model), size = 0.5, show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values = model_colors) +
  scale_x_log10() +
  scale_y_log10(limits = c(0.01, 10)) +     # for quadratic loss, just leave limits away
  xlab("Forecasted mean") +
  ylab("Conditional mean") +
  ggtitle("Reliability Diagram of Daily Forecasts") +
  geom_text(data = collect_stats, mapping = aes(x = 0.08, y = 2.7, label = label),
             size = 6*0.36, hjust = 0) +
  theme_bw() +
  my_theme +
  theme(strip.background = element_blank(), aspect.ratio = 1)

file_path <- file.path(fpath, "Poster_Fig6.pdf")
ggsave(file_path, width = 180, height = 70, unit = "mm")

rm(recal_models, collect_stats)

################################################################################
# Visualize Murphy Diagram
################################################################################

S_elem <- compiler::cmpfun(s_theta)

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

data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.5) +
  scale_color_manual(name = NULL, values = model_colors) +
  xlab(expression(paste("Threshold log", theta))) +
  ylab("Elementary score") +
  ggtitle("Murphy Diagram") +
  theme_bw() +
  my_theme +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1))

file_path <- file.path(fpath, "Poster_Fig7.pdf")
ggsave(file_path, width = 150, height = 100, unit = "mm")

rm(murphy_df)
