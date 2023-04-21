library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

library(gridExtra)        # grid.arrange : put ggplots next to each other
library(grid)             # for grobText : change text size in grid.arrange
library(monotone)         # for fast isotonic regression
library(sf)               # st_as_sf : convert data.frame to geographic format sf
library(rnaturalearth)    # ne_countries - to load country data
library(scales)           # trans_new - custom data transformation
                          # seq_gradient_pal - custom color gradient

source("data_prep.R")

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")

fpath <- "./figures2"

my_theme <- list(
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05))
)

################################################################################
# Figure 1: Distribution of earthquakes
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
  "linearScaleDown",
  function(x) x - 2, function(x) x + 2
)
new_range <- range(size_trans$transform(events$MAG))

eq_map <- ggplot() +
  # plot has a blue background: first draw countries white again
  geom_sf(data = europe, fill = "white", size = 0.2) +
  # before painting them in light green
  geom_sf(data = europe, fill = "#04b50a", color = "black", alpha = 0.3, size = 0.2) +
  geom_tile(data = cells, aes(x = LON, y = LAT), color = "#f5a536", fill = NA) +
  geom_sf(data = st_as_sf(events, coords = c("LON", "LAT"), crs = 4326),
          aes(size = MAG, color = MAG), alpha = 0.6, shape = 13, stroke = 0.3) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL) +
  scale_y_continuous(name = NULL) +
  scale_color_gradient(low = m_low, high = m_high) +
  scale_size(range = new_range, trans = size_trans) +
  guides(color = "none", size = guide_legend(title = "Magnitude",
                                             override.aes = list(color = mag_colors))) +
  ggtitle("M4+ Earthquakes") +
  my_theme +
  theme(panel.background = element_rect(fill = rgb(0.012, 0.663, 0.988, 0.3)))

my_breaks <- c(4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0)
eq_count <- events %>%
  group_by(MAG = sapply(MAG, function(m) sum(m >= my_breaks))) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(MAG = factor(my_breaks[MAG], ordered = T, levels = rev(my_breaks))) %>%
  ggplot() +
  geom_col(aes(y = MAG, x = count, fill = MAG), show.legend = F) +
  scale_fill_manual(values = setNames(mag_colors, my_breaks[-length(my_breaks)])) +
  scale_x_log10(name = "Count [log]") +
  scale_y_discrete(labels = NULL) +
  ylab(NULL) +
  theme_classic(base_size = 6)

combine <- ggplot() +
  # set axis limits
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
  annotation_custom(ggplotGrob(eq_map), xmin = 0, xmax = 0.8, ymin = 0, ymax = 1) +
  annotation_custom(ggplotGrob(eq_count), xmin = 0.75, xmax = 1, ymin = 0.2, ymax = 0.67) +
  theme(panel.background = element_rect("white"))

file_path <- file.path(fpath, "Fig1_Earthquakes.pdf")
ggsave(file_path, width = 140, height = 100, unit = "mm", plot = combine)

rm(eq_map, eq_count, combine)

################################################################################
# Figure 1b: Distribution of forecasts over time and space
################################################################################

# first look at forecasts temporally
pred_by_day <- data.frame(do.call(cbind, lapply(models, rowSums))) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(pred_by_day) <- c(model_names, "X", "earthquake")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

i_time <- 1449  # 3 days later Mag 6.1 earthquake

temp_plot <- pred_by_day %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_vline(data = filter(pred_by_day, earthquake),
             aes(xintercept = X, linetype = "Obs.\nearthquakes"), size = 0.3,
             color = "gray", alpha = 0.3) +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  scale_linetype_manual(name = NULL, values = c("Obs.\nearthquakes" = 1),
                        guide = guide_legend(override.aes = list(alpha = 1, size = 0.5))) +
  xlab(NULL) +
  ylab("Predicted mean") +
  ggtitle("For all of Italy") +
  annotate(geom = "text", label = "(*)", x = i_time, y = 0.15, size = 3,
           color = "black") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", legend.box = "vertical",
        legend.text = element_text(size = 8), legend.background = element_blank(),
        plot.title = element_text(size = 12), plot.margin = margin(5.5, 5.5, 0.0, 3.0))

# now look at forecasts spatially
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

spat_plot <- ggplot() +
  facet_wrap(~Model, nrow = 1) +
  geom_tile(data = pred_by_cell_long, aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", alpha = 0.4,
          size = 0.2, fill = NA) +
  geom_tile(data = events_by_cell, aes(x = LON, y = LAT, color = "Obs.\nearthquakes"), fill = NA) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_viridis_c(name = "Predicted mean",
                       breaks = 10^(c(-9, -7, -5, -3)), labels = paste0("e", c(-9, -7, -5, -3)),
                       trans = "log10", option = "magma",
                       guide = guide_colorbar(title.vjust = 0.5, order = 1)) +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs.\nearthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"),
                                          title.vjust = 0.6, order = 2)) +
  ggtitle("For a single day (*)") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = 8), legend.box.just = "left",
        plot.title = element_text(size = 12), plot.margin = margin(5.5, 15.5, 5.5, 9))

combine <- grid.arrange(temp_plot, spat_plot, nrow = 2, heights = c(0.45, 0.55),
                        top = textGrob("Predicted Mean Number of Events",
                                       gp = gpar(fontsize=15)))

file_path <- file.path(fpath, "Fig1b_Forecasts.pdf")
ggsave(file_path, width = 220, height = 145, unit = "mm", plot = combine)

rm(pred_by_day, pred_by_cell_long, events_by_cell, pred_one_day, temp_plot,
   spat_plot, combine)

################################################################################
# Figure 5: Poisson Score over days for each grid cell
################################################################################

# handle zero forecast specifically
S_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(score)
}

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
diff_scores$Model <- paste(diff_scores$Model, "vs.", cmp_model)

neigh_mat <- function(cells, k, diff = function(x, y) abs(x - y), agg = pmax) {
  # Compute neighborhood matrix for spatial aggregation:
  # binary square matrix where a 1 at position (i,j) indicates that cells i and
  # j are in each others neighborhoods.

  n_cells <- dim(cells)[1]
  # Aggregation will usually be done for small values
  # of k so "sparse = T" makes sense in most cases
  mat <- Matrix(0, nrow = n_cells, ncol = n_cells, sparse = T)
  for (i in 1:n_cells) {
    neighbors <- agg(diff(cells$X[i], cells$X), diff(cells$Y[i], cells$Y)) <= k
    mat[, i] <- as.numeric(neighbors)
  }
  return(mat)
}

m_neigh <- neigh_mat(cells, 2, function(x, y) abs(x - y))
# m_neigh <- neigh_mat(cells, 2, function(x, y) (x - y)^2, function(x, y) sqrt(x + y))
obs_acc <- obs %*% m_neigh

scores_acc <- do.call(cbind, lapply(models, function(X) colMeans(S_pois(X %*% m_neigh, obs_acc))))
scores_acc <- data.frame(scores_acc)
colnames(scores_acc) <- model_names

diff_scores_acc <- scores_acc %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  mutate(LON = cells$LON, LAT = cells$LAT) %>%
  select(LON, LAT, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")
diff_scores_acc$Model <- paste(diff_scores_acc$Model, "vs.", cmp_model)

my_colors <- c("#de7702", muted("red"), "#ffffff", muted("blue"))
limits <- range(diff_scores_acc$value)
col_breaks_1 <- c(limits[1], -limits[2], 0, limits[2])

single_cells <- ggplot(cbind(diff_scores, R = "")) +
  facet_grid(R ~ Model) +
  geom_tile(aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", fill = NA,
          size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "Score\ndifference", colors = my_colors,
                       values = rescale(col_breaks_1)) +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = 8), axis.text.x = element_blank())

limits <- range(diff_scores_acc$value)
col_breaks_2 <- c(limits[1], -limits[2], 0, limits[2])

acc_cells <- ggplot(cbind(diff_scores_acc, R = "Agg. cells with a radius of 2")) +
  facet_grid(R ~ Model) +
  geom_tile(aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", fill = NA,
          size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "", colors = my_colors,
                       values = rescale(col_breaks_1)) +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        strip.text.x = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 12))

combine <- grid.arrange(single_cells, acc_cells, nrow = 2,
                        top = textGrob("Mean Poisson Score Difference",
                                       gp = gpar(fontsize=15)))

file_path <- file.path(fpath, "Fig5_spatialScoreDifferences.pdf")
ggsave(file_path, width = 220, height = 160, unit = "mm", plot = combine)

rm(diff_scores, m_neigh, obs_acc, scores_acc, diff_scores_acc)

################################################################################
# Visualize Poisson Scores over time and score differences spatially
################################################################################

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
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "Score\ndifference",
                       trans = my_trans, colors = my_colors, values = rescale(col_breaks),
                       breaks = my_breaks, labels = my_labels) +
  ggtitle("Average Poisson Score Differences") +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = 8),
        plot.margin = margin(5.5, 5.5, 5.5, 27))

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
             alpha = 0.3, color = "gray", size = 0.3) +
  geom_line(aes(x = X, y = value, color = Model), size = 0.4) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores))) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models),
                     guide = guide_legend(order = 1, label.position = "top")) +
  scale_linetype_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = 1),
                        labels = "",
                        guide = guide_legend(override.aes = list(alpha = 1, size = 0.4),
                                             order = 2)) +
  scale_y_continuous(trans = my_trans2, breaks = my_breaks, labels = my_labels,
                     minor_breaks = minor_breaks) +
  xlab(NULL) +
  ylab("log-transformed score") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "right", legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

combine_plots <- grid.arrange(spat_plot, temp_plot, nrow = 2, heights = c(0.65, 0.35))

file_path <- file.path(fpath, "Fig5_ScoreDifferences.pdf")
ggsave(file_path, width = 220, height = 155, unit = "mm", plot = combine_plots)

rm(scores, diff_scores, combine_plots, spat_plot, temp_plot)
