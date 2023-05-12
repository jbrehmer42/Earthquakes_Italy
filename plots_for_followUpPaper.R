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
source("functions_eval.R")

model_colors <- c("FMC" = "#F8766D", "LG" = "#00BA38", "SMA" = "#619CFF",
                  "LM" = "#DB72FB")

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

fpath <- "./figures2"

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

# handle zero forecast specifically
s_pois <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(score)
}

s_pois_da <- function(X, Y) {
  zero_fcst <- X == 0
  impossible_fcst <- (X == 0) & (Y != 0)
  score <- -Y * log(X) + X
  score[zero_fcst] <- 0
  score[impossible_fcst] <- Inf
  return(sum(score) / n_days)
}

s_quad_da <- function(X, Y) {
  return(sum((X - Y)^2) / n_days)
}

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
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL) +
  scale_color_gradient(low = m_low, high = m_high) +
  scale_size(range = new_range, trans = size_trans) +
  guides(color = "none", size = guide_legend(title = "Magnitude",
                                             override.aes = list(color = mag_colors))) +
  my_theme +
  theme(panel.background = element_rect(fill = rgb(0.012, 0.663, 0.988, 0.3)),
        legend.text = element_text(size = 7))

my_breaks <- c(4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 10.0)
eq_count <- events %>%
  group_by(MAG = sapply(MAG, function(m) sum(m >= my_breaks))) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(MAG = factor(my_breaks[MAG], ordered = T, levels = rev(my_breaks))) %>%
  ggplot() +
  geom_col(aes(y = MAG, x = count, fill = MAG), alpha = 0.95, show.legend = F) +
  scale_fill_manual(values = setNames(mag_colors, my_breaks[-length(my_breaks)])) +
  scale_x_log10(name = "Count [log]") +
  scale_y_discrete(labels = NULL) +
  ylab(NULL) +
  theme_classic(base_size = 8) +
  theme(axis.title = element_text(size = 7), plot.background = element_blank())

combine <- ggplot() +
  # set axis limits
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
  annotation_custom(ggplotGrob(eq_map), xmin = 0, xmax = 0.8, ymin = 0, ymax = 1) +
  annotation_custom(ggplotGrob(eq_count), xmin = 0.75, xmax = 1, ymin = 0.21, ymax = 0.73) +
  theme(panel.background = element_rect("white"))

finish <- grid.arrange(combine,
                        top = textGrob("Earthquakes in Italy",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig1_Earthquakes.pdf")
ggsave(file_path, width = 140, height = 100, unit = "mm", plot = finish)

rm(eq_map, eq_count, combine)

################################################################################
# Figure 2: Distribution of forecasts over time and space
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
             aes(xintercept = X, linetype = "Obs. earthquakes"), size = 0.3,
             color = "gray", alpha = 0.2) +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  scale_linetype_manual(name = NULL, values = c("Obs. earthquakes" = 1),
                        guide = guide_legend(override.aes = list(alpha = 1, size = 0.5))) +
  xlab(NULL) +
  ylab("Predicted mean") +
  labs(subtitle = "For all of Italy") +
  annotate(geom = "text", label = "(*)", x = i_time, y = 0.15, size = 3,
           color = "black") +
  theme_bw() +
  my_theme +
  theme(legend.position = "bottom", legend.margin = margin(0, 5.5, 5.5, 5.5))

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
  geom_tile(data = events_by_cell, aes(x = LON, y = LAT, color = "Obs. earthquakes"), fill = NA) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 12, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_viridis_c(name = "Predicted mean",
                       breaks = 10^(c(-9, -7, -5, -3)),
                       labels = paste0("1e", c(-9, -7, -5, -3)),
                       trans = "log10", option = "magma",
                       guide = guide_colorbar(title.vjust = 0.5, order = 1)) +
  scale_color_manual(name = "Obs. earthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"),
                                          title.vjust = 0.6, order = 2)) +
  labs(subtitle = "For the 7-day Period Following (*)") +
  theme_bw() +
  my_theme +
  theme(legend.position = "bottom", legend.box.just = "left",
        plot.margin = margin(2.5, 5.5, 5.5, 10.5),
        legend.margin = margin(0, 5.5, 5.5, 5.5))

combine <- grid.arrange(temp_plot, spat_plot, nrow = 2, heights = c(0.45, 0.55),
                        top = textGrob("Predicted Mean Number of Events",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig2_Forecasts.pdf")
ggsave(file_path, width = 145, height = 138, unit = "mm", plot = combine)

rm(pred_by_day, pred_by_cell_long, events_by_cell, pred_one_day, temp_plot,
   spat_plot, combine)

################################################################################
# Tables: Calculate Scores
################################################################################

# Table 1: Overall quadratic and Poisson score and its number and spatial component
t1 <- matrix(NA, nrow = length(models), ncol = 4)
rownames(t1) <- model_names
colnames(t1) <- c("quad", "Poisson", "number", "spatial")

for (i in 1:length(models)) {
  x_t <- rowSums(models[[i]])
  t1[i, "quad"] <- s_quad_da(models[[i]], obs)
  t1[i, "Poisson"] <- s_pois_da(models[[i]], obs)
  t1[i, "number"] <- mean(s_pois(x_t, rowSums(obs)))
  t1[i, "spatial"] <- mean(rowSums(s_pois(models[[i]] / x_t, obs)))
}
t1
write.csv(t1, "./../tmp_results/Table1.csv")

# Table 2: Overall quadratic and Poisson score and its MSB, DSC, and UNC component
t2 <- matrix(NA, nrow = length(models), ncol = 8)
rownames(t2) <- model_names
colnames(t2) <- c("quad", "q-MCB", "q-DSC", "q-UNC", "pois", "p-MCB", "p-DSC",
                  "p-UNC")

# we get different scores due to sorting and summing up in a different order?

y <- as.vector(obs)
mean_y <- mean(y)
for (i in 1:length(models)) {
  # recalibrate forecasts x with isotoinc regression from monotone package
  x <- as.vector(models[[i]])
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  j <- 1
  for (scf in list(s_quad_da, s_pois_da)) {
    s <- scf(x, y)
    s_rc <- scf(x_rc, y)
    s_mg <- scf(mean_y, y)
    t2[i, j] <- s
    t2[i, j + 1] <- s - s_rc
    t2[i, j + 2] <- s_mg - s_rc
    t2[i, j + 3] <- s_mg
    j <- j + 4
  }
}
t2
write.csv(t1, "./../tmp_results/Table2.csv")

################################################################################
# Figure 5: Poisson Score over days for each grid cell
################################################################################

scores <- do.call(cbind, lapply(models, function(X) colMeans(s_pois(X, obs))))
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

scores_acc <- do.call(cbind, lapply(models, function(X) colMeans(s_pois(X %*% m_neigh, obs_acc))))
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
  theme(legend.position = "right", axis.text.x = element_blank())

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
  theme(legend.position = "right", strip.text.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 12))

combine <- grid.arrange(single_cells, acc_cells, nrow = 2,
                        top = textGrob("Mean Poisson Score Difference",
                                       gp = gpar(fontsize=15)))

file_path <- file.path(fpath, "Fig5_spatialScoreDifferences.pdf")
ggsave(file_path, width = 190, height = 160, unit = "mm", plot = combine)

rm(diff_scores, m_neigh, obs_acc, scores_acc, diff_scores_acc)

################################################################################
# Visualize Poisson Scores over time and score differences spatially
################################################################################

scores <- do.call(cbind, lapply(models, function(X) colMeans(s_pois(X, obs))))
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
my_labels <- c("-1e-2", "", "-1e-4", "0", "1e-4", "")

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
                       breaks = my_breaks, labels = my_labels, minor_breaks = my_mbreaks,
                       guide = guide_colorbar(barwidth = unit(40, "mm"),
                                              title.vjust = 0.9)) +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  theme_bw() +
  my_theme +
  theme(legend.position = "bottom")

combine_plots <- grid.arrange(spat_plot, nrow = 1,
                              top = textGrob("Average Poisson Score Differences by Grid Cell",
                                             gp = gpar(fontsize = title_size)))
file_path <- file.path(fpath, "Fig6_ScoreDiffSpat.pdf")
ggsave(file_path, width = 145, height = 90, unit = "mm", plot = combine_plots)

# and now temporally

scores <- do.call(cbind, lapply(models, function(X) rowSums(s_pois(X, obs))))
scores <- data.frame(scores) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(scores) <- c(model_names, "X", "earthquake")

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

d2 <- 10^3
my_trans2 <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d2  + 1),
  function(y) sign(y) / d2 * (exp(abs(y)) - 1)
)
my_breaks <- c(-10^(c(2, 0, -2)), 0, 10^(c(-2, 0)))
minor_breaks <- c(-10^(2:-3), 0, 10^(-3:1))
my_labels <- c("-1e2", "-1", "-1e-2", "0", "1e-2", "1")

temp_plot <- ggplot(diff_scores) +
  geom_vline(data = filter(diff_scores, earthquake > 0),
             aes(xintercept = X, linetype = "Obs. earthquakes"),
             alpha = 0.3, color = "gray", size = 0.3) +
  geom_point(aes(x = X, y = value, color = Model), size = 0.3, alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores))) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models),
                     guide = guide_legend(order = 1, direction = "horizontal",
                                          override.aes = list(alpha = 1))) +
  scale_linetype_manual(name = NULL, values = c("Obs. earthquakes" = 1),
                        labels = "Obs. earthquakes",
                        guide = guide_legend(override.aes = list(alpha = 1, size = 0.4),
                                             order = 2, direction = "horizontal")) +
  scale_y_continuous(trans = my_trans2, breaks = my_breaks, labels = my_labels,
                     minor_breaks = minor_breaks) +
  xlab(NULL) +
  ylab("Score difference") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "bottom", legend.box = "horizontal",)

combine_plots <- grid.arrange(temp_plot, nrow = 1,
                              top = textGrob("Poisson Score Differences By Day",
                                             gp = gpar(fontsize = title_size)))
file_path <- file.path(fpath, "Fig4_ScoreDiffTemp.pdf")
ggsave(file_path, width = 145, height = 80, unit = "mm", plot = combine_plots)

rm(scores, diff_scores, combine_plots, spat_plot, temp_plot)

################################################################################
# Poisson Score: Number and spatial component
################################################################################

# score of predicted number of earthquakes in all of Italy
scores_n <- do.call(
  cbind, lapply(models, function(X) s_pois(rowSums(X), rowSums(obs)))
)
scores_n <- data.frame(scores_n)
colnames(scores_n) <- model_names

diff_scores <- scores_n %>%
  mutate(X = 1:nrow(.)) %>%
  filter(rowSums(obs) > 0) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

number_plot <- ggplot(diff_scores) +
  geom_point(aes(x = X, y = value, color = Model), size = 0.75) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = which(new_year), labels = year(times[new_year]),
                     limits = c(0, nrow(scores_n))) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models)) +
  scale_y_continuous(limits = c(-20, NA)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(subtitle = "Number component") +
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

# score of normalized mean predicted number of earthquakes
scores_s <- do.call(
  cbind, lapply(models, function(X) {
    row_sums <- rowSums(X)
    X_norm <- apply(X, 2, function(col) col / row_sums)
    return(rowSums(s_pois(X_norm, obs)))
  })
)
scores_s <- data.frame(scores_s)
colnames(scores_s) <- model_names

diff_scores <- scores_s %>%
  mutate(X = 1:nrow(.)) %>%
  filter(rowSums(obs) > 0) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

spatial_plot <- ggplot(diff_scores) +
  geom_point(aes(x = X, y = value, color = Model), size = 0.75, show.legend = F) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = which(new_year), labels = year(times[new_year]),
                     limits = c(0, nrow(scores_s))) +
  scale_color_manual(name = NULL, values = model_colors, breaks = ana_models,
                     labels = paste(cmp_model, "vs.", ana_models)) +
  scale_y_continuous(limits = c(-20, NA)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(subtitle = "Spatial component") +
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5))

combine_plots <- grid.arrange(number_plot, spatial_plot, nrow = 2,
                              top = textGrob("Daily Poisson Score Differences",
                                             gp = gpar(fontsize = title_size)),
                              left = textGrob("Score difference", rot = 90,
                                              gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "Fig5_NumberSpatials.pdf")
ggsave(file_path, width = 145, height = 135, unit = "mm", plot = combine_plots)

rm(scores_n, diff_scores, number_plot, scores_s, spatial_plot, combine_plots)

################################################################################
# Murphy diagramm
################################################################################

S_elem <- compiler::cmpfun(S_theta)

n_theta <- 100
log_grid <- seq(-24, 4, len = n_theta)
grd <- exp(log_grid)

# use for loops, as otherwise memory demand is too high
# Jonas calculates each row separately, it seems faster? check it!
murphy_df <- matrix(0, nrow = length(grd), ncol = length(models))
for (m in 1:length(models)) {
  print(m)
  for (t in 1:n_theta) {
    cat("*")
    murphy_df[t, m] <- S_elem(models[[m]], obs, grd[t])
  }
}
colnames(murphy_df) <- model_names

# S_theta sums, but we want daily averages
murphy_df <- murphy_df / nrow(obs)

# or read it in
murphy_df <- read.csv("./../tmp_results/murphy_df.csv")

murphy_diag <- data.frame(murphy_df) %>%
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
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0),
        legend.background = element_blank())

combine <- grid.arrange(murphy_diag, nrow = 1,
                        top = textGrob("Murphy Difference Diagram",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig7_MurphyDiag.pdf")
ggsave(file_path, width = 145, height = 80, unit = "mm", plot = combine)

# now look at Murphy diagram of miscalibration and discrimination component
MCB_diag <- DSC_diag <- matrix(0, ncol = n_mods, nrow = n_theta)
for (i in 1:n_mods) {
  print(i)
  decomp <- cell_decomposition(models[[i]], obs, theta = grd)
  # decomp <- day_decomposition(models[[i]], obs, theta = grd) # in helpful_routines.R
  MCB_diag[ ,i] <- rowSums(decomp$MCB)
  DSC_diag[ ,i] <- rowSums(decomp$DSC)
}
UNC_diag <- rowSums(decomp$UNC)

# alternatively: recalibrate everything at once
UNC_diag <- numeric(length(grd))

for (i in 1:n_mods) {
  print(i)
  # recalibrate forecasts x with isotoinc regression from monotone package
  x <- as.vector(models[[i]])
  y <- as.vector(obs)
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  for (t in 1:length(grd)) {
    cat("*")
    s <- S_elem(x, y, grd[t])
    s_rc <- S_elem(x_rc, y, grd[t])
    s_mg <- S_elem(mean(y), y, grd[t])

    MCB_diag[t, i] <- s - s_rc
    DSC_diag[t, i] <- s_mg - s_rc
    if (i == 1) {
      UNC_diag[t] <- s_mg
    }
  }
}

df_collect <- rbind(
  cbind(data.frame(MCB_diag), Type = "MCB"), cbind(data.frame(DSC_diag), Type = "DSC")
)
colnames(df_collect) <- c(model_names, "Type")

df_collect <- df_collect %>%
  mutate(log_theta = rep(log_grid, 2)) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  mutate(Type = ifelse(Type == "MCB", "Miscalibration", "Discrimination"))

# or read precomputed values in
df_collect <- read.csv("./../tmp_results/murphy-MCB-DSC_allRecal.csv")

murphy_score_cmps <- ggplot(df_collect) +
  facet_wrap(~factor(Type, ordered = T, levels = c("Miscalibration", "Discrimination")),
             nrow = 1, scales = "free_y") +
  geom_line(aes(x = log_theta, y = value, color = Model), size = 0.3) +
  scale_color_manual(name = NULL, values = model_colors) +
  scale_x_continuous(breaks = -4:1 * 5) +
  xlab(expression(paste("Threshold log", theta))) +
  ylab(NULL) +
  my_theme +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_blank())

combine <- grid.arrange(murphy_score_cmps, nrow = 1,
                        top = textGrob("Score Components by Elementary Score",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig9_Murphy-MCB-DSC.pdf")
ggsave(file_path, width = 145, height = 80, unit = "mm", plot = combine)

rm(murphy_df, decomp, murphy_diag, df_collect, murphy_score_cmps, combine,
   MCB_diag, DSC_diag, UNC_diag)

################################################################################
# Reliability diagramm
################################################################################

reldiag <- function(x, y, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # we can compress step function by only storing first and last value and jumps!
  filter_jumps <- function(v) {
    n <- length(v)
    return(c(T, v[c(-1, -n)] - v[c(-n+1, -n)] > 0, T))
  }

  score <- s_pois_da

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

# or load already recalibrated values
recal_models <- read.csv("./../tmp_results/recal_models_100_bag-Tilmann.csv")
collect_stats <- read.csv("./../tmp_results/collect_stats.csv")

col_ecdfs <- list()
for (i in 1:length(models)) {
  t <- table(as.vector(models[[i]]))
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / prod(dim(models[[i]])),
                               M = model_names[i])
}

collect_recal_plots <- list()

# for non-aggregated data
breaks_x <- list(c(0, 10^c(-5, 0)), c(0, 10^c(-6, -5, 0)),
                 c(0, 10^c(-6, -5, 0)), c(0, 10^c(-5, 0)))
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

  my_hist <- ggplot(data.table(x = as.vector(models[[i]]))) +
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
    my_theme +
    theme(aspect.ratio = 1)

  collect_recal_plots[[i]] <- main_plot + inset_histograms
}

combine <- grid.arrange(collect_recal_plots[[2]],
                        collect_recal_plots[[3]],
                        collect_recal_plots[[1]],
                        collect_recal_plots[[4]],
                        nrow = 2,
                        top = textGrob("Reliability Diagram",
                                       gp = gpar(fontsize = title_size)),
                        bottom = textGrob("Forecasted mean",
                                       gp = gpar(fontsize = 11)),
                        left = textGrob("Conditional mean", rot = 90,
                                       gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "Fig7_ReliabilityDiagram.pdf")
ggsave(file_path, width = 145, height = 160, unit = "mm", plot = combine)

rm(combine, collect_recal_plots, col_ecdfs, collect_stats, recal_models)
