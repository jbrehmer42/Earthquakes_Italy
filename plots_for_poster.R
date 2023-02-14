library(ggplot2)
library(dplyr)
library(tidyr)

library(sf)               # st_as_sf : convert data.frame to geographic format sf
library(rnaturalearth)    # ne_countries - to load country data
library(scales)           # seq_gradient_pal - custom color gradients
                          # trans_new - custom data transformation

source("data_prep.R")

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
  ggtitle("M4+ Earthquakes duuring Training Period") +
  theme_bw() +
  my_theme +
  theme(panel.background = element_rect(fill = rgb(0.145, 0.588, 0.745, 0.4)))

file_path <- file.path(fpath, "Poster_Fig1.pdf")
ggsave(file_path, width = 140, height = 140, unit = "mm")

rm(europe)

################################################################################
# Visualize predictions over time and over grid cells
################################################################################

# first look at forecasts temporally
pred_by_day <- data.frame(do.call(cbind, lapply(models, rowSums))) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(pred_by_day) <- c(model_names, "X", "earthquake")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0
eq_shape <- 3

pred_by_day %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_hline(data = data.frame(Model = "Climatology", value = sum(clima$RATE)),
             aes(yintercept = value, color = Model), size = 0.3) +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  geom_point(data = filter(pred_by_day, earthquake),
             aes(x = X, y = 0.15, color = "Observed earthquakes"),
             size = 0.7, stroke = 0.4, shape = eq_shape) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(
    name = NULL,
    values = c("Climatology" = "#7d7b7a", model_colors, "Observed earthquakes" = "black"),
    guide = guide_legend(override.aes = list(linetype = c(rep(1, 5), 0),
                                             shape = c(rep(NA, 5), eq_shape)), nrow = 5)
  ) +
  xlab("") +
  ylab("Predicted mean") +
  ggtitle("Predictions for all of Italy") +
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
ggsave(file_path, , width = 200, height = 180, unit = "mm")
