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

rm(eq_map, eq_hist, combine)

################################################################################
# Figure 1: Distribution of forecasts
################################################################################

# first look at forecasts temporally
pred_by_day <- data.frame(do.call(cbind, lapply(models, rowSums))) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(pred_by_day) <- c(model_names, "X", "earthquake")

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0
eq_shape <- 3   # shape of marks for earthquakes

temp_plot <- pred_by_day %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_vline(data = filter(pred_by_day, earthquake),
             aes(xintercept = X, linetype = "Obs.\nearthquakes"), size = 0.3,
             color = "gray", alpha = 0.4) +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(name = NULL, values = model_colors) +
  scale_linetype_manual(name = NULL, values = c("Obs.\nearthquakes" = 1),
                        guide = guide_legend(override.aes = list(alpha = 1))) +
  xlab(NULL) +
  ylab("Predicted mean") +
  ggtitle("For all of Italy") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", legend.box = "vertical",
        legend.background = element_blank(), plot.title = element_text(size = 12))

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
  ggtitle("For a single day") +
  theme_bw() +
  my_theme +
  theme(legend.position = "right", strip.background = element_blank(),
        legend.title = element_text(size = 8), legend.box.just = "left",
        plot.title = element_text(size = 12))

combine <- grid.arrange(temp_plot, spat_plot, nrow = 2, heights = c(0.45, 0.55),
                        top = textGrob("Predicted Mean Number of Events",
                                 gp = gpar(fontsize=15)))

file_path <- file.path(fpath, "Fig1b_Forecasts.pdf")
ggsave(file_path, width = 250, height = 170, unit = "mm", plot = combine)

rm(pred_by_day, pred_by_cell_long, events_by_cell, pred_one_day)

################################################################################
# Figure 2: Poisson Score over days
################################################################################