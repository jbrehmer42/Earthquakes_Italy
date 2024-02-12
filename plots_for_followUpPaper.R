library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

library(gridExtra)        # grid.arrange : put ggplots next to each other
library(grid)             # for grobText : change text size in grid.arrange
library(cowplot)          # for get_legend : extract legend from plot
library(monotone)         # for monotone : fast isotonic mean regression
library(sf)               # for st_as_sf : convert data.frame to geographic sf format
library(rnaturalearth)    # for ne_countries : load country data
library(scales)           # for trans_new : custom data transformation
                          # for seq_gradient_pal : custom color gradient

source("data_prep.R")

# use standard colors of ggplot for discrete variables
model_colors <- c("FCM" = "#F8766D", "LG" = "#00BF7D",  "LM" = "#A3A500",
                  "SMA" = "#00B0F6", "LRWA" = "#E76BF3")

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

fpath <- "./figures8"

title_size <- 13.2  # base size 11 * 1.2 (default for theme_bw())

my_theme <- list(
  theme_bw() +
  theme(panel.grid.major = element_line(size = 0.05),
        panel.grid.minor = element_line(size = 0.05),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(5, "mm"),
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
  return(sum(score) / n_days)       # daily averages
}

s_quad <- function(X, Y) {
  return((X - Y)^2)
}

s_quad_da <- function(X, Y) {
  return(sum((X - Y)^2) / n_days)   # daily averages
}

s_theta_da <- function(x, y, theta) {
  # Elementary scoring function for the mean following Ehm et al. (2016)
  if (sum(y) == 0) {
    s <- theta * sum(x > theta)
  } else {
    s <- sum(pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x))
  }
  return(s / n_days) # daily averages
}

################################################################################
# Figure 1: Distribution of earthquakes
################################################################################

lon_lim <- range(cells$LON)
lat_lim <- range(cells$LAT)

# RAM intensive: PyCharm takes up to 3 GB
europe <- ne_countries(scale = "medium", returnclass = "sf")

m_high <- "#500000"   # old dark red "#5e0306"
m_low <- "#fc474d"
eps <- 0.001
my_breaks <- c(4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0) - eps
n_breaks <- length(my_breaks) - 1
mag_colors <- seq_gradient_pal(m_low, m_high, "Lab")(seq(0,1,length.out=n_breaks))

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
  guides(color = "none", size = "none") +
  my_theme +
  theme(panel.background = element_rect(fill = rgb(0.012, 0.663, 0.988, 0.3)),
        legend.text = element_text(size = 7))

df_hist <- events %>%
  mutate(Bin = paste(sapply(MAG, function(m) sum(m >= my_breaks))))
add_counts <- data.frame(
  x = (my_breaks[-1] + my_breaks[-length(my_breaks)]) / 2,
  counts = as.numeric(table(df_hist$Bin))
) %>%
  mutate(y = counts * 1.3, label = paste(counts))

eq_hist <- ggplot(df_hist) +
  geom_histogram(aes(x = MAG, fill = Bin), breaks = my_breaks,
                 show.legend = F) +
  # last bin has 1 observations that is not displayed due to log10 transformation
  # --> display it manually
  geom_segment(x = my_breaks[n_breaks], xend = my_breaks[n_breaks + 1],
             y = 0.01, yend = 0.01, color = mag_colors[n_breaks], size = 0.3) +
  geom_text(data = add_counts, mapping = aes(x = x, y = y, label = label), size = 6 / .pt,
            color = "black") +
  scale_fill_manual(values = setNames(mag_colors, 1:n_breaks)) +
  scale_y_log10(name = NULL) +
  scale_x_continuous(name = NULL, breaks = my_breaks + eps) +
  theme_classic(base_size = 8) +
  ggtitle("Magnitude") +
  theme(axis.title = element_text(size = 7), plot.background = element_blank(),
        plot.title = element_text(size = 8))

x_start <- 0.7
x_end <- 1
y_start <- 0.3
y_end <- 0.72
annotate_symbols <- data.frame(x = seq(x_start + 0.038, x_end - 0.001,
                                       length.out = n_breaks + 2)[c(-1, -n_breaks - 2)],
                               y = y_start - 0.02, g = factor(1:n_breaks),
                               size = (my_breaks[-1] + my_breaks[-length(my_breaks)]) / 2)

# use scale_color, scale_size to get correct format of annotated points
combine <- ggplot() +
  # set axis limits
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
  geom_point(data = annotate_symbols, aes(x = x, y = y, color = g, size = size),
             shape = 13, stroke = 0.3, show.legend = F) +
  scale_color_manual(values = mag_colors) +
  scale_size(range = new_range, trans = size_trans) +
  annotation_custom(ggplotGrob(eq_map), xmin = 0, xmax = 0.8, ymin = 0, ymax = 1) +
  annotation_custom(ggplotGrob(eq_hist), xmin = x_start, xmax = x_end, ymin = y_start,
                    ymax = y_end) +
  annotate("text", label = "Count", size = 2.5, x = 0.97, y = 0.5, angle = 90) +
  theme(panel.background = element_rect("white"), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

# Warning: Removed 30 rows containing missing values is due to histogram log10
# transformation, as we also fill each bin with a unique color
finish <- grid.arrange(combine,
                       top = textGrob("Earthquakes in Italy",
                                      gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig1_Earthquakes.pdf")
ggsave(file_path, width = 140, height = 100, unit = "mm", plot = finish)

rm(eq_map, eq_hist, combine, finish, annotate_symbols, df_hist)

################################################################################
# Figure 2: Distribution of forecasts over time and space
################################################################################

# first look at forecasts temporally
pred_by_day <- data.frame(do.call(cbind, lapply(models, rowSums))) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(pred_by_day) <- c(model_names, "X", "earthquake")

i_time <- 1449  # 3 days later Mag 6.1 earthquake

temp_plot <- pred_by_day %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_point(data = filter(pred_by_day, earthquake),
             aes(x = X, y = 0.15, shape = "Obs. earthquakes"), size = 0.75,
             color = "gray") +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5), order = 1)) +
  scale_shape_manual(name = NULL, values = c("Obs. earthquakes" = 1)) +
  xlab(NULL) +
  ylab("Expected number") +
  labs(subtitle = "For the entire CSEP Italy region") +
  annotate(geom = "text", label = "*", x = i_time, y = 0.14, size = 5,
           color = "black") +
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

# should we censor data?
pred_by_cell_long <- mutate(pred_by_cell_long, value = pmax(value, 10^-7))

events_filtered <- events %>%
  filter(TS >= times[i_time], TS < times[i_time] + days(7))

# to center align plots in facet_wrap: make plot for each row and combine them
# with grid.arrange
cb_limits <- range(pred_by_cell_long$value)   # colorbar limits

create_row <- function(row, only_legend = F) {
  pred_dt <- filter(pred_by_cell_long, Model %in% row)
  pl <- ggplot() +
    facet_wrap(~factor(Model, ordered = T, levels = row), nrow = 1) +
    geom_tile(data = pred_dt, aes(x = LON, y = LAT, fill = value), alpha = 0.5) +
    geom_sf(data = filter(europe, name == "Italy"), color = "black", alpha = 0.4,
            size = 0.2, fill = NA) +
    geom_sf(data = st_as_sf(events_filtered, coords = c("LON", "LAT"), crs = 4326),
            aes(size = MAG, color = "Obs. earthquakes"), alpha = 0.6, shape = 1,
                stroke = 0.3) +
    coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
    scale_x_continuous(name = NULL, breaks = c(6, 12, 18)) +
    scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
    scale_size(range = new_range, trans = size_trans) +
    scale_fill_viridis_c(name = "Expected\nnumber",
                        breaks = 10^(c(-9, -7, -5, -3)), limits = cb_limits,
                        labels = expression(10^-9, 10^-7, 10^-5, 10^-3),
                        trans = "log10", option = "magma", direction = -1,
                        guide = guide_colorbar(order = 1)) +
    scale_color_manual(name = "Obs.\nearth-\nquakes", values = c("Obs. earthquakes" = "black"),
                       labels = "", guide = guide_legend(order = 2, override.aes = list(size = 4))) +
    my_theme +
    theme(legend.position = "right", plot.margin = margin(5.5, 5.5, 5.5, 11.5))
  if (only_legend) {
    return(get_legend(pl + guides(size = "none")))        # only show legend
  } else {
    return(pl + guides(fill = "none", color = "none", size = "none"))    # don't show legend
  }
}
print(paste("Check:", times[i_time], "= 03.04.2009?"))
# use short dash ("en dash") = \u2013
spat_subtitle <- paste("For the 7\u00adday Period Starting April 3, 2009, marked *")
spat_plot <- grid.arrange(create_row(c("FCM", "LG", "LM")), create_row(c("SMA", "LRWA")),
                          nrow = 2)
spat_plot <- grid.arrange(spat_plot, create_row(model_names, only_legend = T), nrow = 1,
                          widths = c(0.85, 0.15),
                          top = textGrob(spat_subtitle, x = unit(0.09, "npc"),
                                         hjust = 0, gp = gpar(fontsize = 11)))

combine <- grid.arrange(temp_plot, spat_plot, nrow = 2, heights = c(0.37, 0.63),
                        top = textGrob("Expected Number of Events",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig2_Forecasts.pdf")
ggsave(file_path, width = 145, height = 190, unit = "mm", plot = combine)

rm(pred_by_day, pred_by_cell_long, events_filtered, pred_one_day, temp_plot,
   spat_plot, combine, create_row)

################################################################################
# Tables: Calculate Scores
################################################################################

print_tex_table <- function(
  table_matrix,                               # matrix containing table entries
  col_names = list(),                         # column names (mulitple lines possible)
  digits = rep(2, ncol(table_matrix)),        # number of digits in each column
  make_strong = rep("-", ncol(table_matrix)), # rows to emphasize in each columns
  strong_symbol = "\\bf",                     # which symbol to use for emphasis
  sep_after = character(0),             # add separation after specfic rows
  col_align = "c",                            # alginment of text within each column
  file_path = "table.tex"                     # location to save file
) {
  begin <- paste0("\\begin{tabular}{l ",
                  paste(rep(col_align, ncol(table_matrix)), collapse = ""), "}")
  head <- character(0)
  for (add_to_head in col_names) {
    head <- paste(head, paste0(add_to_head, collapse = " & "), "\\\\")
  }
  # Write teX code to file
  write(begin, file_path)
  write("\\toprule", file_path, append = T)
  write(head, file_path, append = T)
  write("\\midrule", file_path, append = T)
  # Write score values
  for (m in rownames(table_matrix)) {
    row_m <- paste(m)
    for (c in 1:ncol(table_matrix)) {
      fmt <- paste0("%.", digits[c], "f")
      add_bold <- ifelse(make_strong[c] == m, paste0(strong_symbol, "{"), "")
      end_bold <- ifelse(make_strong[c] == m, "}", "")
      next_entry <- sprintf(fmt, table_matrix[m, c])
      row_m <- paste(row_m, " & ", add_bold, next_entry, end_bold)
    }
    row_m <- paste0(row_m, " \\\\")
    write(row_m, file_path, append = T)
    if (m %in% sep_after) {
      write("\\midrule", file_path, append = T)
    }
  }
  write("\\bottomrule", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}

# Table 1: Overall quadratic and Poisson score, and specific elementary scores -
t_dom <- matrix(NA, nrow = length(models), ncol = 5)
rownames(t_dom) <- model_names
colnames(t_dom) <- c("Poisson", "Quadratic", "Elem1", "Elem2", "Elem3")

for (i in 1:length(models)) {
  t_dom[i, "Quadratic"] <- s_quad_da(models[[i]], obs)
  t_dom[i, "Poisson"] <- s_pois_da(models[[i]], obs)
  t_dom[i, "Elem1"] <- s_theta_da(models[[i]], obs, exp(-9.5))
  t_dom[i, "Elem2"] <- s_theta_da(models[[i]], obs, exp(-11))
  t_dom[i, "Elem3"] <- s_theta_da(models[[i]], obs, exp(-14))
}
t_dom
write.csv(t_dom, file.path(fpath, "Tab1_dom.csv"))

col_names <- list(
  c("Score", "Poisson", "Quadratic", "Elementary", "Elementary", "Elementary"),
  c("", "", "", "$\\theta = 10^{-9.5}$", "$\\theta = 10^{-11}$", "$\\theta = 10^{-14}$")
)
print_tex_table(
  t_dom,
  col_names = col_names,
  digits = c(2, 4, 4, 4, 4),
  make_strong = rownames(t_dom)[apply(t_dom, 2, which.min)],
  strong_symbol = "\\green",
  sep_after = "LG",
  col_align = "r",
  file_path = file.path(fpath, "Tab1_dom.tex")
)

# Table 3: Overall quadratic and Poisson score and its number and spatial component
t_numspat <- matrix(NA, nrow = length(models), ncol = 4)
rownames(t_numspat) <- model_names
colnames(t_numspat) <- c("quad", "pois", "number", "spatial")

for (i in 1:length(models)) {
  x_t <- rowSums(models[[i]])
  t_numspat[i, "quad"] <- s_quad_da(models[[i]], obs)
  t_numspat[i, "pois"] <- s_pois_da(models[[i]], obs)
  t_numspat[i, "number"] <- mean(s_pois(x_t, rowSums(obs)))
  t_numspat[i, "spatial"] <- mean(rowSums(s_pois(models[[i]] / x_t, obs))) - 1
}
t_numspat
write.csv(t_numspat, file.path(fpath, "Tab4_number-spatial.csv"))

print_tex_table(
  t_numspat,
  col_names = list(c("", "Quadratic", "Poisson", "Number", "Spatial")),
  digits = c(4, 2, 3, 3),
  make_strong = rownames(t_numspat)[apply(t_numspat, 2, which.min)],
  sep_after = "LG",
  file_path = file.path(fpath, "Tab4_number-spatial.tex")
)

# Table 4: Overall quadratic and Poisson score and its MSB, DSC, and UNC component
t_sc_cmps <- matrix(NA, nrow = length(models), ncol = 8)
rownames(t_sc_cmps) <- model_names
colnames(t_sc_cmps) <- c("quad", "q-MCB", "q-DSC", "q-UNC", "pois", "p-MCB", "p-DSC",
                         "p-UNC")

# we could get different scores than Table 1 due to sorting and summing up in a
# different order, but does not seem to be the case here

for (i in 1:length(models)) {
  # recalibrate forecasts x with isotoinc regression from monotone package
  x <- as.vector(models[[i]])
  y <- as.vector(obs)
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  j <- 1
  for (scf in list(s_quad_da, s_pois_da)) {
    s <- scf(x, y)
    s_rc <- scf(x_rc, y)
    s_mg <- scf(mean(y), y)
    t_sc_cmps[i, j] <- s
    t_sc_cmps[i, j + 1] <- s - s_rc
    t_sc_cmps[i, j + 2] <- s_mg - s_rc
    t_sc_cmps[i, j + 3] <- s_mg
    j <- j + 4
  }
}
t_sc_cmps
write.csv(t_sc_cmps, file.path(fpath, "Tab5_score-components.csv"))

make_bold <- rownames(t_sc_cmps)[apply(t_sc_cmps * rep(c(1, 1, -1, 1, 1, 1, -1, 1),
                                                       each = nrow(t_sc_cmps)),
                                       2, which.min)]
make_bold[c(4, 8)] <- "XXXXX" # in column 4 and 8 (UNC), make no Model bold
print_tex_table(
  t_sc_cmps,
  col_names = list(c("", "quad", "MCB", "DSC", "UNC", "pois", "MCB", "DSC", "UNC")),
  digits = c(rep(4, 4), rep(2, 4)),
  make_strong = make_bold,
  sep_after = "LG",
  file_path = file.path(fpath, "Tab5_score-components.tex"))

rm(t_dom, t_numspat, t_sc_cmps, col_names, make_bold, x, y, ord, x_rc, s, s_rc, s_mg, j)

################################################################################
# Statistical tests
################################################################################

csep_test <- function(fcst1, fcst2, y, scf) {
  N <- sum(y)
  diff_logs <- log(fcst1) - log(fcst2)
  I_N_ij <- sum(y * diff_logs) / N - (sum(fcst1) - sum(fcst2)) / N
  test_var <- 1 / (N - 1) * sum((diff_logs - mean(diff_logs))^2)

  test_stat <- I_N_ij / sqrt(test_var) * sqrt(N)
  pval <- 1 - pt(test_stat, N - 1)
  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(test_var), mean_diff = I_N_ij))
}

dm_test <- function(fcst1, fcst2, y, scf) {
  diff_scores <- rowSums(scf(fcst1, y) - scf(fcst2, y))
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = 7, type = "covariance", plot = F)$acf
  dm_var <- auto_covs[1] + 2 * sum(auto_covs[-1])

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(n_days)
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var),
                    mean_diff = mean_diff_score))
}

# 1. DM with Poisson -----------------------------------------------------------
test_func <- dm_test
scf_func <- s_pois
scale <- 1.0
name <- "dm_pois"
n_digits <- 2
model_order <- c(5, 1, 4, 2, 3)
# 2. DM with Quadratic ---------------------------------------------------------
test_func <- dm_test
scf_func <- s_quad
scale <- 1.0
name <- "dm_quad"
n_digits <- 3
model_order <- c(1, 5, 4, 2, 3)
# 3. CSEP T-test ---------------------------------------------------------------
test_func <- csep_test
scf_func <- s_pois
scale <- n_days / sum(obs)
name <- "csep"
n_digits <- 3
# ------------------------------------------------------------------------------

scores <- sapply(models, function(x) sum(scf_func(x, obs)) / n_days)
val_matrix <- matrix(0.0, nrow = length(models), ncol = length(models),
                     dimnames = list(model_names[model_order], model_names[model_order]))

for (i in 1:(length(model_order) - 1)) {
  # diagonal corresponds to Poisson score
  val_matrix[i, i] <- scale * scores[model_order[i]]

  for (j in (i + 1):length(model_order)) {
    test_vals <- test_func(models[[model_order[i]]], models[[model_order[j]]], obs, scf_func)

    val_matrix[j, i] <- test_vals$pval
    val_matrix[i, j] <- test_vals$zval
  }
}
val_matrix[length(models), length(models)] <- scale * scores[model_order[length(models)]]

write.csv(val_matrix, file.path(fpath, paste0("tests_", name, ".csv")))

print_tex_table(
  val_matrix,
  col_names = colnames(val_matrix),
  digits = rep(n_digits, length(models)),
  make_strong = rownames(val_matrix),   # just emphasize diagonals
  strong_symbol = "\\fbox",
  file_path = file.path(fpath, paste0("tests_", name, ".tex"))
)

################################################################################
# Visualize Poisson score (differences) temporally
################################################################################

# 1 - Poisson score ------------------------------------------------------------
scf <- s_pois
add_name <- ""
add_title <- "Poisson"
# transformation for difference plot
d <- 10^3
my_trans <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d  + 1),
  function(y) sign(y) / d * (exp(abs(y)) - 1)
)
my_breaks <- c(-10^(c(2, 0, -2)), 0, 10^(c(-2, 0, 2)))
minor_breaks <- c(-10^(2:-3), 0, 10^(-3:2))
my_labels <- c("-100", "-1", "-0.01", "0", "0.01", "1", "100")
my_limits <- c(-100, 100)

# 2 - Quadratic score ----------------------------------------------------------
scf <- s_quad
add_name <- "_quad"
add_title <- "Quadratic"
# transformation for difference plot
d <- 10^6
my_trans <- trans_new(
  "log", function(x) sign(x) * log(abs(x) * d  + 1),
  function(y) sign(y) / d * (exp(abs(y)) - 1)
)
my_breaks <- c(-10^(c(0, -2, -4)), 0, 10^(c(-4, -2, 0)))
minor_breaks <- c(-10^(1:-7), 0, 10^(-7:1))
my_labels <- my_breaks
my_limits <- c(-10, 10)

# ------------------------------------------------------------------------------

scores <- do.call(cbind, lapply(models, function(X) rowSums(scf(X, obs))))
scores <- data.frame(scores) %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(scores) <- c(model_names, "X", "earthquake")
scores_long <- pivot_longer(scores, cols = all_of(model_names), names_to = "Model")

point_size <- 0.75
point_alpha <- 0.4

score_plot <- ggplot() +
  geom_point(data = filter(scores_long, !earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  # add black borders to triangle used for earthquake days
  geom_point(data = filter(scores_long, earthquake), aes(x = X, y = value),
             color = "black", shape = 2, size = point_size * 1.2, stroke = 0.1,
             alpha = point_alpha) +
  # plot triangles
  geom_point(data = filter(scores_long, earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores))) +
  scale_color_manual(name = "", values = model_colors,
                     guide = guide_legend(order = 1, direction = "horizontal",
                                          override.aes = list(alpha = 1))) +
  scale_shape_manual(name = "earthquake", labels = c("No", "At least one"),
                     values = c("FALSE" = 16, "TRUE" = 17),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 1),
                                          direction = "vertical", title.position = "right")) +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Score") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center")

combine_plots <- grid.arrange(score_plot, nrow = 1,
                              top = textGrob(paste(add_title, "Score by Day"),
                                             gp = gpar(fontsize = title_size)))
file_path <- file.path(fpath, paste0("Fig3_DailyScores", add_name, ".pdf"))
ggsave(file_path, width = 145, height = 83, unit = "mm", plot = combine_plots)

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

temp_plot <- ggplot() +
  # plot points going with no earthquakes first!
  # (that is why we separated geom point in two to define order of drawing groups)
  geom_point(data = filter(diff_scores, !earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  # add black borders to triangle used for earthquake days
  geom_point(data = filter(diff_scores, earthquake), aes(x = X, y = value),
             color = "black", shape = 2, size = point_size * 1.2, stroke = 0.1,
             alpha = point_alpha) +
  # plot triangles
  geom_point(data = filter(diff_scores, earthquake),
             aes(x = X, y = value, color = Model, shape = earthquake),
             size = point_size, alpha = point_alpha) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = scores$X[new_year], labels = year(times[new_year]),
                     limits = c(0, nrow(scores))) +
  scale_y_continuous(trans = my_trans, breaks = my_breaks, labels = my_labels,
                     minor_breaks = minor_breaks, limits = my_limits) +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors, breaks = ana_models,
                     guide = guide_legend(order = 1, direction = "horizontal",
                                          override.aes = list(alpha = 1))) +
  scale_shape_manual(name = "earthquake", labels = c("No", "At least one"),
                     values = c("FALSE" = 16, "TRUE" = 17),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 1),
                                          direction = "vertical", title.position = "right")) +
  xlab(NULL) +
  ylab("Difference") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center")

combine_plots <- grid.arrange(temp_plot, nrow = 1,
                              top = textGrob(paste(add_title, "Score Difference by Day"),
                                             gp = gpar(fontsize = title_size)))
file_path <- file.path(fpath, paste0("Fig4_ScoreDiffTemp", add_name, ".pdf"))
ggsave(file_path, width = 145, height = 83, unit = "mm", plot = combine_plots)

rm(scores, scores_long, diff_scores, combine_plots, score_plot, temp_plot, my_trans)

################################################################################
# Poisson Score: Number and spatial component
################################################################################

# score of predicted number of earthquakes in all of Italy
scores_n <- do.call(
  cbind, lapply(models, function(X) s_pois(rowSums(X), rowSums(obs)))
)
scores_n <- data.frame(scores_n)
colnames(scores_n) <- model_names

diff_scores_n <- scores_n %>%
  mutate(X = 1:nrow(.)) %>%
  filter(rowSums(obs) > 0) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

# score of normalized mean predicted number of earthquakes
# divide every forecast by daily aggregated forecast (row sum)
scores_s <- do.call(
  cbind, lapply(models, function(X) {
    row_sums <- rowSums(X)
    X_norm <- apply(X, 2, function(col) col / row_sums)
    return(rowSums(s_pois(X_norm, obs)) - 1)
  })
)
scores_s <- data.frame(scores_s)
colnames(scores_s) <- model_names

diff_scores_s <- scores_s %>%
  mutate(X = 1:nrow(.)) %>%
  filter(rowSums(obs) > 0) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

number_spatial_plot <- ggplot(rbind(cbind(diff_scores_n, C = "Number component"),
                                    cbind(diff_scores_s, C = "Spatial component"))) +
  facet_wrap(~C, nrow = 2, scales = "free_y") +
  geom_point(aes(x = X, y = value, color = Model), size = 0.75, alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_x_continuous(breaks = which(new_year), labels = year(times[new_year]),
                     limits = c(0, nrow(scores_s))) +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors,
                     breaks = ana_models,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 0.75)))+
  scale_y_continuous(limits = c(-20, NA)) +
  xlab(NULL) +
  ylab(NULL) +
  my_theme +
  theme(legend.position = "bottom")

combine_plots <- grid.arrange(number_spatial_plot,
                              top = textGrob("Poisson Score Difference by Day",
                                             gp = gpar(fontsize = title_size)),
                              left = textGrob("Difference", rot = 90,
                                              gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "Fig5_NumberSpatials.pdf")
ggsave(file_path, width = 145, height = 135, unit = "mm", plot = combine_plots)

# check whether plots fits to overall number and spatial score
colMeans(scores_n)
colMeans(scores_s)

# plot against number of observed earthquakes
diff_scores_n <- scores_n %>%
  mutate(Y = rowSums(obs)) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(Y, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")
diff_scores_s <- scores_s %>%
  mutate(Y = rowSums(obs)) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(Y, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

number_spatial_plot <- ggplot(rbind(cbind(diff_scores_n, C = "Number component"),
                                    cbind(diff_scores_s, C = "Spatial component"))) +
  facet_wrap(~C, nrow = 2, scales = "free_y") +
  geom_point(aes(x = jitter(Y, amount = 0.2), y = value, color = Model), size = 0.75,
             alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors,
                     breaks = ana_models,
                     guide = guide_legend(override.aes = list(alpha = 1, size = 0.75)))+
  # scale_y_continuous(limits = c(-20, NA)) +
  ylab(NULL) +
  xlab("Number of obs. earthquakes") +
  my_theme +
  theme(legend.position = "bottom")

combine_plots <- grid.arrange(number_spatial_plot,
                              top = textGrob("Poisson Score Difference by Day",
                                             gp = gpar(fontsize = title_size)),
                              left = textGrob("Difference", rot = 90,
                                              gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "Fig5_NumberSpatials_2.pdf")
ggsave(file_path, width = 145, height = 135, unit = "mm", plot = combine_plots)

rm(scores_n, diff_scores_n, scores_s, diff_scores_s, number_spatial_plot, combine_plots)

################################################################################
# Visualize Poisson score differences spatially
################################################################################

# use log transform for positive and negative values (see ?modulus_trans)
# but need to scale with d so that log scaling gets active (around zero the
# transformation is the identity, but for large absolute values it is a log
# transform

# 1 - Poisson score ------------------------------------------------------------
scf <- s_pois
add_name <- ""
add_title <- "Poisson"
# transformation for difference plot
d <- 10^5
my_trans <- function(x) sign(x) * log(abs(x) * d  + 1)
my_breaks <- c(-0.01, -0.001, -0.0001, 0, 0.0001, 0.001)
my_labels <- expression(-10^-2, "", -10^-4, 0, 10^-4, "")

# 2 - Quadratic score ------------------------------------------------------------
scf <- s_quad
add_name <- "_quad"
add_title <- "Quadratic"
# transformation for difference plot
d <- 10^8
my_trans <- function(x) sign(x) * log(abs(x) * d  + 1)
my_breaks <- c(-10^c(-3, -4, -5, -6 , -7), 0, 10^c(-7, -6, -5, -4))
my_labels <- expression(-10^-3, "", "", -10^-6, "", 0, "", 10^-6, "", "")

# ------------------------------------------------------------------------------

scores <- do.call(cbind, lapply(models, function(X) colMeans(scf(X, obs))))
scores <- data.frame(scores)
colnames(scores) <- model_names

diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  mutate(LON = cells$LON, LAT = cells$LAT) %>%
  select(LON, LAT, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")
diff_scores$Model <- factor(paste(cmp_model, "vs.", diff_scores$Model), ordered = T,
                            levels = paste(cmp_model, "vs.", names(model_colors)))

my_colors <- c("#800303", "#f51818", "#ffffff", "#057ffa")
limits <- my_trans(range(diff_scores$value))
col_breaks <- c(limits[1], -limits[2], 0, limits[2])

eq_loc <- events %>%
  group_by(N) %>%
  summarise(Count = n(), .groups = "drop")  %>%
  left_join(cells, by = "N") %>%
  select(LON, LAT, Count)

spat_plot <- ggplot() +
  facet_wrap(~Model, nrow = 2) +
  geom_tile(data = diff_scores,
            aes(x = LON, y = LAT, fill = my_trans(value)), alpha = 0.5) +
  geom_tile(data = eq_loc, aes(x = LON, y = LAT, color = "Obs. earthquakes"),
            fill = NA) +
  geom_sf(data = filter(europe, name == "Italy"), color = "black", fill = NA,
          size = 0.2) +
  coord_sf(xlim = lon_lim, ylim = lat_lim, expand = TRUE) +
  scale_x_continuous(name = NULL, breaks = c(6, 10, 14, 18)) +
  scale_y_continuous(name = NULL, breaks = c(36, 40, 44, 48)) +
  scale_fill_gradientn(name = "Score\ndifference",
                       colors = my_colors, values = rescale(col_breaks),
                       breaks = my_trans(my_breaks), labels = my_labels, limits = limits,
                       guide = guide_colorbar(barwidth = unit(40, "mm"),
                                              title.vjust = 0.9, order = 1)) +
  scale_color_manual(name = "Obs.\nearthquakes", values = c("Obs. earthquakes" = "black"),
                     labels = "",
                     guide = guide_legend(keywidth = unit(5, "points"),
                                          keyheight = unit(5, "points"))) +
  my_theme +
  theme(legend.position = "bottom")

combine_plots <- grid.arrange(
  spat_plot, nrow = 1,
  top = textGrob(paste("Average", add_title, "Score Difference by Grid Cell"),
                 gp = gpar(fontsize = title_size))
)
file_path <- file.path(fpath, paste0("Fig6_ScoreDiffSpat", add_name, ".pdf"))
ggsave(file_path, width = 145, height = 160, unit = "mm", plot = combine_plots)

rm(scores, diff_scores, combine_plots, spat_plot, eq_loc)

################################################################################
# Patton diagramm
################################################################################

s_b <- function(x, y, b) {
  b_1 <- b - 1
  if (b == 0) {
    qu <- y / x
    s <- qu - log(qu) - 1
  } else if (b == 1) {
    s <- y * log(y / x) - y + x
  } else {
    s <- 1 / (b * b_1) * (y^b - x^b) - 1 / b_1 * x^b_1 * (y - x)
  }
  return(sum(s) / n_days)
}

# if s_b is the Patton score with parameter b this modified version corresponds to
# s_b_mod(x, y) = s_b(x, y) + s_b(1, y) + 0.5 y^abs(b) - b/2 y + (3-b)/2
# whereby the evaluation was simplified
s_b_mod <- function(x, y, b) {
  if (b == 0) {
    s <- y / x - y + log(x) + 2
  } else if (b == 1) {    # Poisson score
    s <- (-1) * y * log(x) + x
  } else {
    b_1 <- b - 1
    s <- 1 / (b * b_1) * (1 - x^b) - 1 / b_1 * x^b_1 * (y - x) + 0.5 * y^abs(b) +
      (1 / b_1 - b / 2) * y + (3 - b) / 2 - 1 / b_1
  }
  return(sum(s) / n_days)
}

s_patt <- compiler::cmpfun(s_b_mod) # compile function to reduce runtime a bit

n_b <- 100
bs <- seq(-2, 8, length.out = n_b)

patton_df <- matrix(0, nrow = length(bs), ncol = length(models))
for (m in 1:length(models)) {
  print(m)
  for (t in 1:n_b) {
    cat("*")
    patton_df[t, m] <- s_patt(models[[m]], obs, bs[t])
  }
}
colnames(patton_df) <- model_names

patton_df <- data.frame(patton_df) %>% mutate(b = bs)

patton_df <- read.csv("./../tmp_results/patton_df_100.csv")

max_vals <- apply(patton_df[, model_names], 1, max)
min_vals <- apply(patton_df[, model_names], 1, min)

patton_plot <- patton_df %>%
  mutate(across(all_of(model_names),
                 function(vec) (vec - min_vals) / (max_vals - min_vals))) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = b, y = value, color = Model), size = 0.3) +
  geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "darkgray", size = 0.5) +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  xlab("b") +
  ylab("Mean score") +
  annotate("text", x = c(1.3, 2.3), y = c(0.8, 0.7), label = c("b=1", "b=2"),
           size = 3, color = "darkgray") +
  my_theme +
  theme(legend.position = "bottom")

combine <- grid.arrange(patton_plot, nrow = 1,
                        top = textGrob("Mean Score of Patton Scoring Functions",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig10_PattonPlot_4.pdf")
ggsave(file_path, width = 145, height = 75, unit = "mm", plot = combine)

# patton difference plot
patton_plot <- patton_df %>%
  filter(b >= 0.9) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = b, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "darkgray", size = 0.5) +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors,
                     breaks = ana_models,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  xlab("b") +
  ylab("Difference") +
  annotate("text", x = c(1.3, 2.3), y = c(-0.4, -0.4), label = c("b=1", "b=2"),
           size = 3, color = "darkgray") +
  my_theme +
  theme(legend.position = "bottom")

combine <- grid.arrange(patton_plot, nrow = 1,
                        top = textGrob("Score Difference of different Patton Scores",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig10_PattonPlot_4_diff.pdf")
ggsave(file_path, width = 145, height = 75, unit = "mm", plot = combine)


################################################################################
# Murphy diagramm
################################################################################

S_elem <- compiler::cmpfun(s_theta_da) # compile function to reduce runtime a bit

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

# or read it in
murphy_df <- read.csv("./figures7/murphy_df.csv") %>%
  select(all_of(model_names))

x_midpoints <- (log_grid[-1] + log_grid[-n_theta]) / 2
best_models <- apply(as.matrix(murphy_df), 1, which.min)
best_models_bar <- data.frame(Model = colnames(murphy_df)[best_models],
                              xmin = c(log_grid[1], x_midpoints),
                              xmax = c(x_midpoints, log_grid[n_theta]),
                              ymin = 0.32, ymax = 0.325)

murphy_diag <- data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  geom_rect(data = best_models_bar, aes(color = Model, fill = Model, xmin = xmin,
                                        xmax = xmax, ymin = ymin, ymax = ymax),
            show.legend = F) +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  scale_fill_manual(name = NULL, values = model_colors) +
  xlab(expression(paste("Threshold log", (theta)))) +
  ylab("Elementary score") +
  my_theme +
  theme(legend.position = c(0.01, 0.1), legend.justification = c(0, 0))

combine <- grid.arrange(murphy_diag, nrow = 1,
                        top = textGrob("Logarithmic Murphy Diagram",
                                       gp = gpar(fontsize = title_size)))

# murphy difference diagram
murphy_diag <- data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(theta, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors,
                     breaks = ana_models,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  xlab(expression(paste("Threshold log", (theta)))) +
  ylab("Elementary score") +
  my_theme +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0))

combine <- grid.arrange(murphy_diag, nrow = 1,
                        top = textGrob("Murphy Difference Diagram",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig_LogMurphyDiag.pdf")
ggsave(file_path, width = 145, height = 75, unit = "mm", plot = combine)

# now look at Murphy diagram of miscalibration and discrimination component
MCB_diag <- DSC_diag <- matrix(0.0, ncol = n_mods, nrow = n_theta)
for (i in 1:n_mods) {
  print(i)
  # cell_decomposition determines score components for each spatial cell separately
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
  # again divide by n_days to get daily averages
  cbind(data.frame(MCB_diag / n_days), Type = "MCB"),
  cbind(data.frame(DSC_diag / n_days), Type = "DSC")
)
colnames(df_collect) <- c(model_names, "Type")
UNC_diag <- UNC_diag / n_days

df_collect <- df_collect %>%
  mutate(log_theta = rep(log_grid, 2)) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  mutate(Type = ifelse(Type == "MCB", "Miscalibration", "Discrimination"))

# or read precomputed values in
df_collect <- read.csv("./figures7/murphy-MCB-DSC.csv", row.names = 1) %>%
  filter(Model %in% model_names)
murphy_df <- read.csv("./figures7/murphy_df.csv") %>%
  select(all_of(model_names)) %>%
  mutate(log_theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  mutate(Type = "Score")

mcb_fac <- 5
df_plot <- rbind(murphy_df, df_collect) %>%
  mutate(value = ifelse(Type == "Discrimination", value * (-1), value),
         value = ifelse(Type == "Miscalibration", value * mcb_fac, value))

murphy_score_cmps <- ggplot(df_plot) +
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_line(aes(x = log_theta, y = value, color = Model, group = paste(Model, Type)),
            size = 0.3) +
  scale_color_manual(name = NULL, values = model_colors,
                     guide = guide_legend(override.aes = list(size = 0.5))) +
  scale_x_continuous(breaks = -4:1 * 5) +
  scale_y_continuous(breaks = 0:3 * 0.1, minor_breaks = -4:6 * 0.05,
                     sec.axis = sec_axis(~./mcb_fac, name = NULL,
                                         breaks = c(-2:0 * 0.02, 1:3 * 0.01),
                                         labels = c(2:0 * 0.02 * mcb_fac, 1:3 * 0.01))) +
  xlab(expression(paste("Threshold log", (theta)))) +
  ylab(NULL) +
  my_theme +
  annotate("text", x = -Inf, y = 0.12, label = "Average score", angle = 90, vjust = 2) +
  annotate("text", x = Inf, y = -0.12, label = "DSC", angle = -90, vjust = 2) +
  annotate("text", x = Inf, y = 0.09, label = "MCB", angle = -90, vjust = 2) +
  theme(legend.position = "bottom", legend.key.size = unit(4, "mm"))

combine <- grid.arrange(murphy_score_cmps, nrow = 1,
                        top = textGrob("Score Components by Elementary Score",
                                       gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig9_Murphy-MCB-DSC.pdf")
ggsave(file_path, width = 145, height = 110, unit = "mm", plot = combine)

# coherence checks:
# If we integrate murphy diagram with correct measure, do we get score?
# If we integrate MCB / DSC Murphy diagram, do we get MCB and DSC component?
# If we combine MCB, DSC, UNC Murphy diagram, do we get Murphy diagram?
combine_back <- df_collect %>%
  pivot_wider(id_cols = c(log_theta, Model), names_from = Type, values_from = value) %>%
  mutate(MSC_DSC = Miscalibration - Discrimination) %>%
  pivot_wider(id_cols = log_theta, names_from = Model, values_from = MSC_DSC) %>%
  mutate(across(all_of(model_names), function(m) m + UNC_diag)) %>%
  select(all_of(model_names))
colSums(abs(combine_back - murphy_df))

rm(murphy_df, decomp, murphy_diag, df_collect, murphy_score_cmps, combine,
   MCB_diag, DSC_diag, UNC_diag, combine_back)

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
                        sprintf("%.3f", res$stats[1, ]),
                        collapse = "\n"))
  )
}

# or load already recalibrated values
recal_models <- read.csv("./figures7/recal-all_Til-100.csv", row.names = 1) %>%
  filter(Model %in% model_names)
collect_stats <- read.csv("./figures7/collect_stats.csv", row.names = 1) %>%
  filter(Model %in% model_names)

col_ecdfs <- list()
for (i in 1:length(models)) {
  t <- table(as.vector(models[[i]]))
  col_ecdfs[[i]] <- data.table(x = c(0, as.numeric(names(t))),
                               y = c(0, cumsum(as.numeric(t))) / prod(dim(models[[i]])),
                               M = model_names[i])
}

# or just loaded already computed values
for (i in 1:length(models)) {
  col_ecdfs[[i]] <- read.csv(paste0("./../tmp_results/ecdf_", model_names[i], ".csv"), row.names = 1)
}

my_labeller <- function(l) {
  labels <- character(length(l))
  labels[l > 0 & l < 1] <- paste0("1e", log(l[l > 0 & l < 1], base = 10))
  labels[l == 1] <- "1"
  labels[l == 0] <- "0"
  return(labels)
}

# 1 - use averaed empirical CDF transform --------------------------------------
mean_ecdf <- do.call(rbind, col_ecdfs) %>%
  pivot_wider(id_cols = x, names_from = M, values_from = y, values_fill = NA) %>%
  arrange(x) %>%
  setnafill(type = "locf") %>%
  slice(floor(seq(1, nrow(.), length.out = 10^5))) %>%   # sample on grid for plotting
  transmute(x = x, y = rowMeans(select(., all_of(model_names))))

my_trans <- function(x) {
  # first known values, then values we want to evaluate --> fill with last observation
  rbind(cbind(mean_ecdf, a = -1.0), data.table(x = x, y = NA, a = 1:length(x))) %>%
    arrange(x) %>%
    setnafill(type = "locf") %>%
    filter(a > 0) %>%
    arrange(a) %>%
    pull(y)
}
# position of inset histograms
xmin <- 0.7
xmax <- 1.0
ymin <- 0.0
ymax <- 0.25
# axis breaks and labels
breaks <- c(0, 10^c(-7, -6, -5, -4, 0))
t_breaks <- my_trans(breaks)
labels <- c("0", rep("", length(breaks) - 2), "1")
# position of score components
text_x <- 0.02
text_y <- 0.64
# add to name
add_name <- ""

# 2 - use two sided log transform ----------------------------------------------
d <- 10^9
my_trans <- function(x) sign(x) * log(abs(x) * d  + 1)
# position of inset histograms
xmin <- 13.8
xmax <- 20
ymin <- 0.0
ymax <- 6.9
# axis breaks and labels
breaks <- c(0, 10^c(-8, -6, -4, -2, 0))
t_breaks <- my_trans(breaks)
labels <- rep("", length(breaks))
# position of score components
text_x <- 0.02
text_y <- 15
# add to name
add_name <- "_log"

# 3 - use no / standard transform ----------------------------------------------
my_trans <- function(x) x
# position of inset histograms
xmin <- 1.6
xmax <- 2.2
ymin <- 0.05
ymax <- 0.6
# axis breaks and labels
t_breaks <- seq(0, 2, length.out = 6)
labels <- c("0", rep("", length(breaks) - 2), 2)
# position of score components
text_x <- 1.35
text_y <- 1.40
# add to name
add_name <- "_std"
# cut confidence band at plot max
recal_models$upper <- pmin(recal_models$upper, 2.2)

# ------------------------------------------------------------------------------


inset_histograms <- list()
plot_min <- my_trans(min(c(recal_models$x, recal_models$lower)))
plot_max <- my_trans(max(c(recal_models$x, recal_models$upper)))
hist_breaks <- seq(plot_min, plot_max, length.out = 9)
for (i in 1:length(models)) {
  my_hist <- ggplot(data.table(x = as.vector(models[[i]]))) +
    geom_histogram(aes(x = my_trans(x)), fill = "gray", col = "black", size = 0.2,
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

# specify rows separately to center align bottom row
rows <- list(model_names[1:3], model_names[4:5])

create_row <- function(row) {
  dt_recal <- filter(recal_models, Model %in% row)
  dt_stats <- filter(collect_stats, Model %in% row)

  return(
    ggplot(dt_recal, aes(x = my_trans(x))) +
      facet_wrap(~factor(Model, ordered = T, levels = row), nrow = 1) +
      geom_ribbon(aes(ymin = my_trans(lower), ymax = my_trans(upper), fill = Model),
                  alpha = 0.33, show.legend = FALSE) +
      geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
                  linetype = "dashed") +
      geom_step(aes(y = my_trans(x_rc), color = Model), size = 0.3, show.legend = FALSE) +
      scale_color_manual(values = model_colors) +
      scale_fill_manual(values = model_colors) +
      scale_x_continuous(breaks = t_breaks, labels = labels, limits = c(plot_min, plot_max)) +
      scale_y_continuous(breaks = t_breaks, labels = labels, limits = c(plot_min, plot_max)) +
      xlab(NULL) +
      ylab(NULL) +
      ggtitle(NULL) +
      geom_text(data = dt_stats, mapping = aes(x = text_x, y = text_y, label = label),
                size = 7 * 0.36, hjust = 0, vjust = 0) +
      my_theme +
      theme(aspect.ratio = 1, panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), plot.margin = margin(5.5, 3.5, 5.5, 3.5))
  )
}

combine <- grid.arrange(create_row(rows[[1]]) + inset_histograms[which(model_names %in% rows[[1]])],
                        create_row(rows[[2]]) + inset_histograms[which(model_names %in% rows[[2]])],
                        top = textGrob("Reliability Diagram",
                                       gp = gpar(fontsize = title_size)),
                        bottom = textGrob("Forecasted mean",
                                       gp = gpar(fontsize = 11)),
                        left = textGrob("Conditional mean", rot = 90,
                                       gp = gpar(fontsize = 11)))

# create small hist to describe functionality of inset histogram pictogram like
small_hist <- ggplot(data.table(x = rep(hist_breaks[-1], 1:(length(hist_breaks) - 1)))) +
    geom_histogram(aes(x = x), fill = "gray", col = "black", size = 0.2,
                   breaks = hist_breaks, alpha = 0.5) +
    theme_classic(base_size = 5.5) +
    theme(axis.line.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.background = element_blank(),
          panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

finish <- ggplot() +
  # set axis limits
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
  # add reliability diagram
  annotation_custom(combine, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  # add small histogram on the side with describing text
  annotation_custom(ggplotGrob(small_hist), xmin = 0.88, xmax = 0.96, ymin = 0.07, ymax = 0.15) +
  annotate("text", x = 0.92, y = 0.16, label = "Histogram",
           size = 7 * 0.36, lineheight = 0.7, color = "darkgray") +
  annotate("text", x = 0.92, y = 0.06, label = "Forecasted\nmean",
           size = 7 * 0.36, lineheight = 0.7, color = "darkgray") +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        panel.background = element_rect("white"),
        plot.margin = margin(0, 0, 0, 0))

file_path <- file.path(fpath, paste0("Fig8_ReliabilityDiagram", add_name, ".pdf"))
ggsave(file_path, width = 145, height = 125, unit = "mm", plot = finish)

rm(col_ecdfs, recal_models, collect_stats, combine, finish, small_hist,
   inset_histograms, mean_ecdf, my_hist)
