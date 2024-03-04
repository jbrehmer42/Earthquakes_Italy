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
library(geomtextpath)     # for geom_labelabline: draw lines with label

source("data_prep.R")

# use standard colors of ggplot for discrete variables
model_colors <- c("FCM" = "#F8766D", "LG" = "#00BF7D",  "LM" = "#A3A500",
                  "SMA" = "#00B0F6", "LRWA" = "#E76BF3")
model_order <- 1:5

cmp_model <- "LM"
cmp_m <- sym(cmp_model)
ana_models <- model_names[model_names != cmp_model]

new_year <- month(times) == 1 & day(times) == 1 & year(times) %% 2 == 0

fpath <- "./figures9"

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

inf_gain <- function(X1, X2, Y) {
  X1_t <- rowSums(X1)
  X2_t <- rowSums(X2)
  return(
    rowSums(Y * (log(X1) - log(X2))) - (X1_t - X2_t)
  )
}

cum_inf_gain <- function(X1, X2, Y) {
  return(cumsum(inf_gain(X1, X2, Y)))
}

inf_gain_per_event <- function(X1, X2, Y) {
  n_t <- rowSums(Y)
  X1_t <- rowSums(X1)
  X2_t <- rowSums(X2)
  return(
    rowSums(Y * (log(X1) - log(X2))) / n_t - (X1_t - X2_t) / n_t
  )
}

cum_igpe <- function(X1, X2, Y) {
  n_t <- cumsum(rowSums(Y))
  X1_t <- cumsum(rowSums(X1))
  X2_t <- cumsum(rowSums(X2))
  return(
    cumsum(rowSums(Y * (log(X1) - log(X2)))) / n_t - (X1_t - X2_t) / n_t
  )
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
file_path <- file.path(fpath, "Fig1_Earthquakes.pdf")
ggsave(file_path, width = 140, height = 100, unit = "mm", plot = combine)

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
             aes(x = X, y = 0.15, shape = "Observed earthquakes"), size = 0.75,
             color = "gray") +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = pred_by_day$X[new_year], labels = year(times[new_year])) +
  scale_y_log10() +
  scale_color_manual(name = NULL, values = model_colors[model_names[model_order]],
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(name = NULL, values = c("Observed earthquakes" = 1)) +
  xlab(NULL) +
  ylab("Expected number") +
  labs(subtitle = "Entire CSEP Italy region") +
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

create_row <- function(row, only_legend = F, top = F) {
  pred_dt <- filter(pred_by_cell_long, Model %in% row)
  top_margin <- ifelse(top, 5.5, 0.0)
  bottom_margin <- ifelse(top, 0.0, 5.5)
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
    scale_color_manual(name = "Observed\nearth-\nquakes", values = c("Obs. earthquakes" = "black"),
                       labels = "", guide = guide_legend(order = 2, override.aes = list(size = 4))) +
    my_theme +
    theme(legend.position = "right",
          plot.margin = margin(top_margin, 5.5, bottom_margin, 11.5))
  if (only_legend) {
    return(get_legend(pl + guides(size = "none")))        # only show legend
  } else {
    return(pl + guides(fill = "none", color = "none", size = "none"))    # don't show legend
  }
}
print(paste("Check:", times[i_time], "= 03.04.2009?"))
# use short dash ("en dash") = \u2013
spat_subtitle <- paste("7\u00adday Period Starting April 3, 2009*")
spat_plot <- grid.arrange(create_row(c("LM", "FCM", "LG"), top = T),
                          create_row(c("SMA", "LRWA"), top = F),
                          nrow = 2)
spat_plot <- grid.arrange(spat_plot, create_row(model_names, only_legend = T), nrow = 1,
                          widths = c(0.85, 0.15),
                          top = textGrob(spat_subtitle, x = unit(0.09, "npc"),
                                         hjust = 0, gp = gpar(fontsize = 11)))

combine <- grid.arrange(temp_plot, spat_plot, nrow = 2, heights = c(0.37, 0.63))

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
write.csv(t_dom, file.path(fpath, "Tab1_scores.csv"))

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
  file_path = file.path(fpath, "Tab1_Scores.tex")
)

# Table 2: Poisson, Information Gain and Information Gain per event ------------
t_sum <- matrix(NA, nrow = length(models), ncol = 4)
rownames(t_sum) <- model_names
colnames(t_sum) <- c("Poisson", "Pois-diff", "IG", "IGPE")

for (i in 1:length(models)) {
  t_sum[i, "Poisson"] <- s_pois_da(models[[i]], obs)
  # we can do it like that, since LM is the first score we calculate in the loop
  t_sum[i, "Pois-diff"] <- t_sum[i, "Poisson"] - t_sum["LM", "Poisson"]
  t_sum[i, "IG"] <- sum(inf_gain(models[[which(model_names == "LM")]], models[[i]], obs))
  t_sum[i, "IGPE"] <- cum_igpe(models[[which(model_names == "LM")]], models[[i]], obs)[n_days]
}
t_sum
write.csv(t_sum, file.path(fpath, "Tab2_scores_inf.csv"))

col_names <- list(
  c("", "$\\bar{\\myS}$", "$\\bar\\myS - \\bar\\myS^{(\\textrm{LM})}$", "IG", "IGPE")
)
print_tex_table(
  t_sum,
  col_names = col_names,
  digits = c(2, 3, 3, 3),
  sep_after = "LG",
  col_align = "c",
  file_path = file.path(fpath, "Tab2_Scores_Info.tex")
)

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
  diff_scores <- rowSums(scf(fcst2, y) - scf(fcst1, y))
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = 6, type = "covariance", plot = F)$acf
  dm_var <- auto_covs[1] + 2 * sum(auto_covs[-1])

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(n_days)
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var),
                    mean_diff = mean_diff_score))
}

filter_days <- rep(T, n_days)
# 1. DM with Poisson -----------------------------------------------------------
test_func <- dm_test
scf_func <- s_pois
divide <- n_days
name <- "dm_pois"
n_digits <- 2
model_order <- c(5, 1, 4, 2, 3)
# 2. DM with Quadratic ---------------------------------------------------------
test_func <- dm_test
scf_func <- s_quad
divide <- n_days
name <- "dm_quad"
n_digits <- 3
model_order <- c(1, 5, 4, 2, 3)
# 3. CSEP T-test ---------------------------------------------------------------
test_func <- csep_test
scf_func <- s_pois
divide <- sum(obs)
name <- "csep"
n_digits <- 3
model_order <- c(5, 1, 4, 2, 3)
# filter for Mondays -------------
filt_day <- "Mo"
filter_days <- (lubridate::wday(times, label = T) == filt_day)
name <- paste0("csep_", filt_day)
divide <- sum(obs[filter_days, ])
# ------------------------------------------------------------------------------

sum_scores <- sapply(models, function(x) sum(scf_func(x[filter_days,], obs[filter_days,])))
val_matrix <- matrix(0.0, nrow = length(models), ncol = length(models),
                     dimnames = list(model_names[model_order], model_names[model_order]))

for (i in 1:(length(model_order) - 1)) {
  # diagonal corresponds to Poisson score
  val_matrix[i, i] <- sum_scores[model_order[i]] / divide

  for (j in (i + 1):length(model_order)) {
    test_vals <- test_func(models[[model_order[i]]][filter_days,], models[[model_order[j]]][filter_days,], obs[filter_days,],
                           scf_func)

    val_matrix[j, i] <- test_vals$pval
    val_matrix[i, j] <- test_vals$zval
  }
}
val_matrix[length(models), length(models)] <- sum_scores[model_order[length(models)]] / divide

write.csv(val_matrix, file.path(fpath, paste0("tests_", name, ".csv")))

print_tex_table(
  val_matrix,
  col_names = colnames(val_matrix),
  digits = rep(n_digits, length(models)),
  make_strong = rownames(val_matrix),   # just emphasize diagonals
  strong_symbol = "\\fbox",
  file_path = file.path(fpath, paste0("tests_", name, ".tex"))
)

# look at auto-covariance function --------------------------------------------
analyze_autocorrelation_fcn <- function(s_func = s_pois) {
  acf_data <- data.frame()

  for (i in 1:(length(model_order) - 1)) {
    for (j in (i + 1):length(model_order)) {
      score_diffs <- s_func(models[[model_order[i]]], obs) - s_func(models[[model_order[j]]], obs)

      acf_obj <- acf(rowSums(score_diffs), type = "covariance", plot = F)

      acf_data <- rbind(
        acf_data, data.frame(lag = acf_obj$lag, value = acf_obj$acf,
                             cmp = paste(model_names[model_order[i]], "vs.", model_names[model_order[j]]))
      )
    }
  }
  return(acf_data)
}

acf_data <- analyze_autocorrelation_fcn()

acf_plot <- ggplot(acf_data) +
  facet_wrap(~cmp, scales = "free_y", ncol = 2) +
  geom_hline(aes(yintercept = 0), size = 0.2, color = "gray") +
  geom_segment(aes(x = lag, y = value, xend = lag, yend = 0), size = 0.5) +
  xlab("Lag") +
  ylab("Auto-covariance") +
  my_theme

finish <- grid.arrange(acf_plot,
                       top = textGrob("ACF of Poisson score Difference",
                                      gp = gpar(fontsize = title_size)))

file_path <- file.path(fpath, "Fig_DM-acf-values.pdf")
ggsave(file_path, width = 140, height = 160, unit = "mm", plot = finish)

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
  scale_x_continuous(breaks = scores$X[new_year], labels = NULL, limits = c(0, nrow(scores))) +
  scale_color_manual(name = "", values = model_colors) +
  scale_shape_manual(name = "Number of\nearthquakes", labels = c("= 0", "> 0"),
                     values = c("FALSE" = 16, "TRUE" = 17),
                     guide = guide_legend(override.aes = list(alpha = 1, size = 1),
                                          direction = "vertical", title.position = "left")) +
  guides(color = "none") +
  scale_y_log10() +
  xlab(NULL) +
  ylab("Poisson score") +
  ggtitle(NULL) +
  my_theme +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
        legend.margin = margin(0, 0.0, 0, 3.5),
        plot.margin = margin(5.5, 5.5, 0, 5.5), axis.ticks.x = element_blank(),
        legend.background = element_blank())

# score differences ------------------------------------------------------------
diff_scores <- scores %>%
  mutate(across(all_of(ana_models), function(v) !!cmp_m - v)) %>%
  select(X, earthquake, all_of(ana_models)) %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model")

score_diff_plot <- ggplot() +
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
  scale_x_continuous(breaks = scores$X[new_year], labels = NULL, limits = c(0, nrow(scores))) +
  scale_y_continuous(trans = my_trans, breaks = my_breaks, labels = my_labels,
                     minor_breaks = minor_breaks, limits = my_limits) +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors, breaks = ana_models) +
  scale_shape_manual(name = "earthquake", labels = c("No", "At least one"),
                     values = c("FALSE" = 16, "TRUE" = 17)) +
  xlab(NULL) +
  ylab("Score difference") +
  ggtitle(NULL) +
  guides(color = "none", shape = "none") +
  my_theme +
  theme(plot.margin = margin(-3, 5.5, 0, 5.5), axis.ticks.x = element_blank())

# 1 cumulative information gain ------------------------------------------------
df_cum <- do.call(
  cbind, lapply(ana_models, function(i) cum_inf_gain(models[[which(model_names == "LM")]],
                                                     models[[which(model_names == i)]], obs))
) %>%
  as.data.frame() %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(df_cum) <- c(ana_models, "X", "earthquake")

cum_inf_plot <- df_cum %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = (1:n_days)[new_year], labels = NULL, limits = c(1, n_days)) +
  scale_color_manual(name = NULL, values = model_colors[model_names[model_order]],
                     guide = guide_legend(nrow = 1, keywidth = unit(4, "mm"))) +
  xlab(NULL) +
  ylab("Information Gain (IG)") +
  ggtitle(NULL) +
  my_theme +
  theme(plot.margin = margin(-3, 5.5, 0, 5.5), axis.ticks.x = element_blank(),
        legend.position = c(0.01, 0.99), legend.justification = c(0, 1))

# 2 cumulative information gain per event --------------------------------------
df_cum_pe <- do.call(
  cbind, lapply(ana_models, function(i) cum_igpe(models[[which(model_names == "LM")]],
                                                 models[[which(model_names == i)]], obs))
) %>%
  as.data.frame() %>%
  mutate(X = 1:nrow(.), earthquake = rowSums(obs) != 0)
colnames(df_cum_pe) <- c(ana_models, "X", "earthquake")

cum_inf_pe_plot <- df_cum_pe %>%
  pivot_longer(cols = all_of(ana_models), names_to = "Model") %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  geom_line(aes(x = X, y = value, color = Model), size = 0.3) +
  scale_x_continuous(breaks = (1:n_days)[new_year], labels = year(times[new_year]),
                     limits = c(1, n_days)) +
  scale_color_manual(name = paste(cmp_model, "vs."), values = model_colors, breaks = ana_models,
                     guide = guide_legend(order = 1, direction = "horizontal",
                                          override.aes = list(alpha = 0.75))) +
  xlab(NULL) +
  ylab("IG per event") +
  ggtitle(NULL) +
  guides(color = "none") +
  my_theme +
  theme(plot.margin = margin(-3, 5.5, 5.5, 5.5),
        legend.background = element_blank())

aligned <- align_plots(score_plot, score_diff_plot, cum_inf_plot, cum_inf_pe_plot,
                       align = "v")
# first values of each model of cum inf gain per event are NA since division by zero
combine_plots <- grid.arrange(aligned[[1]], aligned[[2]], aligned[[3]], aligned[[4]],
                              ncol = 1)
file_path <- file.path(fpath, "Fig4_score-and-info.pdf")
ggsave(file_path, width = 140, height = 170, unit = "mm", plot = combine_plots)

rm(scores, scores_long, diff_scores, combine_plots, score_plot, score_diff_plot, my_trans,
   df_cum, df_cum_pe, cum_inf_plot, cum_inf_pe_plot, aligned)

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

file_path <- file.path(fpath, paste0("Fig10_ScoreDiffSpat", add_name, ".pdf"))
ggsave(file_path, width = 145, height = 160, unit = "mm", plot = spat_plot)

rm(scores, diff_scores, combine_plots, spat_plot, eq_loc)


################################################################################
# Murphy diagramm
################################################################################

S_elem <- compiler::cmpfun(s_theta_da) # compile function to reduce runtime a bit

daily <- FALSE

n_theta <- 100
if (!daily) {
  log_grid <- seq(-24, 4, len = n_theta)
} else {
  log_grid <- seq(-7, 4, len = n_theta)
}
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
                              ymin = ifelse(daily, 0.24, 0.32),
                              ymax = ifelse(daily, 0.245, 0.325))

murphy_diag <- data.frame(murphy_df) %>%
  mutate(theta = log_grid) %>%
  pivot_longer(cols = all_of(model_names), names_to = "Model") %>%
  ggplot() +
  geom_line(aes(x = theta, y = value, color = Model), size = 0.3) +
  # geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "dashed") +
  geom_rect(data = best_models_bar, aes(color = Model, fill = Model, xmin = xmin,
                                        xmax = xmax, ymin = ymin, ymax = ymax),
            show.legend = F) +
  scale_color_manual(name = NULL, values = model_colors[model_names[model_order]]) +
  scale_fill_manual(name = NULL, values = model_colors) +
  xlab(expression(paste("log", (theta)))) +
  ylab("Elementary score") +
  my_theme +
  theme(legend.position = c(0.01, 0.1), legend.justification = c(0, 0))

file_path <- file.path(fpath, "Fig3_LogMurphyDiag.pdf")
ggsave(file_path, width = 145, height = 75, unit = "mm", plot = murphy_diag)

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
  cbind(data.frame(MCB_diag), Type = "MCB"),
  cbind(data.frame(DSC_diag), Type = "DSC")
)
colnames(df_collect) <- c(model_names, "Type")

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

if (!daily) {
  mcb_fac <- 5
  x_breaks <- -4:1 * 5
  y_breaks <- 0:3 * 0.1
  y_minor_breaks <- -4:6 * 0.05
  sec_breaks <- c(-2:0 / mcb_fac * 0.1, 1:3 / 2 * 1 / mcb_fac * 0.1)
  sec_labels <- c(2:0 * 0.1, 1:3 / 2 * 1 / mcb_fac)
  lower_lim <- NA
} else {
  mcb_fac <- 2
  x_breaks <- -3:2 * 2
  y_breaks <- 0:2 * 0.1
  y_minor_breaks <- -3:5 * 0.05
  sec_breaks <- c(-1:0 / mcb_fac * 0.1, 0:2 / 2 * 1 / mcb_fac * 0.1)
  sec_labels <- c(1:0 * 0.1, 0:2 / 2 * 1 / mcb_fac * 0.1)
  lower_lim <- -0.13
}
df_plot <- rbind(murphy_df, df_collect) %>%
  mutate(value = ifelse(Type == "Discrimination", value * (-1), value),
         value = ifelse(Type == "Miscalibration", value * mcb_fac, value))

murphy_score_cmps <- ggplot(df_plot) +
  geom_rect(data = best_models_bar, aes(color = Model, fill = Model, xmin = xmin,
                                        xmax = xmax, ymin = ymin, ymax = ymax),
            show.legend = F) +
  geom_hline(yintercept = 0.0, color = "black", size = 0.3) +
  geom_line(aes(x = log_theta, y = value, color = Model, group = paste(Model, Type)),
            size = 0.3) +
  scale_color_manual(name = NULL, values = model_colors[model_names[model_order]]) +
  scale_fill_manual(name = NULL, values = model_colors) +
  scale_x_continuous(breaks = x_breaks) +
  scale_y_continuous(breaks = y_breaks, minor_breaks = y_minor_breaks,
                     limits = c(lower_lim, NA),
                     sec.axis = sec_axis(~./mcb_fac, name = NULL,
                                         breaks = sec_breaks, labels = sec_labels)) +
  xlab(expression(paste("log", (theta)))) +
  ylab(NULL) +
  my_theme +
  annotate("text", x = -Inf, y = 0.12, label = "Average score", angle = 90, vjust = 2) +
  annotate("text", x = Inf, y = ifelse(daily, -0.05, -0.12), label = "DSC", angle = -90,
           vjust = 2) +
  annotate("text", x = Inf, y = ifelse(daily, 0.06,0.09), label = "MCB", angle = -90, vjust = 2) +
  theme(legend.position = c(0.01, 0.01), legend.justification = c(0, 0))

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
# Score component plot
################################################################################

dm_test <- function(fcst1, fcst2, y, scf, daily = F) {
  if (daily) {
    diff_scores <- scf(fcst2, y) - scf(fcst1, y)
  } else {
    diff_scores <- rowSums(scf(fcst2, y) - scf(fcst1, y))
  }
  mean_diff_score <- mean(diff_scores)

  auto_covs <- acf(diff_scores, lag.max = 6, type = "covariance", plot = F)$acf
  dm_var <- auto_covs[1] + 2 * sum(auto_covs[-1])

  test_stat <- mean_diff_score / sqrt(dm_var) * sqrt(n_days)
  pval <- 1 - pnorm(test_stat)

  return(data.frame(zval = test_stat, pval = pval, sd = sqrt(dm_var),
                    mean_diff = mean_diff_score))
}

get_score_cmp_plot <- function(results, non_sig_seg, daily = F) {
  scores_wide <- results %>%
    filter(is.finite(value)) %>%
    pivot_wider(id_cols = c(Model, Scoring), names_from = Type, values_from = value) %>%
    mutate(Scoring = ifelse(Scoring == "pois", "Poisson", "Quadratic"))
  unc <- mean(with(scores_wide, Score - MCB + DSC), na.rm = T)    # recover uncertainty component
  pois_panel <- scores_wide$Scoring[1] == "Poisson"

  fmt <- paste0("%.",  ifelse(pois_panel, 2, 3), "f")
  labelline_just <- ifelse(pois_panel, 0.48, 0.32)  # position of scores at lines

  # shift labels away from points to make them readable
  if (!daily) {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 1.2, 0.0, 1.2, -0.2),
                        vjusts = c(1.5, 0.0, -1.0, 1.5, -0.2))
  } else {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 0.5, 1.2, 0.5, 0.0),
                        vjusts = c(1.5, 1.5, 1.5, -1.0, 1.5))
  }

  scf <- ifelse(pois_panel, s_pois_da, s_quad_da)
  cmp_score <- ifelse(daily, scf(rowSums(models[[1]]), rowSums(obs)),
                      scf(models[[1]], obs))
  iso <- data.frame(
    intercept = seq(min(scores_wide$DSC) - max(scores_wide$MCB, na.rm = T),
                    max(scores_wide$DSC, na.rm = T),
                    length.out = 10)) %>%
    mutate(score = unc - intercept,       # miscalibration is 0, intercept corresponds to DSC
           # irrespective of daily or not sum obs is the sum over all earthquakes
           igpe = n_days / sum(obs) * (score - cmp_score),
           label = paste(sprintf(fmt, score), "/", sprintf(fmt, igpe)))

  # mcb_range <- range(scores_wide$MCB)
  dsc_range <- range(scores_wide$DSC)

  pl <- left_join(scores_wide, justs, by = "Model") %>%
    ggplot() +
    facet_wrap(~Scoring) +
    geom_abline(data = iso, aes(intercept = intercept, slope = 1.0), color = "lightgray",
                alpha = 0.5, size = 0.5) +
    geom_labelabline(data = iso, aes(intercept = intercept, slope = 1.0, label = label),
                     color = "gray50", hjust = labelline_just, size = 7 * 0.36, text_only = TRUE,
                     boxcolour = NA, straight = TRUE) +
    geom_segment(data = non_sig_seg, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = "dotted", alpha = 0.5, size = 0.5) +
    geom_point(aes(x = MCB, y = DSC, color = Model), size = 1.5) +
    geom_text(aes(x = MCB, y = DSC, label = Model, hjust = hjusts, vjust = vjusts, color = Model),
              size = 8 * 0.36) +
    scale_color_manual(values = model_colors) +
    coord_cartesian(ylim = (c(dsc_range[1] - 0.1 * (dsc_range[2] - dsc_range[1]), NA))) +
    annotate("label", x = Inf, y = -Inf, label = paste0("UNC = ", sprintf(fmt, unc)),
             hjust = 1.05, vjust = -0.2) +
    my_theme +
    theme(aspect.ratio = 1, legend.position = "none", axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  return(pl)
}

get_score_cmp_plot_2 <- function(results, non_sig_seg, daily = F) {
  scores_wide <- results %>%
    filter(is.finite(value)) %>%
    pivot_wider(id_cols = c(Model, Scoring), names_from = Type, values_from = value) %>%
    mutate(Scoring = ifelse(Scoring == "pois", "Poisson", "Quadratic"))
  unc <- mean(with(scores_wide, Score - MCB + DSC), na.rm = T)    # recover uncertainty component
  pois_panel <- scores_wide$Scoring[1] == "Poisson"

  fmt <- paste0("%.",  ifelse(pois_panel, 2, 3), "f")

  # shift labels away from points to make them readable
  if (!daily) {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 1.2, 0.0, 1.2, -0.2),
                        vjusts = c(1.5, 0.0, -1.0, 1.5, -0.2))
  } else {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 0.5, 1.5, 0.5, 0.0),
                        vjusts = c(1.5, 1.5, 1.5, -1.0, 1.5))
  }

  scf <- ifelse(pois_panel, s_pois_da, s_quad_da)
  cmp_score <- ifelse(daily, scf(rowSums(models[[1]]), rowSums(obs)), scf(models[[1]], obs))

  # define plot ranges manually, so that we exactly know where isolines start and end
  mcb_range <- range(scores_wide$MCB)
  dsc_range <- range(scores_wide$DSC)
  plot_max_mcb <- mcb_range[2] + 0.05 * (mcb_range[2] - mcb_range[1])
  plot_min_mcb <- mcb_range[1] - 0.05 * (mcb_range[2]- mcb_range[1])
  plot_max_dsc <- dsc_range[2] + 0.1 * (dsc_range[2] - dsc_range[1])
  plot_min_dsc <- dsc_range[1] - 0.2 * (dsc_range[2] - dsc_range[1])
  # isolines have slope 1, we can determine everything with linear equations :)
  iso <- data.frame(intercept = seq(plot_min_dsc - plot_max_mcb, plot_max_dsc - plot_min_mcb, length.out = 9)) %>%
    mutate(x_axis_l = intercept + plot_min_mcb < plot_min_dsc,
           x_axis_u = intercept + plot_max_mcb > plot_max_dsc)
  iso <- iso[2:(nrow(iso) - 1), ]   # drop first and last row because they just correspond to the corners

  # scores = 0 - intercept + unc
  # igpe = T / N (scores - comp_score)
  xbreaks <- plot_min_dsc - iso$intercept[iso$x_axis_l]   # MCB values, solve lin equ
  xlabels <- sprintf(fmt, xbreaks - plot_min_dsc + unc)
  ybreaks <- iso$intercept[!iso$x_axis_l] + plot_min_mcb  # DSC values, solve lin equ
  ylabels <- sprintf(fmt, plot_min_mcb - ybreaks + unc)
  x2breaks <- plot_max_dsc - iso$intercept[iso$x_axis_u]   # MCB values, solve lin equ
  x2labels <- sprintf(fmt, n_days / sum(obs) * ((x2breaks - plot_max_dsc + unc) - cmp_score))
  y2breaks <- iso$intercept[!iso$x_axis_u] + plot_max_mcb  # DSC values, solve lin equ
  y2labels <- sprintf(fmt, n_days / sum(obs) * ((plot_max_mcb - y2breaks + unc) - cmp_score))

  col1 <- "darkblue"
  col2 <- "darkred"

  pl <- left_join(scores_wide, justs, by = "Model") %>%
    ggplot() +
    geom_abline(data = iso, aes(intercept = intercept, slope = 1.0), color = "lightgray",
                alpha = 0.5, size = 0.5) +
    geom_segment(data = non_sig_seg, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = "dotted", alpha = 0.5, size = 0.5) +
    geom_point(aes(x = MCB, y = DSC, color = Model), size = 1.5) +
    geom_text(aes(x = MCB, y = DSC, label = Model, hjust = hjusts, vjust = vjusts, color = Model),
              size = 8 * 0.36) +
    scale_color_manual(values = model_colors) +
    scale_x_continuous(limits = c(plot_min_mcb, plot_max_mcb), breaks = xbreaks, labels = xlabels,
                       expand = c(0, 0), name = NULL,
                       sec.axis = sec_axis(~., name = NULL, breaks = x2breaks, labels = x2labels)) +
    scale_y_continuous(limits = c(plot_min_dsc, plot_max_dsc), breaks = ybreaks, labels = ylabels,
                       expand = c(0, 0), name = NULL,
                       sec.axis = sec_axis(~., name = NULL, breaks = y2breaks, labels = y2labels)) +
    annotate("label", x = Inf, y = -Inf, label = paste0("UNC = ", sprintf(fmt, unc)),
             hjust = 1.05, vjust = -0.2) +
    my_theme +
    theme(aspect.ratio = 1, legend.position = "none",
          axis.text = element_text(size = 7),
          axis.text.x.bottom = element_text(color = col1), axis.text.y.left = element_text(color = col1),
          axis.text.y.right = element_text(color = col2), axis.text.x.top = element_text(color = col2),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(pl)
}

get_score_cmp_plot_3 <- function(results, non_sig_seg, daily = F) {
  scores_wide <- results %>%
    filter(is.finite(value)) %>%
    pivot_wider(id_cols = c(Model, Scoring), names_from = Type, values_from = value) %>%
    mutate(Scoring = ifelse(Scoring == "pois", "Poisson", "Quadratic"))
  unc <- mean(with(scores_wide, Score - MCB + DSC), na.rm = T)    # recover uncertainty component
  pois_panel <- scores_wide$Scoring[1] == "Poisson"

  fmt <- paste0("%.",  ifelse(pois_panel, 2, 3), "f")

  # shift labels away from points to make them readable
  if (!daily) {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 1.2, 0.0, 1.2, -0.2),
                        vjusts = c(1.5, 0.0, -1.0, 1.5, -0.2))
  } else {
    justs <- data.frame(Model = c("LM", "FCM", "LG", "SMA", "LRWA"),
                        hjusts = c(0.0, 0.5, 1.5, 0.5, 0.0),
                        vjusts = c(1.5, 1.5, 1.5, -1.0, 1.5))
  }

  scf <- ifelse(pois_panel, s_pois_da, s_quad_da)
  # cmp_score <- ifelse(daily, scf(rowSums(models[[1]]), rowSums(obs)), scf(models[[1]], obs))
  cmp_score <- 2.68

  # define plot ranges manually, so that we exactly know where isolines start and end
  mcb_range <- range(scores_wide$MCB)
  dsc_range <- range(scores_wide$DSC)
  plot_max_mcb <- mcb_range[2] + 0.2 * (mcb_range[2] - mcb_range[1])
  plot_min_mcb <- mcb_range[1] - 0.2 * (mcb_range[2]- mcb_range[1])
  plot_max_dsc <- dsc_range[2] + 0.2 * (dsc_range[2] - dsc_range[1])
  plot_min_dsc <- dsc_range[1] - 0.2 * (dsc_range[2] - dsc_range[1])
  # isolines have slope 1, we can determine everything with linear equations :)
  iso <- data.frame(intercept = seq(plot_min_dsc - plot_max_mcb, plot_max_dsc - plot_min_mcb, length.out = 9)) %>%
    mutate(x_axis_l = intercept + plot_min_mcb < plot_min_dsc,
           x_axis_u = intercept + plot_max_mcb > plot_max_dsc)
  iso <- iso[2:(nrow(iso) - 1), ]   # drop first and last row because they just correspond to the corners

  line_angle <- atan((plot_max_mcb - plot_min_mcb) / (plot_max_dsc - plot_min_dsc)) * 180 / pi

  # scores = mcb - dsc + unc
  # igpe = T / N (scores - comp_score)
  x_shift <- ifelse(daily, 0.15, 1.7)
  y_shift <- ifelse(daily, 1.7, 0.15)
  iso_annotate <- rbind(
    # lower x axis
    data.frame(x = plot_min_dsc - iso$intercept[iso$x_axis_l], y = plot_min_dsc, hjust = -x_shift) %>%
               mutate(text = sprintf(fmt, x - plot_min_dsc + unc)),
    # left y axis
    data.frame(x = plot_min_mcb, y = iso$intercept[!iso$x_axis_l] + plot_min_mcb, hjust = -y_shift) %>%
               mutate(text = sprintf(fmt, plot_min_mcb - y + unc)),
    # upper x axis
    data.frame(x = plot_max_dsc - iso$intercept[iso$x_axis_u], y = plot_max_dsc, hjust = 1 + x_shift) %>%
               mutate(text = sprintf(fmt, n_days / sum(obs) * ((x - plot_max_dsc + unc) - cmp_score))),
    # right y axis
    data.frame(x = plot_max_mcb, y = iso$intercept[!iso$x_axis_u] + plot_max_mcb, hjust = 1 + y_shift) %>%
               mutate(text = sprintf(fmt, n_days / sum(obs) * ((plot_max_mcb - y + unc) - cmp_score)))
  )

  area_width_x <- ifelse(daily, 0.002, 0.01)
  area_width_y <- ifelse(daily, 0.0011, 0.025)
  iso_txt_bg <- data.frame(
    xmin = c(plot_min_mcb, plot_max_mcb - area_width_x, plot_min_mcb, plot_min_mcb),
    xmax = c(plot_min_mcb + area_width_x, plot_max_mcb, plot_max_mcb, plot_max_mcb),
    ymin = c(plot_min_dsc, plot_min_dsc, plot_min_dsc, plot_max_dsc - area_width_y),
    ymax = c(plot_max_dsc, plot_max_dsc, plot_min_dsc + area_width_y, plot_max_dsc)
  )

  pl <- left_join(scores_wide, justs, by = "Model") %>%
    ggplot() +
    geom_abline(data = iso, aes(intercept = intercept, slope = 1.0), color = "lightgray",
                alpha = 0.5, size = 0.5) +
    # draw white area where line annotations will be
    geom_rect(data = iso_txt_bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white") +
    # now draw the text
    geom_text(data = iso_annotate, aes(x = x, y = y, label = text, hjust = hjust), angle = line_angle, size = 7 / .pt,
               color = "lightgray") +
    geom_segment(data = non_sig_seg, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = "dotted", alpha = 0.5, size = 0.5) +
    geom_point(aes(x = MCB, y = DSC, color = Model), size = 1.5) +
    geom_text(aes(x = MCB, y = DSC, label = Model, hjust = hjusts, vjust = vjusts, color = Model),
              size = 8 * 0.36) +
    scale_color_manual(values = model_colors) +
    scale_x_continuous(limits = c(plot_min_mcb, plot_max_mcb), expand = c(0, 0), name = NULL) +
    scale_y_continuous(limits = c(plot_min_dsc, plot_max_dsc), expand = c(0, 0), name = NULL) +
    # annotate("label", x = Inf, y = -Inf, label = paste0("UNC = ", sprintf(fmt, unc)),
    #          hjust = 1.05, vjust = -0.2) +
    my_theme +
    theme(aspect.ratio = 1, legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(pl)
}

get_score_cmps <- function(x, y, name, scf_list = list("pois" = s_pois_da, "quad" = s_quad_da), daily = F) {
  if (!daily) {
    x <- as.vector(x)
    y <- as.vector(y)
  } else {
    x <- rowSums(x)
    y <- rowSums(y)
  }

  scores <- data.frame()
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]
  x_rc <- monotone(y)

  for (scf_name in names(scf_list)) {
    scf <- scf_list[[scf_name]]
    s <- scf(x, y)
    s_rc <- scf(x_rc, y)
    s_mg <- scf(mean(y), y)

    scores <- rbind(
      scores,
      data.frame(Model = name, Scoring = scf_name, Type = c("Score", "MCB", "DSC"),
                 value = c(s, s - s_rc, s_mg - s_rc))
    )
  }
  return(scores)
}

get_non_sig_connections <- function(df_score_cmp, scf, level = 0.1, daily = F) {
  non_sig_con <- data.frame()

  for (i in 1:(length(models) - 1)) {
    for (j in ((i + 1):length(models))) {
      p <- dm_test(models[[i]], models[[j]], obs, scf, daily = daily)$pval
      if ((p > level / 2) & (p < 1 - level / 2)) {
        non_sig_con <- rbind(
          non_sig_con,
          data.frame(x = filter(df_score_cmp, Model == model_names[i], Type == "MCB")$value,
                     y = filter(df_score_cmp, Model == model_names[i], Type == "DSC")$value,
                     xend = filter(df_score_cmp, Model == model_names[j], Type == "MCB")$value,
                     yend = filter(df_score_cmp, Model == model_names[j], Type == "DSC")$value)
        )
      }
    }
  }
  return(non_sig_con)
}

daily <- T
df_score_cmp <- do.call(
  rbind,
  lapply(1:length(models), function(i) get_score_cmps(models[[i]], obs, model_names[i], daily = daily))
)
non_sig_seg_pois <- get_non_sig_connections(filter(df_score_cmp, Scoring == "pois"), s_pois)

df_score_cmp <- read.csv("./figures8/score-cmps.csv")
non_sig_seg_pois <- read.csv("./figures8/score-cmps_seg-pois.csv")
pl_pois <- get_score_cmp_plot_3(filter(df_score_cmp, Scoring == "pois"), non_sig_seg_pois)

non_sig_seg_pois <- read.csv("./figures9/df_score-cmps_seg-pois-daily.csv")
pl_pois_daily <- get_score_cmp_plot_3(filter(df_score_cmp, Scoring == "pois"), non_sig_seg_pois, daily = T)

my_plot <- grid.arrange(pl_pois, pl_pois_daily, nrow = 1,
                        bottom = textGrob("MCB", gp = gpar(fontsize = 11)),
                        left = textGrob("DSC", rot = 90, gp = gpar(fontsize = 11)))

file_path <- file.path(fpath, "Fig9_MCB-DSC-plot-2.pdf")
ggsave(file_path, width = 145, height = 80, unit = "mm", plot = my_plot)

rm(my_plot, pl_pois, pl_pois_daily, df_score_cmp)

################################################################################
# Reliability diagramm
################################################################################

reldiag <- function(x, y, n_resamples = 99, region_level = 0.9) {
  ord <- order(x, y, decreasing = c(FALSE, TRUE))
  x <- x[ord]
  y <- y[ord]

  # compress fit by only storing values around jumps that belong to distinct x-values
  # (if there is a horizontal section of width 0 with more than two data points,
  #  we would have two identical points in the data though we need only one)
  filter_jumps <- function(x, v) {
    n <- length(v)
    jumps <- (v[-1] - v[-n] > 0)
    # filter points around jumps (points immediately before and after a jump)
    next_to_jump <- c(T, jumps) | c(jumps, T)
    x_at_jump <- x[next_to_jump]
    # set to False if corresponding x-values do not change
    next_to_jump[next_to_jump] <- c(T, x_at_jump[-1] - x_at_jump[-length(x_at_jump)] > 0)
    return(next_to_jump)
  }

  score <- s_pois_da

  x_rc <- monotone(y)
  s <- score(x,y)
  s_rc <- score(x_rc, y)
  s_mg <- score(mean(y), y)

  mcb <- s - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg

  jumps <- filter_jumps(x, x_rc)
  rc_fit <- data.table(x = x[jumps], y = x_rc[jumps])

  t <- table(y) / length(y)
  vals <- as.numeric(names(t))
  f_y <- as.numeric(t)
  eps <- sum(f_y[-1] * vals[-1]) / x - 1
  n_split <- 5      # as not everything fits into memory, split
  split <- ceiling(seq(1, length(y), length.out = n_split))
  cmp <- matrix(rep(cumsum(f_y[-length(f_y)]), split[2] - split[1] + 1),
                ncol = length(f_y) - 1, byrow = T)

  collect_vals <- list()
  for (i in 1:n_resamples) {
    for (j in 2:n_split) {
      i_vec <- split[j - 1]:split[j]
      u <- runif(length(i_vec), 0, 1 + eps[i_vec])
      y[i_vec] <- rowSums(u - eps[i_vec] >= cmp[1:pmin(length(i_vec), nrow(cmp)), ])
    }
    ord <- order(x, y, decreasing = c(FALSE, TRUE))
    x_rc <- monotone(y[ord])
    jumps <- filter_jumps(x, x_rc)
    # x is already sorted, so x[ord] would just resort x on sections on which it
    # is already constant
    collect_vals[[i]] <- data.table(x = x[jumps], y = x_rc[jumps], I = paste0("R", i))
  }
  x_grid <- sort(unique(c(do.call(c, lapply(collect_vals, function(vals) vals$x)), rc_fit$x)))
  # linearly interpolate resampled fits on x values of original fit
  resampl_fits <- do.call(cbind, lapply(collect_vals, function(df) {
    approx(x = df$x, y = df$y, xout = x_grid, method = "linear", yleft = NA, yright = NA)$y
  }))
  val_low_and_up <- apply(resampl_fits, 1, function(row) {
      sorted <- sort(row)
      low <- floor(length(sorted) * (1 - region_level) / 2)   # maybe there are NA, so recalculate it
      up <- length(sorted) - low
      low <- pmax(low, 1)
      c(sorted[low], sorted[up])
    }) %>% t() %>%  # apply writes results in columns
    as.data.table()

  results <- cbind(
    data.table(x = x_grid,
               x_rc = approx(x = rc_fit$x, y = rc_fit$y, xout = x_grid, method = "linear", yleft = NA, yright = NA)$y),
    val_low_and_up
  )
  colnames(results) <- c("x", "x_rc", "lower", "upper")
  stats <- data.frame(Score = s, MCB = mcb, DSC = dsc, UNC = unc)
  return(list(results = results, stats = stats))
}

reldiag_cmp <- compiler::cmpfun(reldiag)

recal_models <- data.table()
collect_stats <- data.table()

set.seed(999)

for (i in 1:length(models)) {
  print(paste(model_names[i], Sys.time()))
  res <- reldiag_cmp(as.vector(models[[i]]), as.vector(obs), n_resamples = 100)
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
recal_models <- read.csv("./figures9/df_rel-Til100.csv", row.names = 1) %>%
  filter(Model %in% model_names)
collect_stats <- read.csv("./figures9/df_rel-stats.csv", row.names = 1) %>%
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
    arrange(a) %>%      # restore original sorting
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

# for daily ----------------------------------------
breaks <- c(0, 0.2, 0.3, 0.4, 1)
t_breaks <- my_trans(breaks)
labels <- c("0", rep("", length(breaks) - 2), "1")
add_name <- "_daily-ecdf"

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
xmin <- 0.2
xmax <- 1.0
ymin <- 0.05
ymax <- 0.7
# axis breaks and labels
t_breaks <- seq(0, 2, length.out = 6)
labels <- c("0", rep("", length(breaks) - 2), 2)
# position of score components
text_x <- -1
text_y <- 0.6
# add to name
add_name <- "_std"
# cut confidence band at plot max
recal_models$upper <- pmin(recal_models$upper, 2.2)

# For daily forecasts ----------------------------------------------------------
# 1 - log transform ------------------------------------------------------------
my_trans <- function(x) log(x, base = 10)
# position of inset histograms
xmin <- 0.2
xmax <- 1.2
ymin <- -1.4
ymax <- -0.6
# axis breaks and labels
breaks <- 10^(-3:1)
t_breaks <- my_trans(breaks)
labels <- breaks
# position of score components
text_x <- -1.3
text_y <- 0.22
# add to name
add_name <- "_daily-log"
# we need to take care about plot min
vec <- c(recal_models$x, recal_models$lower)
plot_min <- my_trans(min(vec[vec > 0]))

# ------------------------------------------------------------------------------

show_ticks <- element_blank()
show_ticks <- element_line(colour = "black", size = 0.3)

inset_histograms <- list()
plot_min <- my_trans(min(c(recal_models$x, recal_models$lower)))
plot_max <- my_trans(max(c(recal_models$x, recal_models$upper)))
hist_breaks <- seq(plot_min, plot_max, length.out = 9)
for (i in 1:length(models)) {
  my_hist <- ggplot(data.table(x = as.vector(models[[i]]))) +
    # geom_vline(xintercept = t_breaks, size = 0.2) +
    geom_histogram(aes(x = my_trans(x)), fill = "gray", col = "black", size = 0.2,
                   breaks = hist_breaks) +
    scale_x_continuous(breaks = t_breaks) +
    theme_classic(base_size = 5.5) +
    theme(axis.line.y = element_blank(), axis.text = element_blank(),
          axis.ticks.x = show_ticks, axis.ticks.y = element_blank(),
          axis.ticks.length = unit(1, "mm"),
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

create_row <- function(row, top = TRUE) {
  dt_recal <- filter(recal_models, Model %in% row)
  dt_stats <- filter(collect_stats, Model %in% row)

  top_margin <- ifelse(top, 5.5, 0.0)
  bot_margin <- ifelse(top, 0.0, 5.5)

  df_segments <-  dt_recal %>%
    group_by(Model) %>%
    group_modify(function(df, key) {
      jumps <- (df$x_rc[-1] - df$x_rc[-nrow(df)] > 0)
      df %>%
        mutate(jump_before = c(T, jumps), jump_after = c(jumps, T),
               segment_end_point = xor(jump_before, jump_after)) %>%
        filter(segment_end_point) %>%
        mutate(x_pos = ifelse(jump_before, "x0", "x1")) %>%
        pivot_wider(id_cols = "x_rc", names_from = "x_pos", values_from = "x")
    }, .keep = T)

  return(
    ggplot(dt_recal, aes(x = my_trans(x))) +
      facet_wrap(~factor(Model, ordered = T, levels = row), nrow = 1) +
      geom_ribbon(aes(ymin = my_trans(lower), ymax = my_trans(upper), fill = Model),
                  alpha = 0.33, show.legend = FALSE) +
      geom_abline(intercept = 0 , slope = 1, colour = "grey70", size = 0.3,
                  linetype = "dashed") +
      geom_line(aes(y = my_trans(x_rc), color = Model), size = 0.3, show.legend = FALSE) +
      geom_segment(data = df_segments,
                   aes(x = my_trans(x0), xend = my_trans(x1), y = my_trans(x_rc),
                       yend = my_trans(x_rc), color = Model), size = 0.7, show.legend = FALSE) +
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
            panel.grid.minor = element_blank(),
            plot.margin = margin(top_margin, 3.5, bot_margin, 3.5))
  )
}

combine <- grid.arrange(create_row(rows[[1]], top = T) +
                          inset_histograms[which(model_names %in% rows[[1]])],
                        create_row(rows[[2]], top = F) +
                          inset_histograms[which(model_names %in% rows[[2]])],
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
