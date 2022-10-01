## Checking for mean calibration
## Now excluding the large events

# Path for figures
fpath <- "/media/myData/Plots/"
# Path for r scripts
rpath <- "/media/myData/Doks/Forschung/Code/Earthquakes_Italy"

# Do Data preparation
source(file.path(rpath, "data_prep.R"))
# source functions for scores
source(file.path(rpath, "functions_eval.R"))
# source functions for plotting
source(file.path(rpath, "functions_plot.R"))


## Exclude large earthquakes (first variant)
Mlimit <- 5.5
time_index_large <- events$TI[events$MAG > Mlimit]
time_index_large <- unique(time_index_large)
subset_index <- rep(TRUE, n_days)
for (i in 1:length(time_index_large)) {
  subset_index[(time_index_large[i]-6):time_index_large[i]] <- FALSE
}
n_days - sum(subset_index)

## Exclude large earthquakes (second variant)
Mlimit <- 5.5
time_index_large <- events$TI[events$MAG > Mlimit]
time_index_large <- unique(time_index_large)
subset_index <- rep(TRUE, n_days)
subset_index[c(time_index_large+1, time_index_large + 2)] <- FALSE
n_days - sum(subset_index)

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# Erase days with large events
obs_agg <- obs_agg[subset_index]
for (i in 1:n_mods) models_agg[[i]] <- models_agg[[i]][subset_index]


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_ex1.pdf")
lim1 <- c(0, 4.5)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T, MCB_decomp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_ex1_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T, MCB_decomp = T)
}
dev.off()

# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_new.pdf")
lim1 <- c(0, 6)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T, MCB_decomp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_new_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T, MCB_decomp = T)
}
dev.off()

###############################################
## Plot reliability diagram for each weekday ##

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)

day_names <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
               "Saturday", "Sunday")
day_start <- c(3,4,5,6,7,1,2)

for (i in 1:7) {
  day_ind <- seq(day_start[i], n_days, by = 7)
  models_days <- lapply(models_agg, function(x) x[day_ind])
  obs_days <- obs_agg[day_ind]
  # Make plot
  # Reliability diagram on log scale
  file_path <- file.path(fpath, paste("mean_reldiag", day_names[i], "log.pdf",
                                      sep = "_"))
  lim2 <- c(0.08, 9)
  pdf(file_path, width = 7.5, height = 8)
  par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
  for (j in 1:n_mods) {
    plot_reliability(models_days[[j]], obs_days, model_names[j],
                     model_colors[j], lim = lim2, ln = T, resamp = T,
                     MCB_decomp = T)
  }
  dev.off()
}
