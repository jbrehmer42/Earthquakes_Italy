## Checking for mean calibration 

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


#################################################
#### Part I - Murphy diagram for MCB and DSC ####

## Compute values for MCB and DSC diagrams
n_theta <- 100
log_grid <- seq(-25, 1, len = n_theta)
grd <- exp(log_grid)
MCB_diag <- DSC_diag <- matrix(0, ncol = n_mods, nrow = n_theta)
for (i in 1:n_mods) {
  decomp <- bin_decomposition(models[[i]], obs, theta = grd)
  MCB_diag[ ,i] <- rowSums( decomp$MCB )
  DSC_diag[ ,i] <- rowSums( decomp$DSC )
}
UNC_diag <- rowSums( decomp$UNC )
rm(decomp)

#save(MCB_diag, DSC_diag, UNC_diag, file = paste0(path, '/MCB_etc.RData'))

## Plot MCB diagram
file_path <- file.path(fpath, "plot_MCB_diag.pdf")
plot_elementary_scores(MCB_diag, grd, model_names, model_colors, file_path,
                       "MCB")

## Plot DSC diagram
file_path <- file.path(fpath, "plot_DSC_diag.pdf")
plot_elementary_scores(DSC_diag, grd, model_names, model_colors, file_path,
                       "DSC")

## Plot UNC diagram
file_path <- file.path(fpath, "plot_UNC_diag.pdf")
pdf(file_path, width = 8, height = 5.5)
par(mar = c(4, 4, 0.5, 0.5))
plot(1:n_theta, UNC_diag, ty = "l", xlab = "log(theta)", ylab = "UNC",
     xaxt = "n") 
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(log_grid[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
dev.off()


##################################################
#### Part II - Mean reliability diagrams etc. ####

# Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)

## Compute values for MCB and DSC diagrams (AGGREGATED!)
n_theta <- 100
log_grid <- seq(-6, 2, len = n_theta)
grd <- exp(log_grid)
MCB_diag_agg <- DSC_diag_agg <- matrix(0, ncol = n_mods, nrow = n_theta)
for (i in 1:n_mods) {
  decomp <- bin_decomposition(models_agg[[i]], obs_agg, theta = grd)
  MCB_diag_agg[ ,i] <- decomp$MCB
  DSC_diag_agg[ ,i] <- decomp$DSC
}
UNC_diag_agg <- decomp$UNC


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag.pdf")
lim1 <- c(0, 6)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T)
}
dev.off()


## Plot MCB diagram
file_path <- file.path(fpath, "plot_MCB_diag_agg.pdf")
plot_elementary_scores(MCB_diag_agg, grd, model_names, model_colors, file_path,
                       "MCB")

## Plot DSC diagram
file_path <- file.path(fpath, "plot_DSC_diag_agg.pdf")
plot_elementary_scores(DSC_diag_agg, grd, model_names, model_colors, file_path,
                       "DSC")

## Plot UNC diagram
file_path <- file.path(fpath, "plot_UNC_diag_agg.pdf")
pdf(file_path, width = 8, height = 6)
plot(1:n_theta, UNC_diag_agg, ty = "l", xlab = "log(theta)", ylab = "UNC",
     xaxt = "n")
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(log_grid[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
dev.off()

