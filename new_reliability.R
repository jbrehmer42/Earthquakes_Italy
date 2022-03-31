## Checking for mean calibration
## excluding the large events

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots"
# Path for R code
rpath <- "/home/jbrehmer/Documents/Code/Earthquakes_Italy"

# Do Data preparation
source(file.path(rpath, "data_prep.R"))
# source functions for scores
source(file.path(rpath, "functions_eval.R"))
# source functions for plotting
source(file.path(rpath, "functions_plot.R"))


## Exclude large earthquakes (first variant)
Mlimit <- 5.5
TI_large <- events$TI[events$MAG > Mlimit]
TI_large <- unique(TI_large)
subset_index <- rep(TRUE, ndays)
for (i in 1:length(TI_large)) {
  subset_index[(TI_large[i]-6):TI_large[i]] <- FALSE
}
ndays - sum(subset_index)

## Exclude large earthquakes (second variant)
Mlimit <- 5.5
TI_large <- events$TI[events$MAG > Mlimit]
TI_large <- unique(TI_large)
subset_index <- rep(TRUE, ndays)
subset_index[c(TI_large+1, TI_large + 2)] <- FALSE
ndays - sum(subset_index)

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# Erase days with large events
obs_agg <- obs_agg[subset_index]
for (i in 1:nmods) models_agg[[i]] <- models_agg[[i]][subset_index]


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
filePath <- file.path(fpath, "mean_reldiag_ex1.pdf")
lim1 <- c(0, 4.5)
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
filePath <- file.path(fpath, "mean_reldiag_ex1_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim2, ln = T, resamp = T)
}
dev.off()

# Reliability diagram
filePath <- file.path(fpath, "mean_reldiag_ex2.pdf")
lim1 <- c(0, 4.5)
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
filePath <- file.path(fpath, "mean_reldiag_ex2_m55_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim2, ln = T, resamp = T)
}
dev.off()

###############################################
## Plot reliability diagram for each weekday ##

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)

## Divide by weekday
# for (i in 1:7) {
#   dstring <- paste(times$YY[i], times$MM[i], times$DD[i], sep = "-")
#   print(weekdays(as.Date(dstring)))
# }
day_names <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
day_start <- c(3,4,5,6,7,1,2)

for (i in 1:7) {
  day_ind <- seq(day_start[i], ndays, by = 7)
  models_days <- lapply(models_agg, function(x) x[day_ind])
  obs_days <- obs_agg[day_ind]
  # Make plot
  # Reliability diagram on log scale
  filePath <- file.path(fpath, paste("mean_reldiag", day_names[i], "log.pdf", sep = "_"))
  lim2 <- c(0.08, 9)
  pdf(filePath, width = 7.5, height = 8)
  par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
  for (j in 1:nmods) {
    plotReliability(models_days[[j]], obs_days, mnames[j], mcols[j], lim = lim2, ln = T, resamp = T)
  }
  dev.off()
}
