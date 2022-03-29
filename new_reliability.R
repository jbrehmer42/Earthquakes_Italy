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
Mlimit <- 5
TI_large <- events$TI[events$MAG > Mlimit]
TI_large <- unique(TI_large)
subset_index <- rep(TRUE, ndays)
subset_index[c(TI_large, TI_large + 1)] <- FALSE
ndays - sum(subset_index)

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# Erase days with large events
#bs_agg <- obs_agg[subset_index]
#for (i in 1:nmods) models_agg[[i]] <- models_agg[[i]][subset_index]


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
filePath <- file.path(fpath, "mean_reldiag_ex1.pdf")
lim1 <- c(0, 6)
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
filePath <- file.path(fpath, "mean_reldiag_log.pdf")
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
filePath <- file.path(fpath, "mean_reldiag_ex2_m5_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim2, ln = T, resamp = T)
}
dev.off()
