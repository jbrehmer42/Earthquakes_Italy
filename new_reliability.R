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


#################################################
#### Part I - Murphy diagram for MCB and DSC ####

## Compute values for MCB and DSC diagrams
ntheta <- 100
lgrd <- seq(-25, 1, len = ntheta)
grd <- exp(lgrd)
MCB_diag <- DSC_diag <- matrix(0, ncol = nmods, nrow = ntheta)
for (i in 1:nmods) {
  decomp <- bin_decomp(models[[i]], obs, theta = grd)
  MCB_diag[ ,i] <- rowSums( decomp$MCB )
  DSC_diag[ ,i] <- rowSums( decomp$DSC )
}
UNC_diag <- rowSums( decomp$UNC )
rm(decomp)

#save(MCB_diag, DSC_diag, UNC_diag, file = paste0(path, '/MCB_etc.RData'))

## Plot MCB diagram
filePath <- file.path(fpath, "plot_MCB_diag.pdf")
plotElementary(MCB_diag, grd, mnames, mcols, filePath, "MCB")

## Plot DSC diagram
filePath <- file.path(fpath, "plot_DSC_diag.pdf")
plotElementary(DSC_diag, grd, mnames, mcols, filePath, "DSC")

## Plot UNC diagram
filePath <- file.path(fpath, "plot_UNC_diag.pdf")
pdf(filePath, width = 8, height = 5.5)
par(mar = c(4, 4, 0.5, 0.5))
plot(1:ntheta, UNC_diag, ty = "l", xlab = "log(theta)", ylab = "UNC", xaxt = "n") 
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
dev.off()


##################################################
#### Part II - Mean reliability diagrams etc. ####

# Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)

## Compute values for MCB and DSC diagrams (AGGREGATED!)
ntheta <- 100
lgrd <- seq(-6, 2, len = ntheta)
grd <- exp(lgrd)
MCB_diag_agg <- DSC_diag_agg <- matrix(0, ncol = nmods, nrow = ntheta)
for (i in 1:nmods) {
  decomp <- bin_decomp(models_agg[[i]], obs_agg, theta = grd)
  MCB_diag_agg[ ,i] <- decomp$MCB
  DSC_diag_agg[ ,i] <- decomp$DSC
}
UNC_diag_agg <- decomp$UNC


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
filePath <- file.path(fpath, "mean_reldiag.pdf")
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