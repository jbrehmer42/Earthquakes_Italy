## Checking for mean calibration 

# Do Data preparation
source('~/Documents/_temp/Case/first_steps.R')
# source functions for scores
source('~/Documents/_temp/Case/eval_functions.R')
source('~/Documents/_temp/Case/plot_functions.R')

#################################################
#### Part I - Murphy diagram for MCB and DSC ####

## Compute values for MCB and DSC diagrams
ntheta <- 100
lgrd <- seq(-25, 1, len = ntheta)
grd <- exp(lgrd)
MCB_diag <- DSC_diag <- matrix(0, ncol = 4, nrow = ntheta)
for (i in 1:4) {
  mname <- paste0("model", i)
  decomp <- bin_decomp(get(mname), obs, theta = grd)
  MCB_diag[ ,i] <- rowSums( decomp$MCB )
  DSC_diag[ ,i] <- rowSums( decomp$DSC )
}
UNC_diag <- rowSums( decomp$UNC )
rm(decomp)

save(MCB_diag, DSC_diag, UNC_diag, file = paste0(path, '/MCB_etc.RData'))

## Plot MCB diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_MCB_diag.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, MCB_diag[ ,1], ty = "l", xlab = "log(theta)",
     ylab = "MCB", main = "Miscalibration diagram", xaxt = "n", ylim = c(0,0.104))
# get different x-axis
for (i in 2:4) {
  lines(1:ntheta, MCB_diag[ ,i], col = cols[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.1, mnames, col = cols, lwd = 2)
dev.off()

## Plot DSC diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_DSC_diag.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, DSC_diag[ ,1], ty = "l", xlab = "log(theta)",
     ylab = "DSC", main = "Discrimination diagram", xaxt = "n", ylim = c(0,0.057))
# get different x-axis
for (i in 2:4) {
  lines(1:ntheta, DSC_diag[ ,i], col = cols[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.055, mnames, col = cols, lwd = 2)
dev.off()

## Plot UNC diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_UNC_diag.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, UNC_diag, ty = "l", xlab = "log(theta)",
     ylab = "UNC", main = "Uncertainty diagram", xaxt = "n") 
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
dev.off()


##################################################
#### Part II - Mean reliability diagrams etc. ####

# Aggregate mean forecasts
model1_agg <- rowSums(model1)
model2_agg <- rowSums(model2)
model3_agg <- rowSums(model3)
model4_agg <- rowSums(model4)
obs_agg <- rowSums(obs)

## Compute values for MCB and DSC diagrams (AGGREGATED!)
ntheta <- 100
lgrd <- seq(-6, 2, len = ntheta)
grd <- exp(lgrd)
MCB_diag_agg <- DSC_diag_agg <- matrix(0, ncol = 4, nrow = ntheta)
for (i in 1:4) {
  mname <- paste0("model", i, "_agg")
  decomp <- bin_decomp(get(mname), obs_agg, theta = grd)
  MCB_diag_agg[ ,i] <- decomp$MCB
  DSC_diag_agg[ ,i] <- decomp$DSC
}
UNC_diag_agg <- decomp$UNC


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)
filePath <- "~/Documents/_temp/Case/Plots/mean_reldiag.pdf"
lim1 <- c(0, 6)
pdf(filePath, width = 7.5, height=8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
simple_diag(model1_agg, obs_agg, mnames[1], col = cols[1], lim = lim1, resamp = T)
simple_diag(model2_agg, obs_agg, mnames[2], col = cols[2], lim = lim1, resamp = T)
simple_diag(model3_agg, obs_agg, mnames[3], col = cols[3], lim = lim1, resamp = T)
simple_diag(model4_agg, obs_agg, mnames[4], col = cols[4], lim = lim1, resamp = T)
dev.off()

filePath <- "~/Documents/_temp/Case/Plots/mean_reldiag_log.pdf"
lim2 <- c(0.05, lim1[2])
pdf(filePath, width = 7.5, height=8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
simple_diag(model1_agg, obs_agg, mnames[1], cols[1], lim = lim2, ln = T, resamp = T)
simple_diag(model2_agg, obs_agg, mnames[2], cols[2], lim = lim2, ln = T, resamp = T)
simple_diag(model3_agg, obs_agg, mnames[3], cols[3], lim = lim2, ln = T, resamp = T)
simple_diag(model4_agg, obs_agg, mnames[4], cols[4], lim = lim2, ln = T, resamp = T)
dev.off()



## Plot MCB diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_MCB_diag_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, MCB_diag_agg[ ,1], ty = "l", xlab = "log(theta)",
     ylab = "MCB", main = "Miscalibration diagram (fully aggregated)", xaxt = "n", ylim = c(0,0.02))
# get different x-axis
for (i in 2:4) {
  lines(1:ntheta, MCB_diag_agg[ ,i], col = cols[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.019, mnames, col = cols, lwd = 2)
dev.off()

## Plot DSC diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_DSC_diag_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, DSC_diag_agg[ ,1], ty = "l", xlab = "log(theta)",
     ylab = "DSC", main = "Discrimination diagram (fully aggregated)", xaxt = "n", ylim = c(0,0.04))
# get different x-axis
for (i in 2:4) {
  lines(1:ntheta, DSC_diag_agg[ ,i], col = cols[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.039, mnames, col = cols, lwd = 2)
dev.off()

## Plot UNC diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_UNC_diag_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, UNC_diag_agg, ty = "l", xlab = "log(theta)",
     ylab = "UNC", main = "Uncertainty diagram (fully aggregated)", xaxt = "n")
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
dev.off()

