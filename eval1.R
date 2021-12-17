## Model evaluation as in CSEP

# Do Data preparation
source('~/Documents/_temp/Case/data_prep.R')
# source functions for scores
source('~/Documents/_temp/Case/functions_eval.R')


############################################
#### Part I - Scores and Murphy diagram ####


## Calculate daily scores for quadratic and Poisson loss
scores_pois <- scores_quad <- matrix(0, nrow = ndays, ncol = 4)
for (i in 1:4) {
  mname <- paste0("model", i)
  scores_pois[ ,i] <- rowSums( Spois(get(mname), obs) )
  scores_quad[ ,i] <- rowSums( Squad(get(mname), obs) )
}

save(scores_pois, scores_quad, file = paste0(path, "/scores.RData"))

## Mean scores for full testing period
colMeans(scores_pois)
colMeans(scores_quad)


## Score for multiple bins of Taggart's decomposed
## scoring functions for the mean
binlim <- c(0, exp(-14), exp(-12), exp(-11.2), exp(-10), 20)
scores_bin <- list()

for (i in 1:(length(binlim)-1)) {
  a <- binlim[i]
  b <- binlim[i+1]
  vals <- matrix(0, nrow = ndays, ncol = 4)
  for (k in 1:4) {
    mname <- paste0("model", k)
    cat("\n")
    cat(paste(mname, "bin", i, "of", length(binlim)-1))
    for (j in 1:ndays) {
      vals[j,k] <- Sbin( get(mname)[j, ], obs[j, ], a, b)
    }
  }
  scores_bin[[i]] <- vals
}

save(scores_bin, file = paste0(path, '/scores_bin.RData'))


## Murphy diagrams
# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
ntheta <- 100
lgrd <- seq(-25, 1, len = ntheta)
grd <- exp(lgrd)

tt <- Sys.time()
Mphy_list <- list()
for (j in 1:ntheta) {
  vals <- matrix(0, nrow = ndays, ncol = 4)
  for (i in 1:4) {
    mname <- paste0("model", i)
    for (k in 1:ndays) {
      vals[k,i] <- Sthet(get(mname)[k, ], obs[k, ], grd[j])
    }
  }
  Mphy_list[[j]] <- vals
  cat("\n")
  cat(paste0(j, "/", ntheta))
}
tt <- Sys.time() - tt

save(Mphy_list, file = paste0(path, '/Mphy_list.RData'))

Mphy_diag <- t( sapply(Mphy_list, "colMeans") )



## Time dependent plots of scores
# load(paste(path, '/scores.RData', sep=''))
# prepare time axis
tindex <- (times2$DD == 1) & (times2$MM == 1)
labs <- times2$YY[tindex]
evt_days <- unique(M4events$TI)

## logarithmic
filePath <- "~/Documents/_temp/Case/Plots/plot_pois_time.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, scores_pois[ ,1], ty = "l", ylim = c(-1,70), xlab = "days",
     ylab = "score", main = "Logarithmic score", xaxt = "n")
for (i in 2:4) {
  lines(1:ndays, scores_pois[ ,i], col = cols[i])
}
points(evt_days, rep(-1, length(evt_days)), cex = 0.8)
axis(1, at = which(tindex), labels = labs)
legend(4500, 58, mnames, col = cols, lwd = 2)
dev.off()

## logarithmic (log scale)
filePath <- "~/Documents/_temp/Case/Plots/plot_pois_time_log.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, log(scores_pois[ ,1]), ty = "l", ylim = c(-2,5.4), xlab = "days",
     ylab = "score (log scale)", main = "Logarithmic score (log scale)", xaxt="n")
for (i in 2:4) {
  lines(1:ndays, log(scores_pois[ ,i]), col = cols[i])
}
points(evt_days, rep(-2, length(evt_days)), cex = 0.8)
axis(1, at = which(tindex), labels = labs)
legend(0, 5.3, mnames, col = cols, lwd = 2)
dev.off()

## quadratic
filePath <- "~/Documents/_temp/Case/Plots/plot_quad_time.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, scores_quad[ ,1], ty = "l", ylim = c(0, 50), xlab = "days",
     ylab = "score", main = "Quadratic score", xaxt = "n")
for (i in 2:4) {
  lines(1:ndays, scores_quad[ ,i], col = cols[i])
}
axis(1, at = which(tindex), labels = labs)
legend(4500, 48, mnames, col = cols, lwd = 2)
dev.off()

## quadratic (log scale)
filePath <- "~/Documents/_temp/Case/Plots/plot_quad_time_log.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, log(scores_quad[ ,1]), ty = "l", ylim = c(-11,5.2), xlab = "days",
     ylab = "score (log scale)", main = "Quadratic score (log scale)", xaxt = "n")
for (i in 2:4) {
  lines(1:ndays, log(scores_quad[ ,i]), col = cols[i])
}
axis(1, at = which(tindex), labels = labs)
legend(0, 5.1, mnames, col = cols, lwd = 2)
dev.off()

## bin (lowest)
filePath <- "~/Documents/_temp/Case/Plots/plot_bin_time.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, scores_bin[[1]][ ,1], ty = "l", ylim = c(4.7e-9, 6.2e-9), xlab = "days",
     ylab = "score", main = "Lowest bin score", xaxt = "n")
for (i in 2:4) {
  lines(1:ndays, scores_bin[[1]][ ,i], col = cols[i])
}
axis(1, at = which(tindex), labels = labs)
legend(0, 5.6e-9, mnames, col = cols, lwd = 2)
dev.off()

## Plot of Murphy diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, Mphy_diag[ ,1], ty = "l", ylim = c(0,0.16), xlab = "log(theta)",
     ylab = "score", main = "Murphy diagram", xaxt = "n")
# get different x-axis
for (i in 2:4) {
  lines(1:ntheta, Mphy_diag[ ,i], col = cols[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(2, 0.157, mnames, col = cols, lwd = 2)
dev.off()


##############################################
#### Part II - Fully aggregated forecasts ####

# Aggregate mean forecasts
model1_agg <- rowSums(model1)
model2_agg <- rowSums(model2)
model3_agg <- rowSums(model3)
model4_agg <- rowSums(model4)
obs_agg <- rowSums(obs)


## Calculate daily scores for quadratic and Poisson loss
scores_pois_agg <- scores_quad_agg <- matrix(0, nrow = ndays, ncol = 4)
for (i in 1:4) {
  mname <- paste0("model", i, "_agg")
  scores_pois_agg[ ,i] <- Spois(get(mname), obs_agg) 
  scores_quad_agg[ ,i] <- Squad(get(mname), obs_agg) 
}

colMeans(scores_pois_agg)
colMeans(scores_quad_agg)

save(scores_pois_agg, scores_quad_agg, file = paste0(path, "/scores_agg.RData"))

## Murphy diagrams

# Include climatology: Have to switch
# to reduced testing region
source('~/Documents/_temp/Case/clima_subset.R')

# Aggregate mean forecasts again, as we
# now consider a smaller testing region
model1_agg <- rowSums(model1)
model2_agg <- rowSums(model2)
model3_agg <- rowSums(model3)
model4_agg <- rowSums(model4)
obs_agg <- rowSums(obs)

# PAV-transformed forecasts and climatology
model5_agg <- rep(sum(clima$RATE), ndays)
# attention: this changes the ordering of the forecasts
model5_pav <- isoreg(model5_agg, obs_agg)$yf
model1_pav <- isoreg(model1_agg, obs_agg)$yf
model2_pav <- isoreg(model2_agg, obs_agg)$yf
model3_pav <- isoreg(model3_agg, obs_agg)$yf
model4_pav <- isoreg(model4_agg, obs_agg)$yf

# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
ntheta <- 100
lgrd <- seq(-6, 2, len = ntheta)
grd <- exp(lgrd)

Mphy_diag_agg <- Mphy_diag_pav <- matrix(0, ncol = ntheta, nrow = 5)
for (i in 1:5) {
  mname1 <- paste0("model", i, "_agg")
  mname2 <- paste0("model", i, "_pav")
  Mphy_diag_agg[i, ] <- colMeans( Sthet_vec(get(mname1), obs_agg, grd) )
  Mphy_diag_pav[i, ] <- colMeans( Sthet_vec(get(mname2), obs_agg[order(get(mname1))], grd) )
}


## Time dependent plots of scores

## logarithmic
filePath <- "~/Documents/_temp/Case/Plots/plot_pois_time_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, scores_pois_agg[ ,1], ty = "l", ylim = c(-14, 30), xlab = "days",
     ylab = "score", main = "Logarithmic score (fully aggregated)")
for (i in 2:4) {
  lines(1:ndays, scores_pois_agg[ ,i], col = cols[i])
}
legend(4500, 29, mnames, col = cols, lwd = 2)
dev.off()

## logarithmic (log scale)
filePath <- "~/Documents/_temp/Case/Plots/plot_pois_time_agg_log.pdf"
cutoff <- 0.01
pdf(filePath, width = 8, height=6)
plot(1:ndays, log(pmax(scores_pois_agg[ ,1], cutoff)), ylim = c(log(cutoff)+0.1, 3.5), ty = "l", xlab = "days",
     ylab = "score (log scale)", main = "Logarithmic score (fully aggregated, log scale)")
for (i in 2:4) {
  lines(1:ndays, log(pmax(scores_pois_agg[ ,i], cutoff)), col = cols[i])
}
legend(-85, 3.5, mnames, col = cols, lwd = 2)
dev.off()

## quadratic
filePath <- "~/Documents/_temp/Case/Plots/plot_quad_time_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, scores_quad_agg[ ,1], ty = "l", ylim = c(0,400), xlab = "days",
     ylab = "score", main = "Quadratic score (fully aggregated)")
for (i in 2:4) {
  lines(1:ndays, scores_quad_agg[ ,i], col = cols[i])
}
legend(4600, 398, mnames, col = cols, lwd = 2)
dev.off()

## quadratic (log scale)
filePath <- "~/Documents/_temp/Case/Plots/plot_quad_time_agg_log.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, log(scores_quad_agg[ ,1]), ty = "l", ylim = c(-10.5, 7), xlab = "days",
     ylab = "score (log scale)", main = "Quadratic score (fully aggregated, log scale)")
for (i in 2:4) {
  lines(1:ndays, log(scores_quad_agg[ ,i]), col = cols[i])
}
legend(0, 6.8, mnames, col = cols, lwd = 2)
dev.off()

## Plot of Murphy diagrams
cols2 <- c(cols, "magenta")
mnames2 <- c(mnames, "Clima")

## Murphy diagram 1 (only aggregated forecasts)
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag_agg.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, Mphy_diag_agg[1, ], ty = "l", xlab = "log(theta)",
     ylab = "score", main = "Murphy diagram (fully aggregated)", xaxt = "n", ylim = c(0,0.18))
# get different x-axis
for (i in 2:5) {
  lines(1:ntheta, Mphy_diag_agg[i, ], col = cols2[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.17, mnames2, col = cols2, lwd = 2)
dev.off()

## Murphy diagram 2 (aggregated and recalibrated forecasts)
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag_agg_pav.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ntheta, Mphy_diag_agg[1, ], ty = "l", xlab = "log(theta)",
     ylab = "score", main = "Murphy diagram (fully aggregated)", xaxt = "n", ylim = c(0,0.18))
# get different x-axis
for (i in 1:5) {
  lines(1:ntheta, Mphy_diag_pav[i, ], col = cols2[i], lty = 2)
  if (i == 1) next
  lines(1:ntheta, Mphy_diag_agg[i, ], col = cols2[i])
}
# create log axis
ticks <- axis(1, labels = F, tick = F)
labs <- round(lgrd[pmax(1,ticks)], 1)
axis(1, at = ticks, labels = labs)
legend(3, 0.17, mnames2, col = cols2, lwd = 2, title = "Original")
legend(3, 0.095, mnames2, col = cols2, lwd = 2, lty = 2, title = "Recalibrated")
dev.off()

## Plot of mean forecasts over time
filePath <- "~/Documents/_temp/Case/Plots/plot_mean_forecasts.pdf"
pdf(filePath, width = 8, height=6)
plot(1:ndays, model1_agg, ty = "l", xlab = "days",
     ylab = "mean forecast", main = "Mean forecasts for each day (fully aggregated)")
# get different x-axis
for (i in 2:5) {
  mname <- paste0("model", i, "_agg")
  lines(1:ndays, get(mname), col = cols2[i])
}
legend(10, 4.1, mnames2, col = cols2, lwd = 2)
dev.off()

