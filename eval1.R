## Model evaluation as in CSEP

# Do Data preparation
source('~/Documents/Code/Earthquakes_Italy/data_prep.R')
# source functions for scores
source('~/Documents/Code/Earthquakes_Italy/functions_eval.R')
# source functions for plotting
source('~/Documents/Code/Earthquakes_Italy/functions_plot.R')


# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots"

############################################
#### Part I - Scores and Murphy diagram ####


## Calculate daily scores for quadratic and Poisson loss
scores_pois <- scores_quad <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois[ ,i] <- rowSums( Spois(models[[i]], obs) )
  scores_quad[ ,i] <- rowSums( Squad(models[[i]], obs) )
}

## Mean scores for full testing period
colMeans(scores_pois)
colMeans(scores_quad)


## Murphy diagrams
# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
ntheta <- 100
lgrd <- seq(-25, 1, len = ntheta)
grd <- exp(lgrd)

tt <- Sys.time()
Mphy_list <- list()
for (j in 1:ntheta) {
  vals <- matrix(0, nrow = ndays, ncol = nmods)
  for (i in 1:nmods) {
    for (k in 1:ndays) {
      vals[k,i] <- Sthet(models[[i]][k, ], obs[k, ], grd[j])
    }
  }
  Mphy_list[[j]] <- vals
  cat("\n")
  cat(paste0(j, "/", ntheta))
}
tt <- Sys.time() - tt

save(Mphy_list, file = paste0(dpath, '/Mphy_list.RData'))

Mphy_diag <- t( sapply(Mphy_list, "colMeans") )



## Time dependent plots of scores

# Plot of daily scores (Poisson)
filePath <- paste(fpath, "plot_pois_time.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events = events, logscale = F) 

# Plot of daily scores (Poisson, log)
filePath <- paste(fpath, "plot_pois_time_log.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events = events) 

# Plot of daily scores (Quadratic)
filePath <- paste(fpath, "plot_quad_time.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events = events, logscale = F) 

# Plot of daily scores (Quadratic)
filePath <- paste(fpath, "plot_quad_time_log.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events = events) 



## TO DO:
# Create function for Murphy diagrams

## Plot of Murphy diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag.pdf"
plotElementary(Mphy_diag, grd, mnames, mcols, filePath, "score")


plotElementary <- function(vals, grd, mnames, mcols, filePath, ylab, mltys = NULL,
                           whichmods = 1:4) {
  ntheta <- length(grd)
  lgrd <- log(grd)
  nmods <- length(whichmods)
  if (missing(mltys)) mltys <- rep(1, nmods)
  mnames <- mnames[whichmods]
  mcols <- mcols[whichmods]
  mltys <- mltys[whichmods]
  ylim <- c(min(vals), max(vals))
  # Start plotting
  pdf(filePath, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:ntheta, 1:ntheta, ylim = ylim, xlab = "log(theta)", ylab = ylab,
       xaxt = "n", col = "transparent")
  for (i in 1:nmods) {
    lines(1:ntheta, vals[ ,i], col = mcols[i], lty = mltys[i])
  }
  # create log axis
  ticks <- axis(1, labels = F, tick = F)
  labs <- round(lgrd[pmax(1,ticks)], 1)
  axis(1, at = ticks, labels = labs)
  legend(2, ylim[2], mnames, col = mcols, lty = mltys, lwd = 2)
  dev.off()
}


##############################################
#### Part II - Fully aggregated forecasts ####

# Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] rowSums(models[[i]])
obs_agg <- rowSums(obs)


## Calculate daily scores for quadratic and Poisson loss
scores_pois_agg <- scores_quad_agg <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois_agg[ ,i] <- Spois(models_agg[[i]], obs_agg) 
  scores_quad_agg[ ,i] <- Squad(models_agg[[i]], obs_agg) 
}

# colMeans(scores_pois_agg)
# colMeans(scores_quad_agg)

## Murphy diagrams


## TO DO: 
# This has to be changed to a function if the climate
# subset part is sorted out

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

## Poisson
filePath <- paste(fpath, "plot_pois_time_agg.pdf", sep = "/")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events = events, logscale = F) 

## Poisson (log scale)
filePath <- paste(fpath, "plot_pois_time_agg_log.pdf", sep = "/")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events = events) 

## quadratic
filePath <- paste(fpath, "plot_quad_time_agg.pdf", sep = "/")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events = events, logscale = F) 

## quadratic (log scale)
filePath <- paste(fpath, "plot_quad_time_agg_log.pdf", sep = "/")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events = events) 

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

