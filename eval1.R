## Model evaluation as in CSEP

# Do Data preparation
source('~/Documents/Code/Earthquakes_Italy/data_prep.R')
# source functions for scores
source('~/Documents/Code/Earthquakes_Italy/functions_eval.R')
# source functions for plotting
source('~/Documents/Code/Earthquakes_Italy/functions_plot.R')


# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots_TEST"

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

### --- TEST ---------------------------------->

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
# Interchange rows and columns for Mphy_diag array

## Plot of Murphy diagram
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag.pdf"
plotElementary(Mphy_diag, grd, mnames, mcols, filePath, "score")


##############################################
#### Part II - Fully aggregated forecasts ####

# Aggregate mean forecasts
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
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
# Include climatology: Have to switch
# to reduced testing region
subs <- region_intersect(clima, bins)

## subset forecast models and observations
clima <- clima[subs$clima, ]
obs <- obs[ ,subs$model]
for (i in 1:nmods) models[[i]] models[[i]][ ,subs$model]

# Aggregate mean forecasts again, as we
# now consider a smaller testing region
models_agg <- list()
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# use climatology as last model
models_agg[[nmods+1]] <- rep(sum(clima$RATE), ndays)

# PAV-transformed forecasts and climatology
# attention: this changes the ordering of the forecasts
models_pav <- list()
for (i in 1:(nmods+1)) models_pav <- isoreg(models_agg[[i]], obs_agg)$yf
# model5_pav <- isoreg(model5_agg, obs_agg)$yf
# model1_pav <- isoreg(model1_agg, obs_agg)$yf
# model2_pav <- isoreg(model2_agg, obs_agg)$yf
# model3_pav <- isoreg(model3_agg, obs_agg)$yf
# model4_pav <- isoreg(model4_agg, obs_agg)$yf

# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
ntheta <- 100
lgrd <- seq(-6, 2, len = ntheta)
grd <- exp(lgrd)

Mphy_diag_agg <- Mphy_diag_pav <- matrix(0, ncol = ntheta, nrow = nmods+1)
for (i in 1:(nmods+1)) {
  Mphy_diag_agg[i, ] <- colMeans( Sthet_vec(models_agg[[i]], obs_agg, grd) )
  Mphy_diag_pav[i, ] <- colMeans( Sthet_vec(models_pav[[i]], obs_agg[order(models_agg[[i]])], grd) )
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
mcols2 <- c(mcols, "magenta")
mnames2 <- c(mnames, "Clima")
mnames3 <- c(mnames2, paste0("RC_", mnames2))

## Murphy diagram 1 (only aggregated forecasts)
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag_agg.pdf"
plotElementary(Mphy_diag_agg, grd, mnames2, mcols2, filePath, "score", whichmods = 1:5)


## Murphy diagram 2 (aggregated and recalibrated forecasts)
filePath <- "~/Documents/_temp/Case/Plots/plot_Murphy_diag_agg_pav.pdf"
plotElementary(rbind(Mphy_diag_agg, Mphy_diag_pav), grd, mnames3, c(mcols2, mcols2), filePath,
               "score", mltys = rep(1:2, ea = 5), whichmods = 1:10)

# Remove after testing
# pdf(filePath, width = 8, height=6)
# plot(1:ntheta, Mphy_diag_agg[1, ], ty = "l", xlab = "log(theta)",
#      ylab = "score", main = "Murphy diagram (fully aggregated)", xaxt = "n", ylim = c(0,0.18))
# # get different x-axis
# for (i in 1:5) {
#   lines(1:ntheta, Mphy_diag_pav[i, ], col = cols2[i], lty = 2)
#   if (i == 1) next
#   lines(1:ntheta, Mphy_diag_agg[i, ], col = cols2[i])
# }
# # create log axis
# ticks <- axis(1, labels = F, tick = F)
# labs <- round(lgrd[pmax(1,ticks)], 1)
# axis(1, at = ticks, labels = labs)
# legend(3, 0.17, mnames2, col = cols2, lwd = 2, title = "Original")
# legend(3, 0.095, mnames2, col = cols2, lwd = 2, lty = 2, title = "Recalibrated")
# dev.off()


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

## TO DO:
## Use plotScores function for this plot

