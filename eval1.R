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

# save(Mphy_list, file = file.path(dpath, 'Mphy_list.RData'))

Mphy_diag <- t( sapply(Mphy_list, "colMeans") )



## Time dependent plots of scores

# Plot of daily scores (Poisson)
filePath <- file.path(fpath, "plot_pois_time.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events = events, logscale = F) 

# Plot of daily scores (Poisson, log)
filePath <- file.path(fpath, "plot_pois_time_log.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events = events) 

# Plot of daily scores (Quadratic)
filePath <- file.path(fpath, "plot_quad_time.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events = events, logscale = F) 

# Plot of daily scores (Quadratic)
filePath <- file.path(fpath, "plot_quad_time_log.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events = events) 


## Plot of Murphy diagram
filePath <- file.path(fpath, "plot_Murphy_diag.pdf")
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
subs <- region_intersect(clima, cells)

## subset forecast models and observations
clima <- clima[subs$clima, ]
obs <- obs[ ,subs$model]
for (i in 1:nmods) models[[i]] <- models[[i]][ ,subs$model]

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
for (i in 1:(nmods+1)) models_pav[[i]] <- isoreg(models_agg[[i]], obs_agg)$yf


# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
ntheta <- 100
lgrd <- seq(-6, 2, len = ntheta)
grd <- exp(lgrd)

Mphy_diag_agg <- Mphy_diag_pav <- matrix(0, nrow = ntheta, ncol = nmods+1)
for (i in 1:(nmods+1)) {
  Mphy_diag_agg[ ,i] <- colMeans( Sthet_vec(models_agg[[i]], obs_agg, grd) )
  Mphy_diag_pav[ ,i] <- colMeans( Sthet_vec(models_pav[[i]], obs_agg[order(models_agg[[i]])], grd) )
}


## Time dependent plots of scores

## Poisson
filePath <- file.path(fpath, "plot_pois_time_agg.pdf")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events = events, logscale = F) 

## Poisson (log scale)
filePath <- file.path(fpath, "plot_pois_time_agg_log.pdf")
## obsolete with next plotScores version
scores_pois_agg2 <- pmax(scores_pois_agg, -2) + 3
plotScores(scores_pois_agg2, times, mnames, mcols, filePath, events = events) 

## quadratic
filePath <- file.path(fpath, "plot_quad_time_agg.pdf")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events = events, logscale = F) 

## quadratic (log scale)
filePath <- file.path(fpath, "plot_quad_time_agg_log.pdf")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events = events) 

## Plot of Murphy diagrams
mcols2 <- c(mcols, "magenta")
mnames2 <- c(mnames, "Clima")
mnames3 <- c(mnames2, paste0("RC_", mnames2))

## Murphy diagram 1 (only aggregated forecasts)
filePath <- file.path(fpath, "plot_Murphy_diag_agg.pdf")
plotElementary(Mphy_diag_agg, grd, mnames2, mcols2, filePath, "score", whichmods = 1:5)


## Murphy diagram 2 (aggregated and recalibrated forecasts)
filePath <- file.path(fpath, "plot_Murphy_diag_agg_pav.pdf")
plotElementary(cbind(Mphy_diag_agg, Mphy_diag_pav), grd, mnames3, c(mcols2, mcols2), filePath,
               "score", mltys = rep(1:2, ea = 5), whichmods = 1:10)


## Plot of mean forecasts over time
mean_models <- matrix(unlist(models_agg), ncol = 5)

filePath <- file.path(fpath, "plot_mean_forecasts.pdf")
plotScores(mean_models, times, mnames2, mcols2, filePath, events = events, logscale = F)

filePath <- file.path(fpath, "plot_mean_forecasts_log.pdf")
plotScores(mean_models, times, mnames2, mcols2, filePath, events = events)

