## Model evaluation as in CSEP

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


############################################
#### Part I - Scores and Murphy diagram ####

days0 <- days1 <- 1:ndays
days0[(rowSums(obs) != 0)] <- NA
days1[(rowSums(obs) == 0)] <- NA

## Calculate daily scores for quadratic and Poisson loss
scores_pois <- scores_quad <- scores_spat <- scores_mass <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois[ ,i] <- rowSums( Spois(models[[i]], obs) )
  scores_quad[ ,i] <- rowSums( Squad(models[[i]], obs) )
  scores_spat[ ,i] <- rowSums( Sspat(models[[i]], obs) ) + rowSums(obs) * log(rowSums(models[[i]]))
  # scores_mass considers the fully aggregated forecasts and observations (see below)
  scores_mass[ ,i] <- Spois( rowSums(models[[i]]), rowSums(obs) )
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


## Time dependent plots of scores and score differences

## Poisson Scores
filePath <- file.path(fpath, "plot_pois_time.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "l")

filePath <- file.path(fpath, "plot_pois_time_h.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "h")

filePath <- file.path(fpath, "plot_pois_time0.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- file.path(fpath, "plot_pois_time1.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- file.path(fpath, "plot_pois_time_zoom.pdf")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))


## Poisson differences (to model 1)
## Positive values indicate superior performance
scorediffs <- as.matrix(scores_pois[ ,rep(1, nmods-1)] - scores_pois[ ,2:nmods])

filePath <- file.path(fpath, "plot_pois_time_diff.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- file.path(fpath, "plot_pois_time_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1))

filePath <- file.path(fpath, "plot_pois_time_diff0.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1), days = days0)

filePath <- file.path(fpath, "plot_pois_time_diff1.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-20,20), days = days1)

filePath <- file.path(fpath, "plot_pois_time_diff_trim_zoom.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic Scores
filePath <- file.path(fpath, "plot_quad_time.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "l")

filePath <- file.path(fpath, "plot_quad_time_h.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "h")

filePath <- file.path(fpath, "plot_quad_time0.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- file.path(fpath, "plot_quad_time1.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- file.path(fpath, "plot_quad_time_zoom.pdf")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic differences (to model 1)
## Positive values indicate superior performance
scorediffs <- as.matrix(scores_quad[ ,rep(1, nmods-1)] - scores_quad[ ,2:nmods])

filePath <- file.path(fpath, "plot_quad_time_diff.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- file.path(fpath, "plot_quad_time_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-0.1,0.1))

filePath <- file.path(fpath, "plot_quad_time_diff0.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-0.1,0.1), days = days0)

filePath <- file.path(fpath, "plot_quad_time_diff1.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-1,0.7), days = days1)

filePath <- file.path(fpath, "plot_quad_time_diff_trim_zoom.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-3,0.5), tlim = list(c(1,8,2016), c(31,5,2018)))


## Spatial Poisson scores
filePath <- file.path(fpath, "plot_spat_time1.pdf")
plotScores(scores_spat, times, mnames, mcols, filePath, events, type = "p", days = days1)

## Spatial Poisson differences (to model 1)
scorediffs <- as.matrix(scores_spat[ ,rep(1, nmods-1)] - scores_spat[ ,2:nmods])

filePath <- file.path(fpath, "plot_spat_time_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               days = days1, trim = c(-20,10))

filePath <- file.path(fpath, "plot_spat_time_diff_trim_zoom.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               days = days1, trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))


## Mass Poisson scores
filePath <- file.path(fpath, "plot_mass_time.pdf")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "l")

filePath <- file.path(fpath, "plot_mass_time_h.pdf")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "h")

filePath <- file.path(fpath, "plot_mass_time0.pdf")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- file.path(fpath, "plot_mass_time1.pdf")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- file.path(fpath, "plot_mass_time_zoom.pdf")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))


## Mass Poisson differences (to model 1)
## Positive values indicate superior performance
scorediffs <- as.matrix(scores_mass[ ,rep(1, nmods-1)] - scores_mass[ ,2:nmods])

filePath <- file.path(fpath, "plot_mass_time_diff.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- file.path(fpath, "plot_mass_time_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1))

filePath <- file.path(fpath, "plot_mass_time_diff0.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1), days = days0)

filePath <- file.path(fpath, "plot_mass_time_diff1.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-20,20), days = days1)

filePath <- file.path(fpath, "plot_mass_time_diff_trim_zoom.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))



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


## Time dependent plots of scores and score differences
## (now fully aggregated)

## Poisson
filePath <- file.path(fpath, "plot_pois_time_agg.pdf")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events)

filePath <- file.path(fpath, "plot_pois_time_agg0.pdf")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- file.path(fpath, "plot_pois_time_agg1.pdf")
plotScores(scores_pois_agg, times, mnames, mcols, filePath, events, type = "p", days = days1)

## Poisson differences (to model 1)
## Positive values indicate superior performance
scorediffs <- as.matrix(scores_pois_agg[ ,rep(1, nmods-1)] - scores_pois_agg[ ,2:nmods])

filePath <- file.path(fpath, "plot_pois_time_agg_diff.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- file.path(fpath, "plot_pois_time_agg_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1))

filePath <- file.path(fpath, "plot_pois_time_agg_diff1.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-20,20), days = days1)

## quadratic
filePath <- file.path(fpath, "plot_quad_time_agg.pdf")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events, logscale = F) 

## quadratic (log scale)
filePath <- file.path(fpath, "plot_quad_time_agg_log.pdf")
plotScores(scores_quad_agg, times, mnames, mcols, filePath, events)

## Quadratic differences (to model 1)
## Positive values indicate superior performance
scorediffs <- as.matrix(scores_quad_agg[ ,rep(1, nmods-1)] - scores_quad_agg[ ,2:nmods])

filePath <- file.path(fpath, "plot_quad_time_agg_diff.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- file.path(fpath, "plot_quad_time_agg_diff_trim.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-0.1,0.1))

filePath <- file.path(fpath, "plot_quad_time_agg_diff1.pdf")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-1,0.7), days = days1)



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
## Use plotScores function although this means that ylab is wrong
mean_models <- matrix(unlist(models_agg), ncol = 5)

filePath <- file.path(fpath, "plot_mean_forecasts.pdf")
plotScores(mean_models, times, mnames2, mcols2, filePath, events, logscale = F)

filePath <- file.path(fpath, "plot_mean_forecasts_log.pdf")
plotScores(mean_models, times, mnames2, mcols2, filePath, events)

