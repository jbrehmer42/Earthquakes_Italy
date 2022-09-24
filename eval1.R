## Model evaluation as in CSEP

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


############################################
#### Part I - Scores and Murphy diagram ####

days0 <- days1 <- 1:n_days
days0[(rowSums(obs) != 0)] <- NA
days1[(rowSums(obs) == 0)] <- NA

## Calculate daily scores for quadratic and Poisson loss
scores_pois <- scores_quad <- matrix(0, nrow = n_days, ncol = n_mods)
scores_spat <- scores_mass <- matrix(0, nrow = n_days, ncol = n_mods)
for (i in 1:n_mods) {
  scores_pois[ ,i] <- rowSums( S_pois(models[[i]], obs) )
  scores_quad[ ,i] <- rowSums( S_quad(models[[i]], obs) )
  scores_spat[ ,i] <- rowSums( S_spat(models[[i]], obs) ) +
    rowSums(obs) * log(rowSums(models[[i]]))
  # scores_mass considers the fully aggregated forecasts
  # and observations (see below)
  scores_mass[ ,i] <- S_pois( rowSums(models[[i]]), rowSums(obs) )
}

## Mean scores for full testing period
colMeans(scores_pois)
colMeans(scores_quad)

## Murphy diagrams
# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
n_theta <- 100
log_grid <- seq(-25, 1, len = n_theta)
grd <- exp(log_grid)

Murphy_list <- list()
for (j in 1:n_theta) {
  vals <- matrix(0, nrow = n_days, ncol = n_mods)
  for (i in 1:n_mods) {
    for (k in 1:n_days) {
      vals[k,i] <- S_theta(models[[i]][k, ], obs[k, ], grd[j])
    }
  }
  Murphy_list[[j]] <- vals
  cat("\n")
  cat(paste0(j, "/", n_theta))
}

# save(Murphy_list, file = file.path(dpath, 'Murphy_list.RData'))
Murphy_diag <- t( sapply(Murphy_list, "colMeans") )


## Time dependent plots of scores and score differences

## Poisson Scores
file_path <- file.path(fpath, "plot_pois_time.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events,
            type = "l")

file_path <- file.path(fpath, "plot_pois_time_h.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events,
            type = "h")

file_path <- file.path(fpath, "plot_pois_time0.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events,
            type = "l", days = days0)

file_path <- file.path(fpath, "plot_pois_time1.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events,
            type = "p", days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_pois_time_zoom.pdf")
plot_scores(scores_pois, times, model_names, model_colors, file_path, events,
            type = "p", tlim = list(c(1,8,2016), c(31,5,2018)))


## Poisson differences (to model 1)
## Positive values indicate superior performance
score_diffs <- as.matrix(scores_pois[ ,rep(1, n_mods-1)]
                        - scores_pois[ ,2:n_mods])

file_path <- file.path(fpath, "plot_pois_time_diff.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4)

file_path <- file.path(fpath, "plot_pois_time_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-2,1))

file_path <- file.path(fpath, "plot_pois_time_diff0.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-2,1),
                 days = days0)

file_path <- file.path(fpath, "plot_pois_time_diff1.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4,  trim = c(-20,20),
                 days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_pois_time_diff_trim_zoom.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-20,10),
                 tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic Scores
file_path <- file.path(fpath, "plot_quad_time.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events,
            type = "l")

file_path <- file.path(fpath, "plot_quad_time_h.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events,
            type = "h")

file_path <- file.path(fpath, "plot_quad_time0.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events,
            type = "l", days = days0)

file_path <- file.path(fpath, "plot_quad_time1.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events,
            type = "p", days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_quad_time_zoom.pdf")
plot_scores(scores_quad, times, model_names, model_colors, file_path, events,
            type = "p", tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic differences (to model 1)
## Positive values indicate superior performance
score_diffs <- as.matrix(scores_quad[ ,rep(1, n_mods-1)]
                         - scores_quad[ ,2:n_mods])

file_path <- file.path(fpath, "plot_quad_time_diff.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4)

file_path <- file.path(fpath, "plot_quad_time_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-0.1,0.1))

file_path <- file.path(fpath, "plot_quad_time_diff0.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-0.1,0.1),
                 days = days0)

file_path <- file.path(fpath, "plot_quad_time_diff1.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4, trim = c(-1,0.7),
                 days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_quad_time_diff_trim_zoom.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-3,0.5),
                 tlim = list(c(1,8,2016), c(31,5,2018)))


## Spatial Poisson scores
file_path <- file.path(fpath, "plot_spat_time1.pdf")
plot_scores(scores_spat, times, model_names, model_colors, file_path, events,
            type = "p", days = days1)

## Spatial Poisson differences (to model 1)
score_diffs <- as.matrix(scores_spat[ ,rep(1, n_mods-1)]
                         - scores_spat[ ,2:n_mods])

file_path <- file.path(fpath, "plot_spat_time_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4, days = days1,
                 trim = c(-20,10))

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_spat_time_diff_trim_zoom.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, days = days1,
                 trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))


## Mass Poisson scores
file_path <- file.path(fpath, "plot_mass_time.pdf")
plot_scores(scores_mass, times, model_names, model_colors, file_path, events,
            type = "l")

file_path <- file.path(fpath, "plot_mass_time_h.pdf")
plot_scores(scores_mass, times, model_names, model_colors, file_path, events,
            type = "h")

file_path <- file.path(fpath, "plot_mass_time0.pdf")
plot_scores(scores_mass, times, model_names, model_colors, file_path, events,
            type = "l", days = days0)

file_path <- file.path(fpath, "plot_mass_time1.pdf")
plot_scores(scores_mass, times, model_names, model_colors, file_path, events,
            type = "p", days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_mass_time_zoom.pdf")
plot_scores(scores_mass, times, model_names, model_colors, file_path, events,
            type = "p", tlim = list(c(1,8,2016), c(31,5,2018)))


## Mass Poisson differences (to model 1)
## Positive values indicate superior performance
score_diffs <- as.matrix(scores_mass[ ,rep(1, n_mods-1)]
                         - scores_mass[ ,2:n_mods])

file_path <- file.path(fpath, "plot_mass_time_diff.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4)

file_path <- file.path(fpath, "plot_mass_time_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-2,1))

file_path <- file.path(fpath, "plot_mass_time_diff0.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-2,1),
                 days = days0)

file_path <- file.path(fpath, "plot_mass_time_diff1.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4, trim = c(-20,20),
                 days = days1)

# Plot from August 2016 to May 2018
file_path <- file.path(fpath, "plot_mass_time_diff_trim_zoom.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-20,10),
                 tlim = list(c(1,8,2016), c(31,5,2018)))



## Plot of Murphy diagram
file_path <- file.path(fpath, "plot_Murphy_diag.pdf")
plot_elementary_scores(Murphy_diag, grd, model_names, model_colors, file_path,
                       "score")

## Plot of differences Murphy diagram
Murphy_diag_diff <- Murphy_diag[ ,rep(1, n_mods-1)] - Murphy_diag[ ,2:n_mods]
file_path <- file.path(fpath, "plot_Murphy_diag_diff.pdf")
plot_elementary_scores(Murphy_diag_diff, grd, model_names, model_colors,
                       file_path, "score differences", whichmods = 2:n_mods)



##############################################
#### Part II - Fully aggregated forecasts ####

# Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)


## Calculate daily scores for quadratic and Poisson loss
scores_pois_agg <- scores_quad_agg <- matrix(0, nrow = n_days, ncol = n_mods)
for (i in 1:n_mods) {
  scores_pois_agg[ ,i] <- S_pois(models_agg[[i]], obs_agg) 
  scores_quad_agg[ ,i] <- S_quad(models_agg[[i]], obs_agg) 
}

## Murphy diagrams
# Include climatology: Have to switch
# to reduced testing region
subset_index <- region_intersect(clima, cells)

## subset forecast models and observations
clima <- clima[subset_index$clima, ]
obs <- obs[ ,subset_index$model]
for (i in 1:n_mods) models[[i]] <- models[[i]][ ,subset_index$model]

# Aggregate mean forecasts again, as we
# now consider a smaller testing region
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# use climatology as last model
models_agg[[n_mods+1]] <- rep(sum(clima$RATE), n_days)

# PAV-transformed forecasts and climatology
# Attention: this changes the ordering of the forecasts
models_pav <- list()
for (i in 1:(n_mods+1)) models_pav[[i]] <- isoreg(models_agg[[i]], obs_agg)$yf


# Define suitable grid after checking the quantiles of
# the forecasts models (on logarithmic scale)
n_theta <- 100
log_grid <- seq(-6, 2, len = n_theta)
grd <- exp(log_grid)

Murphy_diag_agg <- Murphy_diag_pav <- matrix(0, nrow = n_theta, ncol = n_mods+1)
for (i in 1:(n_mods+1)) {
  Murphy_diag_agg[ ,i] <- colMeans( S_theta_vec(models_agg[[i]], obs_agg, grd) )
  Murphy_diag_pav[ ,i] <- colMeans( S_theta_vec(models_pav[[i]],
                                                obs_agg[order(models_agg[[i]])],
                                                grd) )
}


## Time dependent plots of scores and score differences
## (now fully aggregated)

## Poisson
file_path <- file.path(fpath, "plot_pois_time_agg.pdf")
plot_scores(scores_pois_agg, times, model_names, model_colors, file_path,
            events)

file_path <- file.path(fpath, "plot_pois_time_agg0.pdf")
plot_scores(scores_pois_agg, times, model_names, model_colors, file_path,
            events, type = "l", days = days0)

file_path <- file.path(fpath, "plot_pois_time_agg1.pdf")
plot_scores(scores_pois_agg, times, model_names, model_colors, file_path,
            events, type = "p", days = days1)

## Poisson differences (to model 1)
## Positive values indicate superior performance
score_diffs <- as.matrix(scores_pois_agg[ ,rep(1, n_mods-1)]
                         - scores_pois_agg[ ,2:n_mods])

file_path <- file.path(fpath, "plot_pois_time_agg_diff.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4)

file_path <- file.path(fpath, "plot_pois_time_agg_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-2,1))

file_path <- file.path(fpath, "plot_pois_time_agg_diff1.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4, trim = c(-20,20),
                 days = days1)

## quadratic
file_path <- file.path(fpath, "plot_quad_time_agg.pdf")
plot_scores(scores_quad_agg, times, model_names, model_colors, file_path,
            events, logscale = F) 

## quadratic (log scale)
file_path <- file.path(fpath, "plot_quad_time_agg_log.pdf")
plot_scores(scores_quad_agg, times, model_names, model_colors, file_path,
            events)

## Quadratic differences (to model 1)
## Positive values indicate superior performance
score_diffs <- as.matrix(scores_quad_agg[ ,rep(1, n_mods-1)]
                         - scores_quad_agg[ ,2:n_mods])

file_path <- file.path(fpath, "plot_quad_time_agg_diff.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4)

file_path <- file.path(fpath, "plot_quad_time_agg_diff_trim.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "l", whichmods = 2:4, trim = c(-0.1,0.1))

file_path <- file.path(fpath, "plot_quad_time_agg_diff1.pdf")
plot_score_diffs(score_diffs, times, model_names, model_colors, file_path,
                 events, type = "p", whichmods = 2:4, trim = c(-1,0.7),
                 days = days1)



## Plot of Murphy diagrams
## Now compared to climatology
model_colors_clima <- c(model_colors, "magenta")
model_names_clima <- c(model_names, "Clima")
model_names_clima_rc <- c(model_names_clima, paste0("RC_", model_names_clima))

## Murphy diagram 1 (only aggregated forecasts)
file_path <- file.path(fpath, "plot_Murphy_diag_agg.pdf")
plot_elementary_scores(Murphy_diag_agg, grd, model_names_clima,
                       model_colors_clima, file_path, "score", whichmods = 1:5)


## Murphy diagram 2 (aggregated and recalibrated forecasts)
file_path <- file.path(fpath, "plot_Murphy_diag_agg_pav.pdf")
plot_elementary_scores(cbind(Murphy_diag_agg, Murphy_diag_pav), grd,
                       model_names_clima_rc,
                       c(model_colors_clima, model_colors_clima), file_path,
                       "score", model_ltys = rep(1:2, ea = 5), whichmods = 1:10)


## Plot of mean forecasts over time
# Use plot_scores function although this means that ylab is wrong in these
# plots
mean_models <- matrix(unlist(models_agg), ncol = 5)

file_path <- file.path(fpath, "plot_mean_forecasts.pdf")
plot_scores(mean_models, times, model_names_clima, model_colors_clima,
            file_path, events, logscale = F)

file_path <- file.path(fpath, "plot_mean_forecasts_log.pdf")
plot_scores(mean_models, times, model_names_clima, model_colors_clima,
            file_path, events)

