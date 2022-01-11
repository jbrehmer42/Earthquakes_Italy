## Spatial plots of skill scores over testing area
## Now including spatial AGGREGATION


max_agg <- 8
#scores_pois_k <- matrix(0, nrow = max_agg + 1, ncol = 4)


for (k in 1:max_agg) {
  # Do Data preparation
  source('~/Documents/_temp/Case/first_steps.R')
  # Load plotting functions
  source('~/Documents/_temp/Case/eval_functions.R')
  # Switch to reduced testing region
  source('~/Documents/_temp/Case/clima_subset.R')
  # Load plotting functions
  source('~/Documents/_temp/Case/plot_functions.R')
  ncols <- 200
  pal <- rev(heat.colors(ncols))
  # compute neighbourhood matrix
  nmat <- neigh_mat(bins, k)
  # Aggregate the forecasts
  model1_agg <- as.matrix( model1 %*% nmat )
  rm(model1)
  gc()
  model2_agg <- as.matrix( model2 %*% nmat )
  rm(model2)
  gc()
  model3_agg <- as.matrix( model3 %*% nmat )
  rm(model3)
  gc()
  model4_agg <- as.matrix( model4 %*% nmat )
  rm(model4)
  gc()
  # climatology
  model5_agg <- as.matrix( matrix(clima$RATE, ncol = nbins, nrow = ndays, byrow = TRUE) %*% nmat )
  gc()
  # events
  obs_agg <- obs %*% nmat
  ## Compute the decomposition
  MCB_map <- DSC_map <- matrix(0, nrow = nbins, ncol = 5)
  for (i in 1:5) {
    mname <- paste0("model", i, "_agg")
    decomp <- bin_decomp(get(mname), obs_agg, scf = Spois2)
    MCB_map[ ,i] <- decomp$MCB
    DSC_map[ ,i] <- decomp$DSC
  }
  UNC_map <- decomp$UNC
  SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = 5, nrow = nbins)
  SKL_map <- ( matrix(SCR_map[ ,5], ncol = 4, nrow = nbins) - SCR_map[ ,1:4] ) / matrix(SCR_map[ ,5], ncol = 4, nrow = nbins)
  ## Create maps for skill scores
  skl_max <- 1
  skl_min <- -10
  filePath <- paste0("~/Documents/_temp/Case/Plots/map_skill_score_high_agg", k, ".pdf")
  source('~/Documents/_temp/Case/map_skill.R')
  # Clear whole workspace
  #if (k < max_agg) {
    rm(events, model1_agg, model2_agg, model3_agg, model4_agg, model5_agg, nmat, obs, obs_agg,
       decomp, DSC_map, MCB_map, SCR_map, SKL_map, times, times2)
    gc()
  #}
}

source('~/Documents/_temp/Case/clima_maps.R')
