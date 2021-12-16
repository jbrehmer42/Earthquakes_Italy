## Spatial plots of scores, calibration, and
## discrimination over testing area
## BUT NOW including spatial aggregation


max_agg <- 10
#scores_pois_k <- matrix(0, nrow = max_agg + 1, ncol = 4)



for (k in 1:max_agg) {
  # Do Data preparation
  source('~/Documents/_temp/Case/first_steps.R')
  # Load evaluation and plotting functions
  source('~/Documents/_temp/Case/eval_functions.R')
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
  model5_agg <- as.matrix( matrix(mean(obs), ncol = nbins, nrow = ndays) %*% nmat )
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
  ## Save the average scores to see how
  ## they vary with aggregation level
  load('~/EQData/scores_pois_k.RData')
  scores_pois_k[k+1, ] <- colSums( SCR_map[ ,1:4] ) / (2*k + 1)^2
  save(scores_pois_k, file = '~/EQData/scores_pois_k.RData')
  ## Create all the maps as in eval3.R
  ## Create maps for the scores
  lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
  file <- "~/Documents/_temp/Case/Plots/map_score_log"
  file <- paste0(file, "_agg", k)
  plot_figure(SCR_map, pal, lims, ncols, file, evts = TRUE)
  ## Create maps for miscalibration
  lims[1] <- min( log(MCB_map) ) - 0.05
  lims[2] <- max( log(MCB_map) ) + 0.05
  file <- "~/Documents/_temp/Case/Plots/map_mcb_log"
  file <- paste0(file, "_agg", k)
  plot_figure(MCB_map, pal, lims, ncols, file, evts = TRUE)
  ## Create maps for discrimination
  offs <- 1e-5
  lims[1] <- log(offs)
  lims[2] <- max( log(DSC_map + offs)) + 0.1
  file <- "~/Documents/_temp/Case/Plots/map_dsc_log"
  file <- paste0(file, "_agg", k)
  plot_figure(DSC_map, pal, lims, ncols, file, offset = offs)
  ## Create maps for score differences
  rootName <- "~/Documents/_temp/Case/Plots/map_score_diff"
  for (i in 1:4) {
    for (j in 1:4) {
      if (i >= j) next
      file <- paste0(rootName, "_", i, j, "_agg", k)
      scr_abs <- 1.01 * max( abs(SCR_map[ ,i] - SCR_map[ ,j]) )
      plot_diff_maps(SCR_map, c(i,j), scr_abs, file)
    }
  }
  ## Create maps for skill scores
  skl_max <- 1
  skl_min <- -10
  filePath <- paste0("~/Documents/_temp/Case/Plots/map_skill_score_agg", k, ".pdf")
  source('~/Documents/_temp/Case/map_skill.R')
  # Clear whole workspace
  if (k < max_agg) {
    rm(list = ls())
    gc()
  }
}


## Do final plots
filePath <- "~/Documents/_temp/Case/Plots/plot_pois_k.pdf"
pdf(filePath, width = 8, height=6)
ylim <- c(min(scores_pois_k), max(scores_pois_k))
m <- dim(scores_pois_k)[1]
plot(0:(m-1), scores_pois_k[ ,1], ty = "l", xlab = "k", ylim = ylim,
     ylab = "score", main = "Scaled logarithmic scores for aggregated forecasts")
for (i in 2:4) lines(0:(m-1), scores_pois_k[ ,i], col = cols[i]) 
for (i in 1:4) points(0:(m-1), scores_pois_k[ ,i], col = cols[i], pch = 20)
legend(8.2, 0.645, mnames, col = cols, lwd = 2)
dev.off()

