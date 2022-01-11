## Spatial plots of scores, calibration, and
## discrimination over testing area
## Now including spatial AGGREGATION
## All maps are created as in eval3.R

# source functions for scores
source('~/Documents/Code/Earthquakes_Italy/functions_eval.R')
# source functions for plotting
source('~/Documents/Code/Earthquakes_Italy/functions_plot.R')
# source functions for data preparation
source('~/Documents/Code/Earthquakes_Italy/functions_prep.R')

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots_TEST"

max_agg <- 10

for (k in 1:max_agg) {
  # Do Data preparation
  source('~/Documents/Code/Earthquakes_Italy/data_prep.R')
  ncols <- 200
  pal <- rev(heat.colors(ncols))
  # compute neighbourhood matrix
  nmat <- neigh_mat(cells, k)
  # Aggregate forecast models and observations
  for (i in 1:nmods) {
    models[[i]] <- as.matrix( models[[i]] %*% nmat )
    gc()
  }
  obs_agg <- obs %*% nmat
  ## Compute values for maps
  MCB_map <- DSC_map <- matrix(0, nrow = ncells, ncol = nmods)
  for (i in 1:nmods) {
    decomp <- bin_decomp(models[[i]], obs_agg, scf = Spois2)
    MCB_map[ ,i] <- decomp$MCB
    DSC_map[ ,i] <- decomp$DSC
  }
  UNC_map <- decomp$UNC
  SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = nmods, nrow = ncells)
  ## Create maps for scores
  lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
  filePath <- paste0(fpath, "/", "map_score_log_agg", k, ".pdf")
  mapComparison(SCR_map, pal, cells, lims, ncols, filePath, evts = events)
    ## Create maps for miscalibration
  lims[1] <- min( log(MCB_map) ) - 0.05
  lims[2] <- max( log(MCB_map) ) + 0.05
  filePath <- paste0(fpath, "/", "map_MCB_log_agg", k, ".pdf")
  mapComparison(MCB_map, pal, cells, lims, ncols, filePath, evts = events)
  ## Create maps for discrimination
  offs <- 1e-5
  lims[1] <- log(offs)
  lims[2] <- max( log(DSC_map + offs)) + 0.1
  filePath <- paste0(fpath, "/", "map_DSC_log_agg", k, ".pdf")
  mapComparison(DSC_map, pal, cells, lims, ncols, filePath, offset = offs)
  ## Create maps for score differences
  rootName <- paste(fpath, "map_score_diff", sep = "/")
  pal <- c(0, 0.66)    # red and blue
  for (i in 1:4) {
    for (j in 1:4) {
      if (i >= j) next
      filePath <- paste0(rootName, "_", i, j, "_agg", k, ".pdf")
      mapDifferences(SCR_map[ ,i] - SCR_map[ ,j], pal, cells, filePath)
    }
  }
  # Clear whole workspace
  if (k < max_agg) {
    rm(models)
    gc()
  }
}

# filePath <- "~/Documents/_temp/Case/Plots/plot_pois_k.pdf"
# pdf(filePath, width = 8, height=6)
# ylim <- c(min(scores_pois_k), max(scores_pois_k))
# m <- dim(scores_pois_k)[1]
# plot(0:(m-1), scores_pois_k[ ,1], ty = "l", xlab = "k", ylim = ylim,
#      ylab = "score", main = "Scaled logarithmic scores for aggregated forecasts")
# for (i in 2:4) lines(0:(m-1), scores_pois_k[ ,i], col = cols[i]) 
# for (i in 1:4) points(0:(m-1), scores_pois_k[ ,i], col = cols[i], pch = 20)
# legend(8.2, 0.645, mnames, col = cols, lwd = 2)
# dev.off()

