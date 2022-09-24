## Spatial plots of scores, calibration, and
## discrimination over testing area
## Now including spatial AGGREGATION
## All maps are created as in eval3.R

# Path for figures
fpath <- "/media/myData/Plots/"
# Path for r scripts
rpath <- "/media/myData/Doks/Forschung/Code/Earthquakes_Italy"


# source functions for scores
source(file.path(rpath, "functions_eval.R"))
# source functions for plotting
source(file.path(rpath, "functions_plot.R"))
# source functions for data preparation
source(file.path(rpath, "functions_prep.R"))

max_agg <- 10

for (k in 1:max_agg) {
  # Do Data preparation
  source(file.path(rpath, "data_prep.R"))
  ncols <- 200
  pal <- rev(heat.colors(ncols))
  # compute neighbourhood matrix
  neighbourhood_matrix <- neigh_mat(cells, k)
  # Aggregate forecast models and observations
  for (i in 1:n_mods) {
    models[[i]] <- as.matrix( models[[i]] %*% neighbourhood_matrix )
    gc()
  }
  obs_agg <- obs %*% neighbourhood_matrix
  ## Compute values for maps
  MCB_map <- DSC_map <- matrix(0, nrow = n_cells, ncol = n_mods)
  for (i in 1:n_mods) {
    decomp <- bin_decomposition(models[[i]], obs_agg, scf = S_pois2)
    MCB_map[ ,i] <- decomp$MCB
    DSC_map[ ,i] <- decomp$DSC
  }
  UNC_map <- decomp$UNC
  SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = n_mods, nrow = n_cells)
  ## Create maps for scores
  lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
  file_path <- file.path(fpath, paste0("map_score_log_agg", k, ".pdf"))
  plot_multi_map(SCR_map, pal, cells, lims, ncols, file_path, evts = events)
  ## Create maps for miscalibration
  lims[1] <- min( log(MCB_map) ) - 0.05
  lims[2] <- max( log(MCB_map) ) + 0.05
  file_path <- file.path(fpath, paste0("map_MCB_log_agg", k, ".pdf"))
  plot_multi_map(MCB_map, pal, cells, lims, ncols, file_path, evts = events)
  ## Create maps for discrimination
  offset <- 1e-5
  lims[1] <- log(offset)
  lims[2] <- max( log(DSC_map + offset)) + 0.1
  file_path <- file.path(fpath, paste0("map_DSC_log_agg", k, ".pdf"))
  plot_multi_map(DSC_map, pal, cells, lims, ncols, file_path, offset = offset)
  ## Create maps for score differences
  root_name <- "map_score_diff"
  pal <- c(0, 0.66)    # red and blue
  for (i in 1:4) {
    for (j in 1:4) {
      if (i >= j) next
      file_path <- file.path(fpath,
                             paste0(root_name, "_", i, j, "_agg", k, ".pdf"))
      plot_diffs_map(SCR_map[ ,i] - SCR_map[ ,j], pal, cells, file_path)
    }
  }
  # Clear whole workspace
  if (k < max_agg) {
    rm(models)
    gc()
  }
}


