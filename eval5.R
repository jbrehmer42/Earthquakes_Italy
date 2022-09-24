## Spatial plots of skill scores over testing area
## Now including spatial AGGREGATION

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


max_agg <- 8

for (k in 1:max_agg) {
  # Do Data preparation
  source(file.path(rpath, "data_prep.R"))
  ## Switch to reduced testing region
  subset_index <- region_intersect(clima, cells)
  # Subset climatology
  clima <- clima[subset_index$clima, ]
  # subset cells, observations, and models
  cells <- cells[subset_index$model, ]
  n_cells <- dim(cells)[1]
  obs <- obs[ ,subset_index$model]
  for (i in 1:n_mods) {
    models[[i]] <- models[[i]][ ,subset_index$model]
    gc()
  }
  # subset events
  events <- filter_region(events, cells)
  # Add climatology as last model
  models[[n_mods+1]] <- matrix(clima$RATE, ncol = n_cells, nrow = n_days,
                               byrow = TRUE)
  # compute neighbourhood matrix
  neighbourhood_matrix <- neigh_mat(cells, k)
  # Aggregate forecast models and observations
  for (i in 1:(n_mods+1)) {
    models[[i]] <- as.matrix( models[[i]] %*% neighbourhood_matrix )
    gc()
  }
  obs_agg <- obs %*% neighbourhood_matrix
  # Compute average scores
  SCR_map <- matrix(0, nrow = n_cells, ncol = (n_mods+1))
  for (i in 1:(n_mods+1)) SCR_map[ ,i] <- colMeans( S_pois(models[[i]], obs) )
  # Compute skills
  SKL_map <- ( matrix(SCR_map[ ,n_mods+1], ncol = n_mods, nrow = n_cells)
               - SCR_map[ ,1:n_mods] ) / matrix(SCR_map[ ,n_mods+1],
                                                ncol = n_mods, nrow = n_cells)
  # Plots
  n_colors <- 200
  pal <- c(0.35, 0)
  lims <- c(-10, 1) 
  file_path <- file.path(fpath, paste0("map_skill_score_agg", k, ".pdf"))
  plot_skill_map(SKL_map, pal, cells, lims, n_colors, file_path, evts = events)
  if (k < max_agg) {
    rm(models, neighbourhood_matrix)
    gc()
  }
}
