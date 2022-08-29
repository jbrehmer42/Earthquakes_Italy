## Spatial plots of skill scores over testing area
## Now including spatial AGGREGATION

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots"
# Path for R code
rpath <- "/home/jbrehmer/Documents/Code/Earthquakes_Italy"

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
  subs <- region_intersect(clima, cells)
  # Subset climatology
  clima <- clima[subs$clima, ]
  # subset cells, observations, and models
  cells <- cells[subs$model, ]
  ncells <- dim(cells)[1]
  obs <- obs[ ,subs$model]
  for (i in 1:nmods) {
    models[[i]] <- models[[i]][ ,subs$model]
    gc()
  }
  # subset events
  events <- filterRegion(events, cells)
  # Add climatology as last model
  models[[nmods+1]] <- matrix(clima$RATE, ncol = ncells, nrow = ndays, byrow = TRUE)
  # compute neighbourhood matrix
  nmat <- neigh_mat(cells, k)
  # Aggregate forecast models and observations
  for (i in 1:(nmods+1)) {
    models[[i]] <- as.matrix( models[[i]] %*% nmat )
    gc()
  }
  obs_agg <- obs %*% nmat
  # Compute average scores
  SCR_map <- matrix(0, nrow = ncells, ncol = (nmods+1))
  for (i in 1:(nmods+1)) SCR_map[ ,i] <- colMeans( Spois(models[[i]], obs) )
  # Compute skills
  SKL_map <- ( matrix(SCR_map[ ,nmods+1], ncol = nmods, nrow = ncells) - SCR_map[ ,1:nmods] ) / matrix(SCR_map[ ,nmods+1], ncol = nmods, nrow = ncells)
  # Plots
  ncols <- 200
  pal <- c(0.35, 0)
  lims <- c(-10, 1) 
  filePath <- file.path(fpath, paste0("map_skill_score_agg", k, ".pdf"))
  mapSkills(SKL_map, pal, cells, lims, ncols, filePath, evts = events)
  if (k < max_agg) {
    rm(models, nmat)
    gc()
  }
}
