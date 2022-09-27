## Spatial plots of scores, calibration,
## discrimination, and skill over testing area

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
# source functions for data preparation
source(file.path(rpath, "functions_prep.R"))


## Compute values for maps
MCB_map <- DSC_map <- matrix(0, nrow = n_cells, ncol = n_mods)
MCB_quad <- DSC_quad <- UNC_quad <- rep(0, n_mods)
for (i in 1:n_mods) {
  # Compute decomposition for maps
  decomp <- cell_decomposition(models[[i]], obs, scf = S_pois2)
  MCB_map[ ,i] <- decomp$MCB
  DSC_map[ ,i] <- decomp$DSC
  # Compute decomposition for quadratic score (table or alternative maps)
  decomp_quad <- cell_decomposition(models[[i]], obs, scf = S_quad2)
  MCB_quad[i] <- sum(decomp_quad$MCB)
  DSC_quad[i] <- sum(decomp_quad$DSC)
  UNC_quad[i] <- sum(decomp_quad$UNC)
}
UNC_map <- decomp$UNC
SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = n_mods, nrow = n_cells)

## Print score decompositions
# Poisson score
(colSums(MCB_map))
(colSums(DSC_map))
(sum(UNC_map))
# Quadratic score
(MCB_quad)
(DSC_quad)
(UNC_quad)

## Do spatial plots
# create color vector from palette
n_colors <- 200
pal <- rev(heat.colors(n_colors))
#pal <- gray.colors(n_colors)

## Create maps for scores
lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
file_path <- file.path(fpath, "map_score_log.pdf")
plot_multi_map(SCR_map, pal, cells, lims, n_colors, file_path, evts = events)

## Create maps for MCB
lims[1] <- min( log(MCB_map) ) - 0.05
lims[2] <- max( log(MCB_map) ) + 0.05
file_path <- file.path(fpath, "map_MCB_log.pdf")
plot_multi_map(MCB_map, pal, cells, lims, n_colors, file_path, evts = events)

## Create maps for DSC
offset <- 1e-5
lims[1] <- log(offset)
lims[2] <- max( log(DSC_map + offset)) + 0.1
file_path <- file.path(fpath, "map_DSC_log.pdf")
plot_multi_map(DSC_map, pal, cells, lims, n_colors, file_path, offset = offset)

## Create maps for score differences
root_name <- "map_score_diff"
pal <- c(0, 0.66)    # red and blue
for (i in 1:4) {
  for (j in 1:4) {
    if (i >= j) next
    file_path <- file.path(fpath, paste0(root_name, "_", i, j, ".pdf"))
    plot_diffs_map(SCR_map[ ,i] - SCR_map[ ,j], pal, cells, file_path)
  }
}

## Create maps for skill scores
## Since we use the climatology the testing
## region changes we have to modify SCR_map

## Include climatology: Have to switch
## to reduced testing region
subset_index <- region_intersect(clima, cells)

# subset forecast models and observations
clima <- clima[subset_index$clima, ]
cells <- cells[subset_index$model, ]
n_cells <- dim(cells)[1]
obs <- obs[ ,subset_index$model]
# subset events
events <- filter_region(events, cells)

# Compute average scores for climatology
SCR_clima <- colMeans( S_pois(matrix(clima$RATE, ncol = n_cells, nrow = n_days,
                                     byrow = TRUE), obs) )

gc()
SCR_map <- SCR_map[subset_index$model, ]
SKL_map <- ( matrix(SCR_clima, ncol = n_mods, nrow = n_cells) - SCR_map ) /
            matrix(SCR_clima, ncol = n_mods, nrow = n_cells)
# Plots
# Specify colours for positive and negative 
# values (0.66 = blue, 0 = red, 0.35 = green)
pal <- c(0.35, 0)
lims <- c(-10, 1)
file_path <- file.path(fpath, "map_skill_score.pdf")
plot_skill_map(SKL_map, pal, cells, lims, n_colors, file_path, evts = events)
