## Spatial plots of scores, calibration,
## discrimination, and skill over testing area

# Do Data preparation
source('~/Documents/Code/Earthquakes_Italy/data_prep.R')
# source functions for scores
source('~/Documents/Code/Earthquakes_Italy/functions_eval.R')
# source functions for plotting
source('~/Documents/Code/Earthquakes_Italy/functions_plot.R')
# source functions for data preparation
source('~/Documents/Code/Earthquakes_Italy/functions_prep.R')

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots_TEST"


## Compute values for maps
MCB_map <- DSC_map <- matrix(0, nrow = ncells, ncol = nmods)
for (i in 1:nmods) {
  decomp <- bin_decomp(models[[i]], obs, scf = Spois2)
  MCB_map[ ,i] <- decomp$MCB
  DSC_map[ ,i] <- decomp$DSC
}
UNC_map <- decomp$UNC
SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = nmods, nrow = ncells)

## Do spatial plots
# create color vector 
#pal <- paste("gray", round(scl * 100), sep = "")   # white = high score
ncols <- 200
pal <- rev(heat.colors(ncols))
#pal <- gray.colors(ncols)

## Create maps for scores
lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
filePath <- paste(fpath, "map_score_log", sep = "/")
mapComparison(SCR_map, pal, cells, lims, ncols, filePath, evts = events)

## Create maps for MCB
lims[1] <- min( log(MCB_map) ) - 0.05
lims[2] <- max( log(MCB_map) ) + 0.05
filePath <- paste(fpath, "map_MCB_log", sep = "/")
mapComparison(MCB_map, pal, cells, lims, ncols, filePath, evts = events)

## Create maps for DSC
offs <- 1e-5
lims[1] <- log(offs)
lims[2] <- max( log(DSC_map + offs)) + 0.1
filePath <- paste(fpath, "map_DSC_log", sep = "/")
mapComparison(DSC_map, pal, cells, lims, ncols, filePath, offset = offs)

## Create maps for score differences
rootName <- paste(fpath, "map_score_diff", sep = "/")
pal <- c(0, 0.66)    # red and blue
for (i in 1:4) {
  for (j in 1:4) {
    if (i >= j) next
    filePath <- paste0(rootName, "_", i, j, ".pdf")
    mapDifferences(SCR_map[ ,i] - SCR_map[ ,j], pal, cells, filePath)
  }
}

## Create maps for skill scores
## Since we use the climatology the testing
## region changes we have to modify SCR_map

## Include climatology: Have to switch
## to reduced testing region
subs <- region_intersect(clima, cells)

# subset forecast models and observations
clima <- clima[subs$clima, ]
cells <- cells[subs$model, ]
ncells <- dim(cells)[1]
obs <- obs[ ,subs$model]
# subset events
events <- filterRegion(events, cells)

# Compute average scores for climatology
scrClima <- colMeans( Spois(matrix(clima$RATE, ncol = ncells, nrow = ndays, byrow = TRUE), obs) )
gc()
SCR_map <- SCR_map[subs$model, ]
SKL_map <- ( matrix(scrClima, ncol = nmods, nrow = ncells) - SCR_map ) / matrix(scrClima, ncol = nmods, nrow = ncells)
# Do plotting
# Specify colours for positive and negative 
# values (0.66 = blue, 0 = red, 0.35 = green)
pal <- c(0.35, 0)
lims <- c(-10, 1)
filePath <- paste(fpath, "map_skill_score.pdf", sep = "/")
mapSkills(SKL_map, pal, cells, lims, ncols, filePath, evts = events)
