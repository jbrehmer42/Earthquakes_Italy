## Different ways of plotting the daily scores

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/figures/"

# Path for data
# The file containing the forecasts and observations should
# be located in this folder
dpath <- "/home/jbrehmer/EQData"
###
### Set to getwd() later

# Set file names (default names)
# The forecast model outputs as arrays with time in rows
# and grid cells in the columns
modelNames <- c("ETAS_LM.txt.xz", "ETES_FMC.txt.xz",
                "STEP_LG.txt.xz", "Bayesian_corr_27_10.txt.xz")
# Time stamps corresponding to model outputs
# (rows of the model output data)
timestampName <- "meta_rows.txt"
# Locations of grid cells corresponding to model outputs
# (columns of the model output data)
cellName <- "meta_column.csv"
# Catalog of observed earthquakes
eventsName <- "meta_catalogo.txt"

# Load necessary packages
library(Matrix)
library(maps)


# Load auxiliary functions
source("/home/jbrehmer/Documents/_temp/JASA_version/Code/case_functions.R") 
# Set model names and their colors
mnames <- c("LM", "FMC", "LG", "SMA")
#mcols <- c("black", "darkgreen", "blue", "red")
mcols <- c("black", "darkgreen", "blue", "red")      # switch color 1 and 3
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
lastday <- list(DD = 20, MM = 5, YY = 2020)


###############
## Load data ##

# Load time stamps for the models
filePath <- paste0(dpath, "/", timestampName)
res <- load_times(filePath, lastday)
times <- res$times

# Load list of model forecasts
filePaths <- paste0(dpath, "/", modelNames)
models <- load_models(filePaths, res$tindex)
nmods <- length(models)

# Load the grid cell data (testing region)
filePath <- paste0(dpath, "/", cellName)
cells <- load_cells(filePath)

# Load data frame of M4+ events
filePath <- paste0(dpath, "/", eventsName)
events <- load_events(filePath, times)

# Filter the M4+ events for testing region
events <- filterRegion(events, cells)

# Convert events data frame to observation matrix of the
# same format as the forecast matrices in models
ncells <- dim(cells)[1]
ndays <- dim(times)[1]
obs <- events2obs(events, ndays, ncells)

days0 <- days1 <- 1:ndays
days0[(rowSums(obs) != 0)] <- NA
days1[(rowSums(obs) == 0)] <- NA

####################
## Plots and maps ##

Sspat <- function(x, y) -y * log(x)

####################
####################

# Compute daily and overall scores for the quadratic
# and the Poisson scoring function, see Section 5.2
scores_pois <- scores_quad <- scores_spat <- scores_mass <- matrix(0, nrow = ndays, ncol = nmods)
for (i in 1:nmods) {
  scores_pois[ ,i] <- rowSums( Spois(models[[i]], obs) )
  scores_quad[ ,i] <- rowSums( Squad(models[[i]], obs) )
  scores_spat[ ,i] <- rowSums( Sspat(models[[i]], obs) ) + rowSums(obs) * log(rowSums(models[[i]]))
  scores_mass[ ,i] <- Spois( rowSums(models[[i]]), rowSums(obs) )
}

######################
######################

## Poisson Scores
filePath <- paste(fpath, "plot_pois_time.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "l")

filePath <- paste(fpath, "plot_pois_time_h.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "h")

filePath <- paste(fpath, "plot_pois_time0.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- paste(fpath, "plot_pois_time0_h.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "h", days = days0)

filePath <- paste(fpath, "plot_pois_time1.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- paste(fpath, "plot_pois_time_zoom.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))

filePath <- paste(fpath, "plot_pois_time_zoom_h.pdf", sep = "/")
plotScores(scores_pois, times, mnames, mcols, filePath, events, type = "h",
           tlim = list(c(1,8,2016), c(31,5,2018)))



## Poisson differences (to model 1)
scorediffs <- as.matrix(scores_pois[ ,rep(1, nmods-1)] - scores_pois[ ,2:nmods])

filePath <- paste(fpath, "plot_pois_time_diff.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- paste(fpath, "plot_pois_time_diff_trim.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1))

filePath <- paste(fpath, "plot_pois_time_diff0.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1), days = days0)

filePath <- paste(fpath, "plot_pois_time_diff1.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-20,20), days = days1)

filePath <- paste(fpath, "plot_pois_time_diff_trim_zoom.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic Scores
filePath <- paste(fpath, "plot_quad_time.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "l")

filePath <- paste(fpath, "plot_quad_time_h.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "h")

filePath <- paste(fpath, "plot_quad_time0.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "l", days = days0)

filePath <- paste(fpath, "plot_quad_time0_h.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "h", days = days0)

filePath <- paste(fpath, "plot_quad_time1.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- paste(fpath, "plot_quad_time_zoom.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))

filePath <- paste(fpath, "plot_quad_time_zoom_h.pdf", sep = "/")
plotScores(scores_quad, times, mnames, mcols, filePath, events, type = "h",
           tlim = list(c(1,8,2016), c(31,5,2018)))


## Quadratic differences (to model 1)
scorediffs <- as.matrix(scores_quad[ ,rep(1, nmods-1)] - scores_quad[ ,2:nmods])

filePath <- paste(fpath, "plot_quad_time_diff.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- paste(fpath, "plot_quad_time_diff_trim.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-0.1,0.1))

filePath <- paste(fpath, "plot_quad_time_diff0.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-0.1,0.1), days = days0)

filePath <- paste(fpath, "plot_quad_time_diff1.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-1,0.7), days = days1)

filePath <- paste(fpath, "plot_quad_time_diff_trim_zoom.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-3,0.5), tlim = list(c(1,8,2016), c(31,5,2018)))




## Spatial Poisson scores
filePath <- paste(fpath, "plot_spat_time1.pdf", sep = "/")
plotScores(scores_spat, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- paste(fpath, "plot_spat_time1_h.pdf", sep = "/")
plotScores(scores_spat, times, mnames, mcols, filePath, events, type = "h", days = days1)

## Spatial Poisson differences (to model 1)
scorediffs <- as.matrix(scores_spat[ ,rep(1, nmods-1)] - scores_spat[ ,2:nmods])

filePath <- paste(fpath, "plot_spat_time_diff_trim.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               days = days1, trim = c(-20,10))

filePath <- paste(fpath, "plot_spat_time_diff_trim_zoom.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               days = days1, trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))


## Mass Poisson scores
filePath <- paste(fpath, "plot_mass_time.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "l")

filePath <- paste(fpath, "plot_mass_time_h.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "h")

filePath <- paste(fpath, "plot_mass_time0.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "l", days = days0)

# filePath <- paste(fpath, "plot_mass_time0_h.pdf", sep = "/")
# plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "h", days = days0)

filePath <- paste(fpath, "plot_mass_time1.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "p", days = days1)

filePath <- paste(fpath, "plot_mass_time_zoom.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "p",
           tlim = list(c(1,8,2016), c(31,5,2018)))

filePath <- paste(fpath, "plot_mass_time_zoom_h.pdf", sep = "/")
plotScores(scores_mass, times, mnames, mcols, filePath, events, type = "h",
           tlim = list(c(1,8,2016), c(31,5,2018)))

## Mass Poisson differences (to model 1)
scorediffs <- as.matrix(scores_mass[ ,rep(1, nmods-1)] - scores_mass[ ,2:nmods])

filePath <- paste(fpath, "plot_mass_time_diff.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4)

filePath <- paste(fpath, "plot_mass_time_diff_trim.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1))

filePath <- paste(fpath, "plot_mass_time_diff0.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-2,1), days = days0)

filePath <- paste(fpath, "plot_mass_time_diff1.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "p", whichmods = 2:4,
               trim = c(-20,20), days = days1)

filePath <- paste(fpath, "plot_mass_time_diff_trim_zoom.pdf", sep = "/")
plotScoreDiffs(scorediffs, times, mnames, mcols, filePath, events, type = "l", whichmods = 2:4,
               trim = c(-20,10), tlim = list(c(1,8,2016), c(31,5,2018)))



