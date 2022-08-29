## First steps for case study
## Data preparation


# Path for r scripts
rpath <- "/home/jbrehmer/Documents/Code/Earthquakes_Italy"

# Path for data
# The files containing the forecasts and observations should
# be located in this folder
dpath <- "/home/jbrehmer/EQData"

# Set file names (default names)
# The forecast model outputs are arrays with time in rows
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

# Load auxiliary functions
source(file.path(rpath, "functions_prep.R"))

# Set model names and their colors
mnames <- c("LM", "FMC", "LG", "SMA")
mcols <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
lastday <- list(DD = 20, MM = 5, YY = 2020)


###############
## Load data ##

# Load time stamps for the models
filePath <- file.path(dpath, timestampName)
res <- load_times(filePath, lastday)
times <- res$times

# Load list of model forecasts
filePaths <- file.path(dpath, modelNames)
models <- load_models(filePaths, res$tindex)
nmods <- length(models)

# Load the grid cell data (testing region)
filePath <- file.path(dpath, cellName)
cells <- load_cells(filePath)

# Load data frame of M4+ events
filePath <- file.path(dpath, eventsName)
events <- load_events(filePath, times)

# Filter the M4+ events for testing region
events <- filterRegion(events, cells)

# Convert events data frame to observation matrix of the
# same format as the forecast matrices in models. Thus 
# any scoring function S can be applied to the components
# of these matrices. For example S( models[[i]], obs )
# gives a matrix of scores for all grid cells and days
ncells <- dim(cells)[1]
ndays <- dim(times)[1]
obs <- events2obs(events, ndays, ncells)

# Load climatological model (constant in time)
cfile <- file.path(dpath, "rate_clima.txt")
clima <- read.table(cfile, header = F, col.names = c("LON", "LAT", "RATE"))
#evts_per7 <- 25.95 * 7/365      # See Mail by Warner (08.09.21)
evts_per7 <- 12 * 7/365      # See Mail by Warner (08.09.21)
#evts_per7 <- 16.97 * 7/365      # See events per year (calculated)
clima$RATE <- clima$RATE * evts_per7

## Clean up
rm(filePath, cfile)
