## First steps for case study
## Data preparation


# Path for r scripts
rpath <- "./"

# Path for data
# The files containing the forecasts and observations should
# be located in this folder
dpath <- "./../data"

# Set file names (default names)
# The forecast model outputs are arrays with time in rows
# and grid cells in the columns
model_files <- c("forecast_ETAS_LM_FP32.dat.xz",
                 "forecast_ETES_FCM_FP32.dat.xz",
                 "forecast_STEP_LG_FP32.dat.xz",
                 "forecast-ensemble_SMA_FP32.dat.xz",
                 "forecast-ensemble_LR-WA_M3+s2-1d_FP32.dat.xz")
# Time stamps corresponding to model outputs
# (rows of the model output data)
time_stamps_file <- "meta_rows_dates.csv"
# Locations of grid cells corresponding to model outputs
# (columns of the model output data)
cell_file <- "meta_columns_cells.csv"
# Catalog of observed earthquakes
events_file <- "catalog.csv"

# Load necessary packages
library(Matrix)

# Load auxiliary functions
source(file.path(rpath, "functions_prep.R"))

# Set model names and their colors
model_names <- c("LM", "FCM", "LG", "SMA", "LRWA", "LRWAsb")
model_colors <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
last_day <- list(DD = 20, MM = 5, YY = 2020)

select_models <- c(1, 2, 3, 4, 5)
model_files <- model_files[select_models]
model_names <- model_names[select_models]


###############
## Load data ##

# Load time stamps for the models
file_path <- file.path(dpath, time_stamps_file)
res <- load_times(file_path, last_day)
times <- res$times

# Load list of model forecasts
file_paths <- file.path(dpath, model_files)
models <- load_models(file_paths, res$time_index)
n_mods <- length(models)

# Load the grid cell data (testing region)
file_path <- file.path(dpath, cell_file)
cells <- load_cells(file_path)

# Load data frame of M4+ events
file_path <- file.path(dpath, events_file)
events <- load_events(file_path, times)

# Filter the M4+ events for testing region
events <- filter_region(events, cells)

# Convert events data frame to observation matrix of the
# same format as the forecast matrices in models. Thus 
# any scoring function S can be applied to the components
# of these matrices. For example S( models[[i]], obs )
# gives a matrix of scores for all grid cells and days
n_cells <- dim(cells)[1]
n_days <- length(times)
obs <- observation_matrix(events, times, n_cells)

# Load climatological model (constant in time)
# clima_file <- file.path(dpath, "rate_clima.txt")
# clima <- read.table(clima_file, header = F, col.names = c("LON", "LAT", "RATE"))

# The climatology file contains only normalized rates, i.e. spatial
# distribution of events. We have to multiply it by the value of events
# per time period. Different choices (depending on how we include after
# shocks) are possible
#events_per7 <- 25.95 * 7/365      # See Mail by Warner (08.09.21)
# events_per7 <- 12 * 7/365      # See Mail by Warner (08.09.21)
#events_per7 <- 16.97 * 7/365      # Events per year (calculated from events file)
# clima$RATE <- clima$RATE * events_per7

## Clean up
rm(file_path)
