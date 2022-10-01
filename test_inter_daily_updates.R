## Produce mean reliability diagrams including inter-daily updates
## Have to modify the data loading and preparation such that multiple model runs
## on one day are no longer ignored.

# Path for figures
fpath <- "/media/myData/Plots/"
# Path for r scripts
rpath <- "/media/myData/Doks/Forschung/Code/Earthquakes_Italy"


# Load auxiliary functions
source(file.path(rpath, "functions_prep.R"))
# source functions for plotting
source(file.path(rpath, "functions_plot.R"))


# Path for data
# The files containing the forecasts and observations should
# be located in this folder
dpath <- "/media/myData/EQData"

# Set file names (default names)
# The forecast model outputs are arrays with time in rows
# and grid cells in the columns
model_files <- c("ETAS_LM.txt.xz",
                 "ETES_FMC.txt.xz",
                 "STEP_LG.txt.xz",
                 "Bayesian_corr_27_10.txt.xz")
# Time stamps corresponding to model outputs
# (rows of the model output data)
time_stamps_file <- "meta_rows.txt"
# Locations of grid cells corresponding to model outputs
# (columns of the model output data)
cell_file <- "meta_column.csv"
# Catalog of observed earthquakes
events_file <- "meta_catalogo.txt"



# Set model names and their colors
model_names <- c("LM", "FMC", "LG", "SMA")
model_colors <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
last_day <- list(DD = 20, MM = 5, YY = 2020)

## Difference to other evaluations: Do not ignore multiple model runs on one day
# Load time stamps for the models
load_times_new <- function(file_path, last_day) {
  # Load time stamps of the model outputs
  #
  # Input values:
  # file_path  - File path for the time stamps file
  # last_day   - Last day for which forecasts can be
  #            evaluated.
  # Read time stamps from file
  times <- read.table(file_path, col.names = c("DD", "MM", "YY", "H", "M", "S"))
  n <- dim(times)[1]
  time_index <- rep(T, n)
  # Adjust for last day for which data is available. All
  # days after this day will be deleted
  match_last_day <- (times$DD == last_day$DD) &
                      (times$MM == last_day$MM) &
                      (times$YY == last_day$YY)
  if (any(match_last_day)) time_index[max(which(match_last_day)+1):n] <- FALSE
  times <- times[time_index, ]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, time_index = time_index))
}
file_path <- file.path(dpath, time_stamps_file)
res <- load_times_new(file_path, last_day)
times <- res$times

# Load list of model forecasts
file_paths <- file.path(dpath, model_files)
models <- load_models(file_paths, res$time_index)
n_mods <- length(models)

# Load the grid cell data (testing region)
file_path <- file.path(dpath, cell_file)
cells <- load_cells(file_path)

# Difference to other evaluations: Do not need a time index
# Load data frame of M4+ events
load_events_new <- function(file_path, times) {
  # Input values:
  # file_path - File path for the events file
  # times     - Time stamps of testing period
  # Read events file
  events <- read.table(file_path,
                       col.names = c("YY", "MM", "DD", "H", "M", "S", "LAT",
                                     "LON", "DEP", "MAG"))
  # Use only events with magnitude M >= 4
  # Using M >= 4 is equivalent to using M >= 3.95
  M4ind <- (events$MAG >= 4)
  events <- events[M4ind, ]
  n <- dim(events)[1]
  time_index <- rep(-1, n)
  for (i in 1:n) {
    DD <- events$DD[i]
    MM <- events$MM[i]
    YY <- events$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) time_index[i] <- 1
  }
  # Erase events which do not occur during the testing period
  events <- events[(time_index > 0), ]
  return(events)
}
file_path <- file.path(dpath, events_file)
events <- load_events_new(file_path, times)

# Filter the M4+ events for testing region
events <- filter_region(events, cells)

## NEW PART ##
# Count events in 7-day period after forecasts
# Produce aggregated observations (ignore spatial behavior)
n_days <- dim(times)[1]
obs_agg <- rep(0, n_days)
event_days <- as.Date(paste(events$YY, events$MM, events$DD, sep = "-"))
for (i in 1:n_days) {
  # Add seven days to current day
  day <- as.Date(paste(times$YY[i], times$MM[i], times$DD[i], sep = "-"))
  day_end <- day + 7
  # find days
  ind <- (event_days >= day) & (event_days < day_end)
  events_select <- events[event_days == day_end, ]
  if (dim(events_select)[1] > 0) {
    # check time
    ind_time <- (events_select$H < times$H[i]) | 
      ((events_select$H == times$H[i]) & (events_select$M < times$M[i]) )
    obs_agg[i] <- sum(ind) + sum(ind_time)
  } else {
    obs_agg[i] <- sum(ind)
  }
}


# Aggregate models
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])


# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_inter1.pdf")
lim1 <- c(0, 27)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T, MCB_decomp = T)
}
dev.off()

# Reliability diagram (version 2)
file_path <- file.path(fpath, "mean_reldiag_inter2.pdf")
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = c(0, 6), resamp = T, MCB_decomp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_inter_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T, MCB_decomp = T)
}
dev.off()
