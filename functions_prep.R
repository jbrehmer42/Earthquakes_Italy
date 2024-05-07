############################################
## Auxiliary functions - Data preparation ##

library(lubridate)    # work with dates

load_times <- function(file_path, last_day) {
  # Load time stamps of the model outputs
  #
  # Input values:
  # file_path  - File path for the time stamps file
  # last_day   - Last day for which forecasts can be
  #            evaluated.
  # Read time stamps from file
  times <- read.csv(file_path, col.names = "T")
  times <- ymd_hms(times$T)
  # Filter out days with multiple model runs: Compute
  # index of days with only one model run. Only the
  # first time stamp of the model runs is retained
  time_index <- (hour(times) == 0 & minute(times) == 0 & second(times) == 0)
  # Adjust for last day for which data is available. All
  # days after this day will be deleted
  before_last_day <- (times <= ymd(with(last_day, paste(YY, MM, DD))))

  times_filter <- time_index & before_last_day
  times <- times[times_filter]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, time_index = times_filter))
}


load_models <- function(file_paths, time_index) {
  # Load forecast values of the models
  #
  # Input values:
  # file_paths   - File paths for the model outputs
  # time_index   - Index of model runs (rows) to be used
  # Assume that the file time_index applies to all
  # forecast model files.
  # Read forecast values from files
  n_mods <- length(file_paths)
  models <- list()
  for (i in 1:n_mods) {
    # Rows are days
    # Columns are grid cells
    models[[i]] <- as.matrix( read.table( xzfile(file_paths[i]) ) )
    attr(models[[i]], "dimnames") <- NULL
    # delete rows corresponding to multiple model runs on one day
    models[[i]] <- models[[i]][time_index, ]
  }
  # Return list of model output matrices
  return(models)
}


load_cells <- function(file_path) {
  # Load grid cell information (specifies the testing region)
  #
  # Input values:
  # file_path  -  File path for the grid cells file
  # Read grid cell file
  cells <- read.csv(file_path, header = T, col.names = c("LON", "LAT"))
  # Start numbering the grid cells with 1 (just for
  # convenience and interpretability)
  cells$N <- 1:nrow(cells)
  # Add x-y-coordinates for cells
  xvals <- sort(unique(cells$LON))
  yvals <- sort(unique(cells$LAT))
  cells$X <- sapply(cells$LON, function(x) which(x == xvals))
  cells$Y <- sapply(cells$LAT, function(x) which(x == yvals))
  return(cells)
}


load_events <- function(file_path, times) {
  # Load catalog/events data frame and treat dates properly
  #
  # Input values:
  # file_path - File path for the events file
  # times     - Time stamps of testing period
  # Read events file
  events <- read.csv(file_path, col.names = c("TS", "LAT", "LON", "DEP", "MAG"))
  events$TS <- ymd_hms(events$TS)
  # Use only events with magnitude M >= 4
  # Using M >= 4 is equivalent to using M >= 3.95
  M4ind <- (events$MAG >= 4)
  events <- events[M4ind, ]
  # Assign the day number (time index TI) to events. Days
  # of the testing period (times) are consecutively numbered
  n <- dim(events)[1]
  time_index <- rep(-1, n)

  for (i in 1:n) {
    # filter for events lying in the 7-day period of a prediction
    gind <- (events$TS[i] >= times) & (events$TS[i] < times + days(7))
    # pick pred date closest to event date as time index
    if (any(gind)) time_index[i] <- max(which(gind))
  }
  # Erase events which do not occur during the testing period
  events <- events[(time_index > 0), ]
  # Add column time index (TI) to the events data frame
  events$TI <- time_index[time_index > 0]
  return(events)
}


filter_region <- function(events, cells) {
  # Filter out events which do not fall into the 
  # testing region. Assign a grid cell number to
  # the remaining events
  # In accordance with pycsep: "Bins are inclusive on lower bound
  # and exclusive on upper bound"
  #
  # Input values:
  # events - Data frame of observed events
  # cells  - Data frame of grid cells
  # Specify edge length of grid cells
  size_LON <- 0.1
  size_LAT <- 0.1
  # Assign cell number to every event. Cell numbers
  # start with 1. If an event does not fall into the
  # testing region its cell number is -1 (and it is
  # filtered out)
  n <- dim(events)[1]
  cell_number <- rep(0, n)
  # Set right, left, lower, and upper limits of cells
  cell_ri <- cells$LON + 0.5 * size_LON
  cell_le <- cells$LON - 0.5 * size_LON
  cell_lo <- cells$LAT - 0.5 * size_LAT
  cell_up <- cells$LAT + 0.5 * size_LAT
  # need tolerance due to floating point inaccuracies (some events lie exactly on boundary which is not catched)
  tol <- 1e-12
  for (i in 1:n) {
    isLON <- (cell_le < events$LON[i] + tol) & (events$LON[i] + tol < cell_ri)
    isLAT <- (cell_lo < events$LAT[i] + tol) & (events$LAT[i] + tol < cell_up)
    if (any(isLON & isLAT) ) {
      # Assign cell number inside testing region
      cell_number[i] <- cells$N[isLON & isLAT]
    } else {
      # Set to -1 outside testing region
      cell_number[i] <- -1
    }
  }
  # Filter events and add cell number N
  events <- events[(cell_number > 0), ]
  events$N <- cell_number[cell_number > 0]
  return(events)
}


observation_matrix <- function(events, times, n_cells) {
  # Create an observation matrix from the events
  # data frame
  #
  # Input values:
  # events   - Data frame of observed events
  # times    - Data frame of dates of forecasts
  # n_cells  - Number of cells in testing region
  # Create an observation matrix which can be directly
  # compared to the forecast model output matrices.
  # Rows are days, columns are grid cells.
  obs <- Matrix(0, ncol = n_cells, nrow = length(times), sparse = T)

  for (i in 1:length(times)) {
    # Collect events in a 7-day period
    is_in_period <- (events$TS >= times[i]) & (events$TS < times[i] + days(7))
    if (any(is_in_period)) {
      obs[i, ] <- tabulate(events$N[is_in_period], nbins = n_cells)
    }
  }
  return(obs)
}

count_missing_days <- function(times) {
  # Count how many dates are missing, assuming times spans a contagious timespan
  # of days
  #
  # Input values:
  # times - Data frame of dates of forecasts
  diffs <- times[-1] - times[-nrow(times)]
  return(sum(diffs > days(1)))
}

recycle_forecasts <- function(models, times) {
  # Recycle forecasts of previous days to impute missing forecasts
  #
  # Input values:
  # models - List of forecast matrices (time x space)
  # times  - Data frame of dates of forecasts
  #
  cont_dates <- times[1] + days(1):(times[length(times)] - times[1])
  new_n <- length(cont_dates)

  map_to_last <- sapply(1:new_n, function(i) max(which(cont_dates[i] >= times)))

  for (i in 1:n_mods) {
    models[[i]] <- models[[i]][map_to_last, ]
  }
  return(list(models = models, times = cont_dates))
}
