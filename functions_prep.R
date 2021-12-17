############################################
## Auxiliary functions - Data preparation ##


# Load time stamps of the model outputs
load_times <- function(filePath, lastday) {
  # Input values:
  # filePath - File path for the time stamps file
  # lastday  - Last day for which forecasts can be
  #            evaluated.
  # Read time stamps from file
  times <- read.table(filePath, col.names = c("DD", "MM", "YY", "H", "M", "S"))
  # Filter out days with multiple model runs: Compute
  # index of days with only one model run. Only the
  # first time stamp of the model runs is retained
  n <- dim(times)[1]
  tindex <- rep(T, n)
  for (i in 2:n) {
    # Check whether previous day agrees with current day
    tindex[i] <- !all(times[i-1, 1:3] == times[i, 1:3])
  }
  # Adjust for last day for which data is available. All
  # days after this day will be deleted
  matchlastday <- (times$DD == lastday$DD) & (times$MM == lastday$MM) & (times$YY == lastday$YY)
  if (any(matchlastday)) tindex[max(which(matchlastday)+1):length(tindex)] <- FALSE
  times <- times[tindex, ]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, tindex = tindex))
}

# Load forecast values of the models
load_models <- function(filePaths, tindex) {
  # Input values:
  # filePaths - File paths for the model outputs
  # tindex    - Index of model runs (rows) to be used
  # Assume that the file tindex applies to all
  # forecast model files
  # Read forecast values from files
  nmod <- length(filePaths)
  models <- list()
  for (i in 1:nmod) {
    # Rows are days
    # Columns are grid cells
    models[[i]] <- as.matrix( read.table( xzfile(filePaths[i]) ) )
    attr(models[[i]], "dimnames") <- NULL
    # delete rows corresponding to multiple model runs on one day
    models[[i]] <- models[[i]][tindex, ]
  }
  # Return list of model output matrices
  return(models)
}


# Load grid cells (testing region)
load_cells <- function(filePath) {
  # Input values:
  # filePath - File path for the grid cells file
  # Read grid cell file
  cells <- read.csv(filePath, header = F, col.names = c("LON", "LAT", "N"))
  # Start numbering the grid cells with 1 (just for
  # convenience and interpretability)
  cells$N <- cells$N + 1
  # Add x-y-coordinates for cells
  xvals <- sort(unique(cells$LON))
  yvals <- sort(unique(cells$LAT))
  cells$X <- sapply(cells$LON, function(x) which(x == xvals))
  cells$Y <- sapply(cells$LAT, function(x) which(x == yvals))
  return(cells)
}


# Load events data.frame
load_events <- function(filePath, times) {
  # Input values:
  # filePath - File path for the events file
  # times    - Time stamps of testing period
  # Read events file
  events <- read.table(filePath,
                       col.names = c("YY", "MM", "DD", "H", "M", "S", "LAT", "LON", "DEP", "MAG"))
  # Use only events with magnitude M >= 4
  # Using M >= 4 is equivalent to using M >= 3.95
  M4ind <- (events$MAG >= 4)
  events <- events[M4ind, ]
  # Assign the day number (time index TI) to events. Days
  # of the testing period (times) are consecutively numbered
  n <- dim(events)[1]
  Tindex <- rep(-1, n)
  for (i in 1:n) {
    DD <- events$DD[i]
    MM <- events$MM[i]
    YY <- events$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) Tindex[i] <- which(gind)
  }
  # Erase events which do not occur during the testing period
  events <- events[(Tindex > 0), ]
  # Add column time index (TI) to the events data frame
  events$TI <- Tindex[Tindex > 0]
  return(events)
}


# Filter out events which do not fall into the 
# testing region. Assign a grid cell number to
# the remaining events
filterRegion <- function(events, cells) {
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
  ind <- rep(0, n)
  for (i in 1:n) {
    cell_ri <- cells$LON + 0.5 * size_LON
    cell_le <- cells$LON - 0.5 * size_LON
    cell_lo <- cells$LAT - 0.5 * size_LAT
    cell_up <- cells$LAT + 0.5 * size_LAT
    isLON <- (cell_le < events$LON[i]) & (events$LON[i] <= cell_ri)
    isLAT <- (cell_lo < events$LAT[i]) & (events$LAT[i] <= cell_up)
    if (any(isLON & isLAT) ) {
      # Assign cell number inside testing region
      ind[i] <- cells$N[isLON & isLAT]
    } else {
      # Set to -1 outside testing region
      ind[i] <- -1
    }
  }
  # Filter events and add cell number N
  events <- events[(ind > 0), ]
  events$N <- ind[ind > 0]
  return(events)
}


# Create an observation matrix from the events
# data frame
events2obs <- function(events, ndays, ncells) {
  # Input values:
  # events - Data frame of observed events
  # ndays  - Number of days in testing period
  # ncells - Number of cells in testing region
  # Create an observation matrix which can be directly
  # compared to the forecast model output matrices
  # Rows are days
  # Columns are grid cells
  obs <- Matrix(0, ncol = ncells, nrow = ndays, sparse = T)
  for (i in 1:ndays) {
    # Collect events in a 7-day period
    ind <- (events$TI >= i) & (events$TI < i + 7)
    if (any(ind)) {
      obs[i, ] <- tabulate(events$N[ind], nbins = ncells)
    } else next
  }
  return(obs)
}