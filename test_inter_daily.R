## Produce mean reliability diagrams including interdaily updates

# Path for figures
fpath <- "/media/myData/Plots/"
# Path for r scripts
rpath <- "/media/myData/Doks/Forschung/Code/Earthquakes_Italy"


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

# Load auxiliary functions
source(file.path(rpath, "functions_prep.R"))

# Set model names and their colors
model_names <- c("LM", "FMC", "LG", "SMA")
model_colors <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
last_day <- list(DD = 20, MM = 5, YY = 2020)

## Difference to other evaluations: Do not ignore multiple model runs on one day
# Load time stamps for the models
load_times_new <- function(file_path, last_day) {
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
  # times    - Time stamps of testing period
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


## Plot mean forecasts over time period
mean_models <- matrix(unlist(models_agg), ncol = 4)
file_path <- file.path(fpath, "plot_mean_forecasts_inter.pdf")
plot_scores(mean_models, times, model_names, model_colors, file_path, events,
            logscale = F)
file_path <- file.path(fpath, "plot_observations.pdf")
pdf(file_path, width = 7.5, height = 6)
plot(1:n_days, obs_agg)
dev.off()

############################################
## Create rel diagrams using new function ##

## Define modified version of reliability diagram function
# function to create a reliability diagram
plot_reliability <- function(aggr, y, txt = "", col = "black", lim = NULL,
                             ln = F, resamp = NULL) {
  # Plots a mean reliability diagram similar
  # to the one in Gneiting and Resin (2021)
  # Resampling is done, but this is not as
  # straightforward, since all forecasts have
  # to be nonnegative!
  yf <- isoreg(aggr, y)$yf
  x <- sort(aggr)
  res <- y - aggr
  # calculate score values
  score <- function(x,y) mean((x-y)^2)
  s <- score(aggr, y)
  s_rc <- score(yf, y[order(aggr)])
  s_mg <- score(mean(y), y)
  c_rc_ucond = optim(par = 0, fn = function(c) score(aggr + c, y), 
                     method = "Brent", lower = min(res), upper = max(res))$par
  s_rc_ucond = score(aggr + c_rc_ucond, y)
  # calculate decomposition
  MCB <- s - s_rc
  uMCB = s - s_rc_ucond
  cMCB = s_rc_ucond - s_rc
  DSC <- s_mg - s_rc
  UNC <- s_mg
  # Do resampling
  if (!missing(resamp)) {
    res <- y[order(aggr)] - x
    if (is.numeric(resamp)) n <- round(resamp) + 1 else n <- 5000
    lower_limit <- floor(n * (1 - 0.9)/2)
    upper_limit <- n - lower_limit
    # Resample residuals under assumption of iid residuals
    resamples <- sapply(1:n, function(z) x + sample(res, length(y)))
    # use rounding for integer data 
    resamples <- round(pmax(resamples, 0))
    yf_resamples <- apply(resamples, 2, function(z) isoreg(x, z)$yf )
    # sort resamples, include observed values, and correct bias (shift by mean residual)
    yf_resamples_sorted <- apply(cbind(yf, yf_resamples), 1, sort) - mean(res)
    # Compute limits for plotting
    pind <- (x < yf_resamples_sorted[upper_limit, ]) & (x > yf_resamples_sorted[lower_limit, ])
    # Compute MCB p-value
    MCB_resamples <- sapply(1:n, function(i) score(x, resamples[ ,i]) - score(yf_resamples[ ,i], resamples[ ,i]))
    rank_obs <- rank(c(MCB_resamples, MCB))[n+1]
    pval <- 1 - (rank_obs - 1)/(n + 1)
  }
  # Prepare plots
  if (missing(lim)) {
    lim <- c(min(x), max(x))
    lim <- lim + c(-1,1) * diff(lim) * 0.08
  }
  ttl <- paste("Mean reliability", txt)
  # Do plotting
  if (ln) {
    # Modify for log-log plot
    yf_line <- yf[ (yf > 0) & pind ]
    x_line <- x[ (yf > 0) & pind ]
    if (missing(lim))  lim[1] <- max(1e-6, min(yf_line))
    plot(NULL, xlim = lim, ylim = lim, main = paste(ttl, "(log scale)"),
         xlab = "", ylab = "", log = "xy")
  } else {
    # Standard plot
    plot(NULL, xlim = lim, ylim = lim, main = ttl, xlab = "", ylab = "")
    x_line <- x[pind]
    yf_line <- yf[pind]
  }
  if (!missing(resamp)) {
    # Add consistency band
    polygon(c(x[pind], rev(x[pind])),
            c(yf_resamples_sorted[upper_limit, pind],
              rev(yf_resamples_sorted[lower_limit, pind])),
            border = NA, col = "indianred1")
  }
  lines(x_line, yf_line, col = col, lwd = 2)
  abline(a = 0, b = 1, col = "grey40", lty = 2)
  lim <- par("usr")[1:2]
  # switch back to "normal" coordinates
  par(xlog = F)
  par(ylog = F)
  title(xlab = "means", ylab = "recalibrated means", line = 2)
  # add decomposition values
  offs <- 0.05 * diff(lim)
  text(x = lim[1] + offs, y = lim[2] - offs, adj = c(0,1),
       labels = paste0(c("uMCB", "cMCB", "DSC ", "UNC", "SCR"), collapse = "\n"))
  text(x = lim[1] + 4 * offs, y = lim[2] - offs, adj = c(0,1),
       labels = paste0(round(c(uMCB, cMCB, DSC, UNC, s), digits = 3), collapse = "\n"))
  if (!missing(resamp)) {
    text(x = lim[1] + 7 * offs, y = lim[2] - offs, adj = c(0,1),
         paste0("[p = ", round(pval, 3), "]"))
  }
  return(invisible(yf))
}


# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_inter.pdf")
lim1 <- c(0, 27)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram (version 2)
file_path <- file.path(fpath, "mean_reldiag_inter2.pdf")
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = c(0, 6), resamp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_inter_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T)
}
dev.off()
