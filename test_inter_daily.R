## Produce mean reliability diagrams including interdaily updates

# Path for figures
fpath <- "/home/jbrehmer/Documents/_temp/Case/Plots"
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

# Load auxiliary functions
source(file.path(rpath, "functions_prep.R"))

# Set model names and their colors
mnames <- c("LM", "FMC", "LG", "SMA")
mcols <- c("black", "darkgreen", "blue", "red")
# Define last day where model evaluation is possible (needed
# because we treat 7-day periods)
lastday <- list(DD = 20, MM = 5, YY = 2020)

## NEW: Do not ignore multiple model runs on one day
# Load time stamps for the models
load_times_new <- function(filePath, lastday) {
  # Read time stamps from file
  times <- read.table(filePath, col.names = c("DD", "MM", "YY", "H", "M", "S"))
  n <- dim(times)[1]
  tindex <- rep(T, n)
  # Adjust for last day for which data is available. All
  # days after this day will be deleted
  matchlastday <- (times$DD == lastday$DD) & (times$MM == lastday$MM) & (times$YY == lastday$YY)
  if (any(matchlastday)) tindex[max(which(matchlastday)+1):n] <- FALSE
  times <- times[tindex, ]
  # Return new time stamps and index of corresponding rows
  return(list(times = times, tindex = tindex))
}
filePath <- file.path(dpath, timestampName)
res <- load_times_new(filePath, lastday)
times <- res$times

# Load list of model forecasts
filePaths <- file.path(dpath, modelNames)
models <- load_models(filePaths, res$tindex)
nmods <- length(models)

# Load the grid cell data (testing region)
filePath <- file.path(dpath, cellName)
cells <- load_cells(filePath)

# NEW: Do not need a time index
# Load data frame of M4+ events
load_events_new <- function(filePath, times) {
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
  n <- dim(events)[1]
  Tindex <- rep(-1, n)
  for (i in 1:n) {
    DD <- events$DD[i]
    MM <- events$MM[i]
    YY <- events$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) Tindex[i] <- 1
  }
  # Erase events which do not occur during the testing period
  events <- events[(Tindex > 0), ]
  return(events)
}
filePath <- file.path(dpath, eventsName)
events <- load_events_new(filePath, times)

# Filter the M4+ events for testing region
events <- filterRegion(events, cells)

## NEW PART ##
# Count events in 7-day period after forecasts
# Produce aggregated observations (ignore spatial)
ndays <- dim(times)[1]
obs_agg <- rep(0, ndays)
event_days <- as.Date(paste(events$YY, events$MM, events$DD, sep = "-"))
for (i in 1:ndays) {
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
for (i in 1:nmods) models_agg[[i]] <- rowSums(models[[i]])


## Plot mean forecasts over time period
mean_models <- matrix(unlist(models_agg), ncol = 4)
filePath <- file.path(fpath, "plot_mean_forecasts_inter.pdf")
plotScores(mean_models, times, mnames, mcols, filePath, events, logscale = F)

filePath <- file.path(fpath, "plot_observations.pdf")
pdf(filePath, width = 7.5, height = 6)
plot(1:ndays, obs_agg)
dev.off()

############################################
## Create rel diagrams using new function ##

## Define modified version of reliability diagram function
# function to create a reliability diagram
plotReliability <- function(aggr, y, txt = "", col = "black", lim = NULL, ln = F, resamp = NULL) {
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
  mcb <- s - s_rc
  umcb = s - s_rc_ucond
  cmcb = s_rc_ucond - s_rc
  dsc <- s_mg - s_rc
  unc <- s_mg
  # Do resampling
  if (!missing(resamp)) {
    res <- y[order(aggr)] - x
    if (is.numeric(resamp)) n <- round(resamp) + 1 else n <- 5000
    low <- floor(n * (1 - 0.9)/2)
    #low2 <- floor(n * (1 - 0.7)/2)
    up <- n - low
    #up2 <- n - low2
    # Resample residuals under assumption of iid residuals
    resamples <- sapply(1:n, function(z) x + sample(res, length(y)))
    # use rounding for integer data 
    resamples <- round(pmax(resamples, 0))
    yf_resamples <- apply(resamples, 2, function(z) isoreg(x, z)$yf )
    # sort resamples, include observed values, and correct bias (shift by mean residual)
    yf_resamples_sorted <- apply(cbind(yf, yf_resamples), 1, sort) - mean(res)
    # Compute limits for plotting
    pind <- (x < yf_resamples_sorted[up, ]) & (x > yf_resamples_sorted[low, ])
    # Compute MCB p-value
    mcb_resamples <- sapply(1:n, function(i) score(x, resamples[ ,i]) - score(yf_resamples[ ,i], resamples[ ,i]))
    rank_obs <- rank(c(mcb_resamples, mcb))[n+1]
    pval <- 1 - (rank_obs - 1)/(n + 1)
  }
  # Prepare plots
  if (missing(lim)) {
    lim <- c(min(x), max(x))
    lim <- lim + c(-1,1) * diff(lim) * 0.08
  }
  ttl <- paste("Mean reliability", txt)
  if (ln) {
    # Modify for log-log plot
    yf_line <- yf[ (yf > 0) & pind ]
    x_line <- x[ (yf > 0) & pind ]
    if (missing(lim))  lim[1] <- max(1e-6, min(yf_line))
    plot(NULL, xlim = lim, ylim = lim, main = paste(ttl, "(log scale)"), xlab = "",
         ylab = "", log = "xy")
  } else {
    # Standard plot
    plot(NULL, xlim = lim, ylim = lim, main = ttl, xlab = "", ylab = "")
    x_line <- x[pind]
    yf_line <- yf[pind]
  }
  if (!missing(resamp)) {
    # Add consistency band
    polygon(c(x[pind], rev(x[pind])), c(yf_resamples_sorted[up,pind], rev(yf_resamples_sorted[low,pind])),
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
       labels = paste0(round(c(umcb, cmcb, dsc, unc, s), digits = 3), collapse = "\n"))
  if (!missing(resamp)) {
    text(x = lim[1] + 7 * offs, y = lim[2] - offs, adj = c(0,1),
         paste0("[p = ", round(pval, 3), "]"))
  }
  return(invisible(yf))
}


# Reliability diagram
filePath <- file.path(fpath, "mean_reldiag_inter.pdf")
lim1 <- c(0, 27)
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram (version 2)
filePath <- file.path(fpath, "mean_reldiag_inter2.pdf")
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = c(0, 6), resamp = T)
}
dev.off()

# Reliability diagram on log scale
filePath <- file.path(fpath, "mean_reldiag_inter_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(filePath, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:nmods) {
  plotReliability(models_agg[[i]], obs_agg, mnames[i], mcols[i], lim = lim2, ln = T, resamp = T)
}
dev.off()
