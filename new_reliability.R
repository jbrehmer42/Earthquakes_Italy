## Checking for mean calibration
## Now excluding the large events

# Path for figures
fpath <- "/media/myData/Plots/"
# Path for r scripts
rpath <- "/media/myData/Doks/Forschung/Code/Earthquakes_Italy"

## Define modified version of reliability diagram function
# function to create a reliability diagram
plot_reliability <- function(aggr, y, txt = "", col = "black", lim = NULL, ln = F, resamp = NULL) {
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



# Do Data preparation
source(file.path(rpath, "data_prep.R"))
# source functions for scores
source(file.path(rpath, "functions_eval.R"))
# source functions for plotting
source(file.path(rpath, "functions_plot.R"))


## Exclude large earthquakes (first variant)
Mlimit <- 5.5
time_index_large <- events$TI[events$MAG > Mlimit]
time_index_large <- unique(time_index_large)
subset_index <- rep(TRUE, n_days)
for (i in 1:length(time_index_large)) {
  subset_index[(time_index_large[i]-6):time_index_large[i]] <- FALSE
}
n_days - sum(subset_index)

## Exclude large earthquakes (second variant)
Mlimit <- 5.5
time_index_large <- events$TI[events$MAG > Mlimit]
time_index_large <- unique(time_index_large)
subset_index <- rep(TRUE, n_days)
subset_index[c(time_index_large+1, time_index_large + 2)] <- FALSE
n_days - sum(subset_index)

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)
# Erase days with large events
obs_agg <- obs_agg[subset_index]
for (i in 1:n_mods) models_agg[[i]] <- models_agg[[i]][subset_index]


## Reliability diagrams for aggregated number
## of earthquakes (sum over all bins)

# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_ex1.pdf")
lim1 <- c(0, 4.5)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_ex1_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T)
}
dev.off()

# Reliability diagram
file_path <- file.path(fpath, "mean_reldiag_new.pdf")
lim1 <- c(0, 6)
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim1, resamp = T)
}
dev.off()

# Reliability diagram on log scale
file_path <- file.path(fpath, "mean_reldiag_new_log.pdf")
lim2 <- c(0.05, lim1[2])
pdf(file_path, width = 7.5, height = 8)
par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
for (i in 1:n_mods) {
  plot_reliability(models_agg[[i]], obs_agg, model_names[i], model_colors[i],
                   lim = lim2, ln = T, resamp = T)
}
dev.off()

###############################################
## Plot reliability diagram for each weekday ##

## Aggregate mean forecasts
models_agg <- list()
for (i in 1:n_mods) models_agg[[i]] <- rowSums(models[[i]])
obs_agg <- rowSums(obs)

day_names <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
               "Saturday", "Sunday")
day_start <- c(3,4,5,6,7,1,2)

for (i in 1:7) {
  day_ind <- seq(day_start[i], n_days, by = 7)
  models_days <- lapply(models_agg, function(x) x[day_ind])
  obs_days <- obs_agg[day_ind]
  # Make plot
  # Reliability diagram on log scale
  file_path <- file.path(fpath, paste("mean_reldiag", day_names[i], "log.pdf",
                                      sep = "_"))
  lim2 <- c(0.08, 9)
  pdf(file_path, width = 7.5, height = 8)
  par(mfrow = c(2,2), mar = c(3.6, 3, 2.2, 0.3))
  for (j in 1:n_mods) {
    plot_reliability(models_days[[j]], obs_days, model_names[j],
                     model_colors[j], lim = lim2, ln = T, resamp = T)
  }
  dev.off()
}
