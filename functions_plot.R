## Auxilary functions for plotting of
## Italian earthquake data

# Load package for maps
library(maps)

# function to create a reliability diagram
simple_diag <- function(aggr, y, txt = "", col = "black", lim = NULL, ln = F, resamp = NULL) {
  # Plots a mean reliability diagram similar
  # to the one in Gneiting and Resin (2021)
  # Resampling is done, but this is not as
  # straightforward, since all forecasts have
  # to be nonnegative!
  yf <- isoreg(aggr, y)$yf
  x <- sort(aggr)
  # calculate decomposition
  score <- function(x,y) mean((x-y)^2)
  s <- score(aggr, y)
  s_rc <- score(yf, y[order(aggr)])
  s_mg <- score(mean(y), y)
  mcb <- s - s_rc
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
       labels = paste0(c("MCB ","DSC ","UNC", "SCR"), collapse = "\n"))
  text(x = lim[1] + 4 * offs, y = lim[2] - offs, adj = c(0,1),
       labels = paste0(round(c(mcb, dsc, unc, s), digits = 3), collapse = "\n"))
  if (!missing(resamp)) {
    text(x = lim[1] + 7 * offs, y = lim[2] - offs, adj = c(0,1),
         paste0("[p = ", round(pval, 3), "]"))
  }
  return(invisible(yf))
}



# xlim <- c(min(bins$LON), max(bins$LON))
# ylim <- c(min(bins$LAT), max(bins$LAT))
plot_map <- function(vals, main = "", evts = NULL, borders = T) {
  border_col <- rgb(0, 0, 0, alpha = 0.4) 
  plot(1, 1, xlim = xlim, ylim = ylim, col = "white", asp = 1.3, xaxt = "n",
       yaxt = "n", xlab = "", ylab = "", main = main)
  points(bins$LON, bins$LAT, pch = 15, col = vals, cex = 0.4)
  if (borders) map('world', fill = F, add = T, col = border_col)
  if (!is.null(evts)) {
    binN <- unique(evts$N)
    ind <- apply(as.matrix(bins$N), 1, function(x) any(x == binN))
    points(bins$LON[ind], bins$LAT[ind], pch = 5, col = border_col, 
           lwd = 0.5, cex = 0.4)
  }
}


plot_figure <- function(vals, pal, lims, ncols, file, offset = 0, evts = F) {
  if (evts) evts <- M4events else evts <- NULL
  ylen <- ncols/10
  nticks <- 7
  filePath <- paste0(file, ".pdf")
  pdf(filePath, width = 8, height = 7)
  layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE), widths = c(0.42, 0.42, 0.16))
  ## start with maps
  par(mar = c(2/3,2/3,1.5,1/3))
  for (i in 1:4) {
    scl <- (log(vals[ ,i] + offset) - lims[1]) / (lims[2] - lims[1])
    col_scl <- pal[round(scl * (ncols-1)) + 1]
    # main <- paste0("Model ", i)
    main <- mnames[i]
    plot_map(col_scl, main = main, evts = evts)
  }
  ## add color bar
  par(mar = c(4,1,1,4), mgp = c(3,0,-1))
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(1,ylen), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labs <- sprintf("%.e", exp(seq(lims[1], lims[2], len = nticks)) - offset)
  if (round(exp(lims[1]) - offset, 15) == 0) labs[1] <- "0"
  axis(4, at = seq(1, ylen, len = nticks), labels = labs, las = 1)
  for (i in 1:(ncols-1)) {
    rect(0, 1 + (i-1) * (ylen - 1)/(ncols - 1), 1, 1 + i * (ylen - 1)/(ncols - 1), 
         col = pal[i], border = NA)
  }
  dev.off()
}


plot_diff_maps <- function(vals, mdiff, scr_abs, file, evts = F) {
  if (evts) evts <- M4events else evts <- NULL
  m1 <- mdiff[1]
  m2 <- mdiff[2]
  ylen <- 5
  nticks <- 7
  main <- paste0("M", m1, " - ", "M", m2, "  (", mnames[m1], " vs. ", mnames[m2], ")")
  neg_col <- 0.66       # spec via hue in hsv colors
  pos_col <- 0          # spec via hue in hsv colors
  scl <- (vals[ ,m1] - vals[ ,m2] + scr_abs) / (2 * scr_abs)
  scl <- 2 * scl - 1
  col_scl <- rep("", nbins)
  col_scl[scl >= 0] <- hsv(pos_col, scl[scl >= 0])
  col_scl[scl < 0]  <- hsv(neg_col, -scl[scl < 0])
  ##
  filePath <- paste0(file, ".pdf")
  pdf(filePath, width = 5, height = 14/3)
  layout(matrix(1:2, 1, 2, byrow = TRUE), widths = c(0.78, 0.22))
  par(mar = c(2/3,2/3,2,2/3))
  plot_map(col_scl, main = main, evts = evts)
  ## add color bar
  par(mar = c(1/3,0,2.3,2), mgp = c(-1,-2.4,-3.4), cex = 0.6)
  kk <- ylen * 20
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(-ylen, ylen), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labs <- sprintf("%.1e", seq(- scr_abs, scr_abs, len = nticks))
  labs[(nticks + 1)/2] <- "0"
  axis(4, at = seq(-ylen, ylen, len = nticks), labels = labs, las = 1)
  for (l in 1:(kk-1)) {
    sat <- l/(kk-1)
    rect(0, (l-1) * ylen/(kk - 1), 1/2,  l * ylen/(kk - 1), 
         col = hsv(pos_col, sat), border = NA)
    rect(0, - (l-1) * ylen/(kk - 1), 1/2,  -l * ylen/(kk - 1), 
         col = hsv(neg_col, sat), border = NA)
  }
  dev.off()
}


# Plot daily scores of the models
plotScores <- function(scores, times, mnames, mcols, filePath, events = NULL, logscale = T) {
  # Input values:
  # scores   - Matrix of daily scores
  # times    - Time stamps of the scores
  # mnames   - Model names for the legend
  # mcols    - Array of colors for the lines
  # filePath - File path for .pdf file
  # events   - Data frame of events (optional)
  # logscale - Should scores be on log scale?
  ndays <- dim(times)[1]
  ylab <- "score"
  if(logscale) {
    scores <- log(scores)
    ylab <- paste(ylab, "(log scale)")
  }
  # Determine ylim
  ylim <- range(scores)
  # Ajust ylim if events are added to the plot
  if (hasArg(events)) {
    # Determine event circle coordinates
    evt_days <- unique(events$TI)
    evt_y <- ylim[1] - 3/96 * diff(ylim)
    ylim[1] <- ylim[1] - 4/96 * diff(ylim)
  }
  # Determine x-axis ticks and labels
  january1s <- (times$DD == 1) & (times$MM == 1)
  xlabels <- times$YY[january1s]
  xats <- which(january1s)
  pdf(filePath, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:ndays, rep(0, ndays), ylim = ylim, xlab = "days", col = "white",
       ylab = ylab, main = "", xaxt="n")
  for (i in 1:ncol(scores)) {
    lines(1:ndays, scores[ ,i], col = mcols[i])
  }
  # Add event circles if possible
  if (hasArg(events))  points(evt_days, rep(evt_y, length(evt_days)), cex = 0.8)
  axis(1, at = xats, labels = xlabels)
  legend(-10, ylim[2], mnames, col = mcols, lwd = 2)
  dev.off()
}


# Plot scores, MCB, etc. with respect to elementary scoring functions
plotElementary <- function(vals, grd, mnames, mcols, filePath, ylab, mltys = NULL,
                           whichmods = 1:4) {
  ntheta <- length(grd)
  lgrd <- log(grd)
  nmods <- length(whichmods)
  if (missing(mltys)) mltys <- rep(1, nmods)
  mnames <- mnames[whichmods]
  mcols <- mcols[whichmods]
  mltys <- mltys[whichmods]
  ylim <- c(min(vals), max(vals))
  # Start plotting
  pdf(filePath, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:ntheta, 1:ntheta, ylim = ylim, xlab = "log(theta)", ylab = ylab,
       xaxt = "n", col = "transparent")
  for (i in 1:nmods) {
    lines(1:ntheta, vals[i, ], col = mcols[i], lty = mltys[i])
  }
  # create log axis
  ticks <- axis(1, labels = F, tick = F)
  labs <- round(lgrd[pmax(1,ticks)], 1)
  axis(1, at = ticks, labels = labs)
  legend(2, ylim[2], mnames, col = mcols, lty = mltys, lwd = 2)
  dev.off()
}

# Plot a map of score differences by individually
# coloring the grid cells
diffMap <- function(vals, cells, filePath, events = NULL, borders = T) {
  # Input values:
  # vals     - Numeric values corresponding to
  #            the grid cells
  # cells    - Data frame of grid cells
  # filePath - File path for .pdf file
  # events   - Data frame of events (optional)
  # borders  - Plot national borders? Requires 
  #            R-package "maps"
  # Set graphical paramters e.g. colors
  ylen <- 5
  nticks <- 7
  neg_col <- 0.66       # spec via hue in hsv colors
  pos_col <- 0          # spec via hue in hsv colors
  border_col <- rgb(0, 0, 0, alpha = 0.4) 
  # Transform values in vals into colors: Scale them to values in the
  # interval [-1, 1] and interpret these values as saturation
  val_abs <- 1.01 * max(abs(vals))
  scl <- (vals + val_abs) / (2 * val_abs)
  scl <- 2 * scl - 1
  col_scl <- rep("", length(vals))
  col_scl[scl >= 0] <- hsv(pos_col, scl[scl >= 0])
  col_scl[scl < 0]  <- hsv(neg_col, -scl[scl < 0])
  # Create pdf file with two parts
  pdf(filePath, width = 5, height = 4.4)
  layout(matrix(1:2, 1, 2, byrow = T), widths = c(0.78, 0.22))
  # Plot the map
  par(mar = 2/3 * rep(1,4))
  plot(1, 1, xlim = range(cells$LON), ylim = range(cells$LAT), col = "white",
       asp = 1.3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "")
  points(cells$LON, cells$LAT, pch = 15, col = col_scl, cex = 0.4)
  # Add national borders if wanted
  if (borders) map("world", fill = F, add = T, col = border_col)
  # Add locations of events if possible
  if (!is.null(events)) {
    binN <- unique(events$N)
    ind <- apply(as.matrix(cells$N), 1, function(x) any(x == binN))
    points(cells$LON[ind], cells$LAT[ind], pch = 5, col = border_col, 
           lwd = 0.5, cex = 0.4)
  }
  # Add color bar to the right of the map
  par(mar = c(1/3,0,1/3,2), mgp = c(-1,-2.4,-3.4), cex = 0.6)
  kk <- ylen * 20
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(-ylen, ylen), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labs <- sprintf("%.1e", seq(- val_abs, val_abs, len = nticks))
  labs[(nticks + 1)/2] <- "0"
  axis(4, at = seq(-ylen, ylen, len = nticks), labels = labs, las = 1)
  for (l in 1:(kk-1)) {
    sat <- l/(kk-1)
    rect(0, (l-1) * ylen/(kk - 1), 1/2,  l * ylen/(kk - 1), 
         col = hsv(pos_col, sat), border = NA)
    rect(0, - (l-1) * ylen/(kk - 1), 1/2,  -l * ylen/(kk - 1), 
         col = hsv(neg_col, sat), border = NA)
  }
  dev.off()
}

