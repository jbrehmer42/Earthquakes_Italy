

plotScores <- function(scores, times, mnames, mcols, filePath, events = NULL, logscale = T,
                       type = "l", days = NULL, ylim = NULL, tlim = NULL, whichmods = 1:4) {
  # Check arguments
  haslogscale <- hasArg(logscale)
  if ( !any(type == c("l", "h", "p")) ) {
    warning(paste0("Invalid plot type \'", type, "\'"))
    type <- "l"
  }
  if (missing(whichmods)) whichmods <- 1:(dim(as.matrix(scores))[2])
  if (missing(days)) {
    days <- 1:dim(times)[1]
  } else {
    stopifnot(length(days) == dim(times)[1])
  }
  anyneg <- any(scores <= 0, na.rm = TRUE)
  # Set parameters
  ndays <- length(days)
  nmods <- length(whichmods)
  mnames <- mnames[whichmods]
  mcols <- mcols[whichmods]
  ylab <- "score"
  # Checks if logarithms should be plotted
  if(logscale & anyneg) {
    logscale <- FALSE
    if (haslogscale) warning("Cannot use logarithmic scale. There are negative values")
  }
  if (logscale) {
    scores <- log(scores)
    ylab <- paste(ylab, "(log scale)")
  }
  # Determine ylim
  rs <- range(scores[days, ], na.rm = TRUE)
  if (missing(ylim) || ylim[2] < rs[1] || ylim[1] > rs[2]) {
    ylim <- rs
  } else {
    ylim <- c( max(rs[1], ylim[1]), min(rs[2], ylim[2]) )
  }
  ylow <- rs[1]
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
  if (type == "l") {
    for (i in 1:nmods) lines(days, scores[ ,i], col = mcols[i])
  }
  if (type == "p") {
    for (i in 1:nmods) points(days, scores[ ,i], col = mcols[i], cex = 1/2)
  }
  if (type == "h") {
    for (j in days) {
      if (is.na(j)) next
      vals <- sort(scores[j, ], decreasing = T)
      col <- mcols[order(scores[j, ], decreasing = T)]
      for (i in 1:nmods) lines(c(j,j), c(ylow,vals[i]), col = col[i], lwd = 0.3)
    }
  }
  # Add event circles if possible
  if (hasArg(events))  points(evt_days, rep(evt_y, length(evt_days)), cex = 0.8)
  axis(1, at = xats, labels = xlabels)
  # Add legend
  legend(-10, ylim[2], mnames, col = mcols, lwd = 2)
  dev.off()
}
