

getTI <- function(DD, MM, YY, times) {
  matchday <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY) 
  if (any(matchday)) {
    TI <- which(matchday)
  } else {
    TI <- NA
  }
  return(TI)
}

plotScores <- function(scores, times, mnames, mcols, filePath, events = NULL, logscale = T,
                       type = "l", days = NULL, ylim = NULL, tlim = NULL, whichmods = 1:4) {
  # Check arguments
  if (hasArg(tlim)) {
    tbeg <- getTI(tlim[[1]][1], tlim[[1]][2], tlim[[1]][3], times)
    tend <- getTI(tlim[[2]][1], tlim[[2]][2], tlim[[2]][3], times)
    stopifnot(tbeg < tend)
    scores <- as.matrix( scores[tbeg:tend, ] )
    times <- times[tbeg:tend, ]
    if (hasArg(events)) {
      events <- events[(events$TI >= tbeg) & (events$TI <= tend), ]
      events$TI <- events$TI - tbeg + 1
    }
    if (hasArg(days)) {
      days <- days[tbeg:tend]
      days <- days - tbeg + 1
    }
  }
  haslogscale <- hasArg(logscale)
  if ( !any(type == c("l", "h", "p")) ) {
    warning(paste0("Invalid plot type \'", type, "\'"))
    type <- "l"
  }
  if (missing(whichmods)) whichmods <- 1:(dim(scores)[2])
  if (missing(days)) {
    days <- 1:dim(times)[1]
  } else {
    stopifnot(length(days) == dim(times)[1])
  }
  anyneg <- any(scores[days, ] <= 0, na.rm = TRUE)
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
  if (hasArg(tlim) & ndays < 1000) {
    xdays <- (times$DD == 1)
    xlabels <- paste0(sprintf("%2d", times$MM[xdays]), "\'", times$YY[xdays] %% 100)
  } else {
    xdays <- (times$DD == 1) & (times$MM == 1)
    xlabels <- times$YY[xdays]
  }
  xats <- which(xdays)
  # Start plotting
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
  axis(1, at = xats, labels = xlabels, gap.axis = 0.9)
  # Add legend
  offsetleft <- -10
  leg <- legend(offsetleft, ylim[2], mnames, col = mcols, lwd = 2, plot = FALSE)$rect
  ind_leg <- days[1:(round(leg$w) + offsetleft)]
  overplot_leg <- sum(scores[ind_leg, ] > leg$top - leg$h, na.rm = TRUE)
  if (overplot_leg > 5) {
    legend(ndays - offsetleft - leg$w, ylim[2], mnames, col = mcols, lwd = 2, bg = "white")
  } else {
    legend(offsetleft, ylim[2], mnames, col = mcols, lwd = 2)
  }
  dev.off()
}


plotScoreDiffs <- function(scores, times, mnames, mcols, filePath, events = NULL, type = "l",
                           days = NULL, ylim = NULL, trim = NULL, tlim = NULL, whichmods = 1:4) {
  # Check arguments
  if (hasArg(tlim)) {
    tbeg <- getTI(tlim[[1]][1], tlim[[1]][2], tlim[[1]][3], times)
    tend <- getTI(tlim[[2]][1], tlim[[2]][2], tlim[[2]][3], times)
    stopifnot(tbeg < tend)
    scores <- as.matrix( scores[tbeg:tend, ] )
    times <- times[tbeg:tend, ]
    if (hasArg(events)) {
      events <- events[(events$TI >= tbeg) & (events$TI <= tend), ]
      events$TI <- events$TI - tbeg + 1
    }
    if (hasArg(days)) {
      days <- days[tbeg:tend]
      days <- days - tbeg + 1
    }
  }
  if (hasArg(trim)) scores <- pmax( pmin(scores, trim[2]), trim[1] )
  if ( !any(type == c("l", "p")) ) {
    warning(paste0("Invalid plot type \'", type, "\'"))
    type <- "l"
  }
  if (missing(whichmods)) whichmods <- 1:(dim(scores)[2])
  if (missing(days)) {
    days <- 1:dim(times)[1]
  } else {
    stopifnot(length(days) == dim(times)[1])
  }
  # Set parameters
  ndays <- length(days)
  nmods <- length(whichmods)
  mnames <- mnames[whichmods]
  mcols <- mcols[whichmods]
  ylab <- "score"
  # Determine ylim
  rs <- range(scores[days, ], na.rm = TRUE)
  if (missing(ylim) || ylim[2] < rs[1] || ylim[1] > rs[2]) {
    ylim <- rs
  } else {
    ylim <- c( max(rs[1], ylim[1]), min(rs[2], ylim[2]) )
  }
  # Ajust ylim if events are added to the plot
  if (hasArg(events)) {
    # Determine event circle coordinates
    evt_days <- unique(events$TI)
    evt_y <- ylim[1] - 3/96 * diff(ylim)
    ylim[1] <- ylim[1] - 4/96 * diff(ylim)
  }
  # Determine x-axis ticks and labels
  if (hasArg(tlim) & ndays < 1000) {
    xdays <- (times$DD == 1)
    xlabels <- paste0(sprintf("%2d", times$MM[xdays]), "\'", times$YY[xdays] %% 100)
  } else {
    xdays <- (times$DD == 1) & (times$MM == 1)
    xlabels <- times$YY[xdays]
  }
  xats <- which(xdays)
  # Start plotting
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
  abline(h = 0)
  # Add event circles if possible
  if (hasArg(events))  points(evt_days, rep(evt_y, length(evt_days)), cex = 0.8)
  axis(1, at = xats, labels = xlabels, gap.axis = 0.9)
  # Add legend
  offsetleft <- -10
  leg <- legend(offsetleft, ylim[2], mnames, col = mcols, lwd = 2, plot = FALSE)$rect
  ind_leg <- days[1:(round(leg$w) + offsetleft)]
  overplot_leg <- sum(scores[ind_leg, ] > leg$top - leg$h, na.rm = TRUE)
  if (overplot_leg > 5) {
    legend(ndays - offsetleft - leg$w, ylim[2], mnames, col = mcols, lwd = 2, bg = "white")
  } else {
    legend(offsetleft, ylim[2], mnames, col = mcols, lwd = 2)
  }
  dev.off()
}
