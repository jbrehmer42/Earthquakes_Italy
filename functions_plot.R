####################################
## Auxiliary functions - Plotting ##

# Load package for maps
library(maps)


get_time_index <- function(DD, MM, YY, times) {
  # Compute time index (i.e. data row) from given date
  #
  # Input values:
  # DD     -  day of date as integer
  # MM     -  month of date as integer
  # YY     -  year of date as integer
  # times  -  time stamps of the model outputs
  match_day <- times == ymd(paste(YY, MM, DD))
  if (any(match_day)) {
    time_index <- which(match_day)
  } else {
    time_index <- NA
  }
  return(time_index)
}


plot_reliability <- function(aggr, y, txt = "", col = "black", lim = NULL,
                             ln = F, resamp = NULL, MCB_decomp = F) {
  # Create a mean reliability diagram
  # Based on code by Johannes Resin.
  #
  # Input values:
  # aggr       - forecast vector
  # y          - observation vector
  # txt        - title for the diagram
  # col        - color of the isotonic regression curve
  # lim        - limits for both axes
  # ln         - using log to transform both axes
  # resamp     - number of samples if resampling should be done
  # MCB_decomp - Differentiate between conditional and
  #              unconditional miscalibration
  # Plots a mean reliability diagram similar
  # to the one in Gneiting and Resin (2021).
  # Resampling is done, but this is not as
  # straightforward, since all forecasts have
  # to be nonnegative.
  yf <- isoreg(aggr, y)$yf
  x <- sort(aggr)
  # calculate score values
  score <- function(x,y) mean((x-y)^2)
  s <- score(aggr, y)
  s_rc <- score(yf, y[order(aggr)])
  s_mg <- score(mean(y), y)
  # Decompose MCB into conditional and unconditional part
  if (MCB_decomp)  {
    res <- y - aggr
    c_rc_ucond = optim(par = 0, fn = function(c) score(aggr + c, y), 
                       method = "Brent", lower = min(res), upper = max(res))$par
    s_rc_ucond = score(aggr + c_rc_ucond, y)
    uMCB = s - s_rc_ucond
    cMCB = s_rc_ucond - s_rc
  }
  # calculate decomposition
  MCB <- s - s_rc
  DSC <- s_mg - s_rc
  UNC <- s_mg
  # Do resampling
  if (!missing(resamp)) {
    res <- y[order(aggr)] - x
    if (is.numeric(resamp)) n <- round(resamp) + 1 else n <- 5000
    lower_limit <- floor(n * (1 - 0.9)/2)
    #lower_limit <- floor(n * (1 - 0.7)/2)
    upper_limit <- n - lower_limit
    #upper_limit <- n - lower_limit
    # Resample residuals under assumption of iid residuals
    resamples <- sapply(1:n, function(z) x + sample(res, length(y)))
    # use rounding for integer data 
    resamples <- round(pmax(resamples, 0))
    yf_resamples <- apply(resamples, 2, function(z) isoreg(x, z)$yf )
    # sort resamples, include observed values, and correct bias
    # (shift by mean residual)
    yf_resamples_sorted <- apply(cbind(yf, yf_resamples), 1, sort) - mean(res)
    # Compute limits for plotting
    pind <- (x < yf_resamples_sorted[upper_limit, ]) &
              (x > yf_resamples_sorted[lower_limit, ])
    # Compute MCB p-value
    MCB_resamples <- sapply(1:n, function(i) score(x, resamples[ ,i])
                            - score(yf_resamples[ ,i], resamples[ ,i]))
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
            c(yf_resamples_sorted[upper_limit,pind],
              rev(yf_resamples_sorted[lower_limit,pind])),
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
  if (MCB_decomp) {
    text(x = lim[1] + offs, y = lim[2] - offs, adj = c(0,1),
         labels = paste0(c("uMCB", "cMCB", "DSC ", "UNC", "SCR"),
                         collapse = "\n"))
    text(x = lim[1] + 4 * offs, y = lim[2] - offs, adj = c(0,1),
         labels = paste0(round(c(uMCB, cMCB, DSC, UNC, s), digits = 3),
                         collapse = "\n"))
  } else {
    text(x = lim[1] + offs, y = lim[2] - offs, adj = c(0,1),
         labels = paste0(c("MCB ","DSC ","UNC", "SCR"), collapse = "\n"))
    text(x = lim[1] + 4 * offs, y = lim[2] - offs, adj = c(0,1),
         labels = paste0(round(c(MCB, DSC, UNC, s), digits = 3),
                         collapse = "\n"))
  }
  if (!missing(resamp)) {
    text(x = lim[1] + 7 * offs, y = lim[2] - offs, adj = c(0,1),
         paste0("[p = ", round(pval, 3), "]"))
  }
  return(invisible(yf))
}


plot_map <- function(color_vals, cells, main = "", evts = NULL, borders = T) {
  # Plot map with colored cells
  #
  # Input values:
  # color_vals  - color values for the cells given by cells
  # cells       - data frame with longitudes and latitudes of the cells
  #               which should be colored
  # main        - title for the map
  # evts        - data frame of events. Added to the plot as small 
  #               diamonds (optional)
  # borders     - whether to add state borders to the map
  # This function can be used to illustrate the spatial behavior of 
  # different forecasts, e.g. calibration. The values for each cell
  # have to be transformed to color values first.
  border_col <- rgb(0, 0, 0, alpha = 0.4)
  xlim <- c(min(cells$LON), max(cells$LON))
  ylim <- c(min(cells$LAT), max(cells$LAT))
  plot(1, 1, xlim = xlim, ylim = ylim, col = "white", asp = 1.3, xaxt = "n",
       yaxt = "n", xlab = "", ylab = "", main = main)
  points(cells$LON, cells$LAT, pch = 15, col = color_vals, cex = 0.4)
  if (borders) map('world', fill = F, add = T, col = border_col)
  if (!is.null(evts)) {
    cell_num <- unique(evts$N)
    plot_indic <- apply(as.matrix(cells$N), 1, function(x) any(x == cell_num))
    points(cells$LON[plot_indic], cells$LAT[plot_indic], pch = 5,
           col = border_col, lwd = 0.5, cex = 0.4)
  }
}


plot_multi_map <- function(vals, pal, cells, lims, n_colors, file_path,
                           offset = 0, evts = NULL) {
  # Plot maps for four models in one graph
  # 
  # Input values:
  # vals       -  values which are transformed to colors
  # pal        -  collection of colors to use for the maps
  # cells      -  data frame with longitudes and latitudes of the
  #               cells which should be colored
  # lims       -  limits to cut off very large or small values
  # n_colors   -  number of colors to use for mapping the values
  #               to a color scale (has to be length(pal), obsolete)
  # file_path  -  name of .pdf file where the graph is saved
  # offset     -  small number to shift the values such that log
  #               does not attain extremely negative values
  # evts       -  data frame of events. Added to the plots as small 
  #               diamonds (optional)
  # Uses plot_map to plot the values given by vals (e.g.
  # calibration) in color on four different maps. Adds a color
  # bar on the right to indicate large and small values.
  y_len <- n_colors/10
  n_ticks <- 7
  n_mods <- dim(vals)[2]
  pdf(file_path, width = 8, height = 7)
  layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE),
         widths = c(0.42, 0.42, 0.16))
  ## start with maps
  par(mar = c(2/3,2/3,1.5,1/3))
  for (i in 1:n_mods) {
    scaled_vals <- (log(vals[ ,i] + offset) - lims[1]) / (lims[2] - lims[1])
    color_scale <- pal[round(scaled_vals * (n_colors-1)) + 1]
    plot_map(color_scale, cells, main = model_names[i], evts = evts)
  }
  ## add color bar
  par(mar = c(4,1,1,4), mgp = c(3,0,-1))
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(1,y_len), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labels <- sprintf("%.e", exp(seq(lims[1], lims[2], len = n_ticks)) - offset)
  if (round(exp(lims[1]) - offset, 15) == 0) labels[1] <- "0"
  axis(4, at = seq(1, y_len, len = n_ticks), labels = labels, las = 1)
  for (i in 1:(n_colors-1)) {
    rect(0, 1 + (i-1) * (y_len - 1)/(n_colors - 1), 1, 1 +
           i * (y_len - 1)/(n_colors - 1), 
         col = pal[i], border = NA)
  }
  dev.off()
}


plot_skill_map <- function(vals, pal, cells, lims, n_colors, file_path,
                           evts = NULL) {
  # Plot skill maps for four models in one graph.
  # Similar to plot_multi_map, but with different treatment of
  # values and colors and an asymmetric color bar.
  #
  # Input values:
  # vals       -  skill values which are transformed to colors
  # pal        -  collection of colors to use for the maps
  # cells      -  data frame with longitudes and latitudes of the
  #               cells which should be colored
  # lims       -  limits to cut off very large or small values
  # n_colors   -  number of colors to use for mapping the values
  #               to a color scale (has to be length(pal), obsolete)
  # file_path  -  name of .pdf file where the graph is saved
  # evts       -  data frame of events. Added to the plots as small 
  #               diamonds (optional)
  # Uses plot_map to plot the skill valueson four different maps.
  # Adds a color bar on the right to indicate positive and negative
  # values.
  n_mods <- dim(vals)[2]
  n_cells <- dim(cells)[1]
  pdf(file_path, width = 8, height = 7)
  layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE),
         widths = c(0.42, 0.42, 0.16))
  ## start with maps
  par(mar = c(2/3,2/3,1.5,1/3))
  for (i in 1:n_mods) {
    skills <- pmax(vals[ ,i], lims[1])
    # Use different scaling and colors for positive and negative
    # values because skills are in (-\infty , 1]
    scaled_vals_pos <-  skills[skills >= 0]
    scaled_vals_neg <-  - (skills[skills < 0] - lims[1]) / lims[1]
    color_scale <- rep("", n_cells)
    color_scale[skills >= 0] <- hsv(pal[1], scaled_vals_pos,
                                    v = 1 - 0.4 * scaled_vals_pos)
    color_scale[skills < 0]  <- hsv(pal[2], 1 - scaled_vals_neg)
    plot_map(color_scale, cells, main = model_names[i], evts = evts)
  }
  ## add color bar
  par(mar = c(4,1,1,4), mgp = c(3,0,-1))
  len_pos <- (n_colors/10 - 1)/3
  len_neg <- 2 * len_pos
  n_ticks <- 3
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(-len_neg, len_pos), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labels <- round(seq(0, lims[2], len = n_ticks+1), 2)
  labels <- c(round(seq(lims[1], 0, len = 2*n_ticks), 2),
              labels[2:length(labels)])
  ats <- seq(0, len_pos, len = n_ticks+1)
  ats <- c(seq(-len_neg, 0, len = 2*n_ticks), ats[2:length(ats)])
  labels[1] <- paste0("< ", lims[1])
  axis(4, at = ats, labels = labels, las = 1)
  # Set number of colors for color bar
  kk <- 50
  # positive part
  for (l in 1:(kk-1)) {
    rect(0, (l-1) * len_pos/(kk - 1), 1,  l * len_pos/(kk - 1), 
         col = hsv(pal[1], l/(kk-1), v = 1 - 0.4 * l/(kk-1)), border = NA)
  }
  kk <- 2*kk
  # negative part
  for (l in 1:(kk-1)) {
    rect(0, - (l-1) * len_neg/(kk - 1), 1,  - l * len_neg/(kk - 1), 
         col = hsv(pal[2], l/(kk-1)), border = NA)
  }
  dev.off()
}


plot_diffs_map <- function(vals, pal, cells, file_path, evts = NULL,
                           borders = T, main = "") {
  # Plot map with colored cells and color bar on the right.
  # Like plot_map, but designed to illustrate value differences.
  #
  # Input values:
  # vals       -  values which are transformed to colors
  # pal        -  collection of colors to use for the maps
  # cells      -  data frame with longitudes and latitudes of the
  #               cells which should be colored
  # file_path  -  name of .pdf file where the graph is saved
  # evts       -  data frame of events. Added to the plots as small 
  #               diamonds (optional)
  # borders    -  whether to add state borders to the map
  # main       -  title of the graph
  # Set graphical paramters e.g. colors
  y_len <- 5
  n_ticks <- 7
  # Transform values in vals into colors: Scale them to values in the
  # interval [-1, 1] and interpret these values as saturation
  abs_vals <- 1.01 * max(abs(vals))
  scaled_vals <- (vals + abs_vals) / (2 * abs_vals)
  scaled_vals <- 2 * scaled_vals - 1
  color_scale <- rep("", length(vals))
  color_scale[scaled_vals >= 0] <- hsv(pal[1], scaled_vals[scaled_vals >= 0])
  color_scale[scaled_vals < 0]  <- hsv(pal[2], -scaled_vals[scaled_vals < 0])
  # Create pdf file with two parts
  pdf(file_path, width = 5, height = 4.4)
  layout(matrix(1:2, 1, 2, byrow = T), widths = c(0.78, 0.22))
  # Plot the map
  if (main == "") {
    par(mar = 2/3 * rep(1,4))
  } else {
    par(mar = c(2/3, 2/3, 1.5, 2/3))
  }
  plot_map(color_scale, cells, main = main, evts = evts, borders = borders)
  # Add color bar to the right of the map
  par(mar = c(1/3,0,1/3,2), mgp = c(-1,-2.4,-3.4), cex = 0.6)
  kk <- y_len * 20
  plot(1,1, col = "white", xlim = c(0,1), ylim = c(-y_len, y_len), asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  labels <- sprintf("%.1e", seq(- abs_vals, abs_vals, len = n_ticks))
  labels[(n_ticks + 1)/2] <- "0"
  axis(4, at = seq(-y_len, y_len, len = n_ticks), labels = labels, las = 1)
  for (l in 1:(kk-1)) {
    sat <- l/(kk-1)
    # positive part
    rect(0, (l-1) * y_len/(kk - 1), 1/2,  l * y_len/(kk - 1), 
         col = hsv(pal[1], sat), border = NA)
    # negative part
    rect(0, - (l-1) * y_len/(kk - 1), 1/2,  -l * y_len/(kk - 1), 
         col = hsv(pal[2], sat), border = NA)
  }
  dev.off()
}



plot_scores <- function(scores, times, model_names, model_colors, file_path,
                        events = NULL, logscale = T, type = "l", days = NULL,
                        ylim = NULL, tlim = NULL, whichmods = 1:4) {
  # Plot daily scores of the forecast models
  #
  # Input values:
  # scores        -  list of scores (or other values) of the models
  # times         -  time stamps of the model forecasts
  # model_names   -  names of the models, depicted in the legend
  # model_colors  -  colors of the lines/points of the models
  # file_path     -  name of .pdf file where the graph is saved
  # events        -  data frame of events. Added to the plot as small
  #                  circles below the time axis (optional)
  # logscale      -  Whether to transform the y axis with log
  # type          -  How to plot the score values. Possible are:
  #                   "l" for lines
  #                   "p" for points
  #                   "h" for vertical bars
  # days          -  Increasing vector of integers for which to plot
  #                  the values. Set a position to NA to omit the 
  #                  value on this day (optional)
  # ylim          -  limits for the y axis (optional)
  # tlim          -  limits for the time axis (optional)
  # whichmods     -  indices of model forecasts in scores. E.g. if
  #                  only model 1 and 3 are given, set to c(1, 3)
  # Allows for a lot of (maybe too much) customization of the plot of
  # realized scores, but is not 'finished'.
  # Check arguments
  if (hasArg(tlim)) {
    t_start <- get_time_index(tlim[[1]][1], tlim[[1]][2], tlim[[1]][3], times)
    t_end <- get_time_index(tlim[[2]][1], tlim[[2]][2], tlim[[2]][3], times)
    stopifnot(t_start < t_end)
    scores <- as.matrix( scores[t_start:t_end, ] )
    times <- times[t_start:t_end]
    if (hasArg(events)) {
      events <- events[(events$TI >= t_start) & (events$TI <= t_end), ]
      events$TI <- events$TI - t_start + 1
    }
    if (hasArg(days)) {
      days <- days[t_start:t_end]
      days <- days - t_start + 1
    }
  }
  has_logscale <- hasArg(logscale)
  if ( !any(type == c("l", "h", "p")) ) {
    warning(paste0("Invalid plot type \'", type, "\'"))
    type <- "l"
  }
  if (missing(whichmods)) whichmods <- 1:(dim(scores)[2])
  if (missing(days)) {
    days <- 1:length(times)
  } else {
    stopifnot(length(days) == length(times))
  }
  negative_scores <- any(scores[days, ] <= 0, na.rm = TRUE)
  # Set parameters
  n_days <- length(days)
  n_mods <- length(whichmods)
  model_names <- model_names[whichmods]
  model_colors <- model_colors[whichmods]
  ylab <- "score"
  # Checks if logarithms should be plotted
  if(logscale & negative_scores) {
    logscale <- FALSE
    if (has_logscale) {
      warning("Cannot use logarithmic scale. There are negative values")
    }
  }
  if (logscale) {
    scores <- log(scores)
    ylab <- paste(ylab, "(log scale)")
  }
  # Determine ylim
  range_scores <- range(scores[days, ], na.rm = TRUE)
  if (missing(ylim) || ylim[2] < range_scores[1] || ylim[1] > range_scores[2]) {
    ylim <- range_scores
  } else {
    ylim <- c( max(range_scores[1], ylim[1]), min(range_scores[2], ylim[2]) )
  }
  ylow <- range_scores[1]
  # Ajust ylim if events are added to the plot
  if (hasArg(events)) {
    # Determine event circle coordinates
    event_days <- unique(events$TI)
    event_y <- ylim[1] - 3/96 * diff(ylim)
    ylim[1] <- ylim[1] - 4/96 * diff(ylim)
  }
  # Determine x-axis ticks and labels
  if (hasArg(tlim) & n_days < 1000) {
    x_days <- (day(times) == 1)
    x_labels <- paste0(sprintf("%2d", month(times[x_days])), "\'",
                       year(times[x_days]) %% 100)
  } else {
    x_days <- (day(times) == 1) & (month(times) == 1)
    x_labels <- year(times[x_days])
  }
  x_ats <- which(x_days)
  # Start plotting
  pdf(file_path, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:n_days, rep(0, n_days), ylim = ylim, xlab = "days", col = "white",
       ylab = ylab, main = "", xaxt="n")
  if (type == "l") {
    for (i in 1:n_mods) lines(days, scores[ ,i], col = model_colors[i])
  }
  if (type == "p") {
    for (i in 1:n_mods) {
      points(days, scores[ ,i], col = model_colors[i], cex = 1/2)
    }
  }
  if (type == "h") {
    for (j in days) {
      if (is.na(j)) next
      vals <- sort(scores[j, ], decreasing = T)
      col <- model_colors[order(scores[j, ], decreasing = T)]
      for (i in 1:n_mods) {
        lines(c(j,j), c(ylow, vals[i]), col = col[i], lwd = 0.3)
      }
    }
  }
  # Add event circles if possible
  if (hasArg(events)) {
    points(event_days, rep(event_y, length(event_days)), cex = 0.8)
  }
  axis(1, at = x_ats, labels = x_labels, gap.axis = 0.9)
  # Add legend
  offset_left <- -10
  legend_box <- legend(offset_left, ylim[2], model_names, col = model_colors,
                       lwd = 2, plot = FALSE)$rect
  ind_leg <- days[1:(round(legend_box$w) + offset_left)]
  # Check how many values are blocked by legend
  sum_vals_blocked <- sum(scores[ind_leg, ] > legend_box$top - legend_box$h,
                      na.rm = TRUE)
  # If more than 5 values are blocked, shift legend to the right
  if (sum_vals_blocked > 5) {
    legend(n_days - offset_left - legend_box$w, ylim[2], model_names,
           col = model_colors, lwd = 2, bg = "white")
  } else {
    legend(offset_left, ylim[2], model_names, col = model_colors, lwd = 2)
  }
  dev.off()
}


plot_score_diffs <- function(scores, times, model_names, model_colors,
                             file_path, events = NULL, type = "l", days = NULL,
                             ylim = NULL, trim = NULL, tlim = NULL,
                             whichmods = 1:4) {
  # Plot daily score differences of the forecast models
  # Very similar to plot_scores
  #
  # Input values:
  # scores        -  list of score differences of the models
  # times         -  time stamps of the model forecasts
  # model_names   -  names of the models, depicted in the legend
  # model_colors  -  colors of the lines/points of the models
  # file_path     -  name of .pdf file where the graph is saved
  # events        -  data frame of events. Added to the plot as small
  #                  circles below the time axis (optional)
  # type          -  How to plot the score differences Possible are:
  #                   "l" for lines
  #                   "p" for points
  # days          -  Increasing vector of integers for which to plot
  #                  the values. Set a position to NA to omit the 
  #                  value on this day (optional)
  # ylim          -  limits for the y axis (optional)
  # trim          -  limits to truncate the score differences (optional)
  # tlim          -  limits for the time axis (optional)
  # whichmods     -  indices of model forecasts in scores. E.g. if
  #                  only model 1 and 3 are given, set to c(1, 3)
  # Check arguments
  if (hasArg(tlim)) {
    t_start <- get_time_index(tlim[[1]][1], tlim[[1]][2], tlim[[1]][3], times)
    t_end <- get_time_index(tlim[[2]][1], tlim[[2]][2], tlim[[2]][3], times)
    stopifnot(t_start < t_end)
    scores <- as.matrix( scores[t_start:t_end, ] )
    times <- times[t_start:t_end]
    if (hasArg(events)) {
      events <- events[(events$TI >= t_start) & (events$TI <= t_end), ]
      events$TI <- events$TI - t_start + 1
    }
    if (hasArg(days)) {
      days <- days[t_start:t_end]
      days <- days - t_start + 1
    }
  }
  if (hasArg(trim)) scores <- pmax( pmin(scores, trim[2]), trim[1] )
  if ( !any(type == c("l", "p")) ) {
    warning(paste0("Invalid plot type \'", type, "\'"))
    type <- "l"
  }
  if (missing(whichmods)) whichmods <- 1:(dim(scores)[2])
  if (missing(days)) {
    days <- 1:length(times)
  } else {
    stopifnot(length(days) == length(times))
  }
  # Set parameters
  n_days <- length(days)
  n_mods <- length(whichmods)
  model_names <- model_names[whichmods]
  model_colors <- model_colors[whichmods]
  ylab <- "score"
  # Determine ylim
  range_scores <- range(scores[days, ], na.rm = TRUE)
  if (missing(ylim) || ylim[2] < range_scores[1] || ylim[1] > range_scores[2]) {
    ylim <- range_scores
  } else {
    ylim <- c( max(range_scores[1], ylim[1]), min(range_scores[2], ylim[2]) )
  }
  # Ajust ylim if events are added to the plot
  if (hasArg(events)) {
    # Determine event circle coordinates
    event_days <- unique(events$TI)
    event_y <- ylim[1] - 3/96 * diff(ylim)
    ylim[1] <- ylim[1] - 4/96 * diff(ylim)
  }
  # Determine x-axis ticks and labels
  if (hasArg(tlim) & n_days < 1000) {
    x_days <- (day(times) == 1)
    x_labels <- paste0(sprintf("%2d", month(times[x_days])), "\'",
                       year(times[x_days]) %% 100)
  } else {
    x_days <- (day(times) == 1) & (month(times) == 1)
    x_labels <- year(times[x_days])
  }
  x_ats <- which(x_days)
  # Start plotting
  pdf(file_path, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:n_days, rep(0, n_days), ylim = ylim, xlab = "days", col = "white",
       ylab = ylab, main = "", xaxt="n")
  if (type == "l") {
    for (i in 1:n_mods) lines(days, scores[ ,i], col = model_colors[i])
  }
  if (type == "p") {
    for (i in 1:n_mods) {
      points(days, scores[ ,i], col = model_colors[i], cex = 1/2)
    }
  }
  abline(h = 0)
  # Add event circles if possible
  if (hasArg(events)) {
    points(event_days, rep(event_y, length(event_days)), cex = 0.8)
  }
  axis(1, at = x_ats, labels = x_labels, gap.axis = 0.9)
  # Add legend
  offset_left <- -10
  legend_box <- legend(offset_left, ylim[2], model_names, col = model_colors,
                       lwd = 2, plot = FALSE)$rect
  ind_leg <- days[1:(round(legend_box$w) + offset_left)]
  # Check how many values are blocked by legend
  sum_vals_blocked <- sum(scores[ind_leg, ] > legend_box$top - legend_box$h,
                          na.rm = TRUE)
  # If more than 5 values are blocked, shift legend to the right
  if (sum_vals_blocked > 5) {
    legend(n_days - offset_left - legend_box$w, ylim[2], model_names,
           col = model_colors, lwd = 2, bg = "white")
  } else {
    legend(offset_left, ylim[2], model_names, col = model_colors, lwd = 2)
  }
  dev.off()
}


plot_elementary_scores <- function(vals, grd, model_names, model_colors,
                                   file_path, ylab,
                                   model_ltys = rep(1, max(whichmods)),
                                   whichmods = 1:4) {
  # Plot elemenatry scores, MCB etc. for each theta
  #
  # Input values:
  # vals          -  list of values of the models
  # grd           -  grid of theta values for which values are
  #                  given
  # model_names   -  names of the models, depicted in the legend
  # model_colors  -  colors of the lines of the models
  # file_path     -  name of .pdf file where the graph is saved
  # ylab          -  title for y axis
  # model_ltys    -  line types of the models (optional)
  # whichmods     -  indices of model forecasts in scores. E.g. if
  #                  only model 1 and 3 are given, set to c(1, 3)
  n_theta <- length(grd)
  log_grid <- log(grd)
  n_mods <- length(whichmods)
  model_names <- model_names[whichmods]
  model_colors <- model_colors[whichmods]
  model_ltys <- model_ltys[whichmods]
  ylim <- c(min(vals), max(vals))
  # Start plotting
  pdf(file_path, width = 8, height = 5.5)
  par(mar = c(4, 4, 0.5, 0.5))
  plot(1:n_theta, 1:n_theta, ylim = ylim, xlab = "log(theta)", ylab = ylab,
       xaxt = "n", col = "transparent")
  for (i in 1:n_mods) {
    lines(1:n_theta, vals[ ,i], col = model_colors[i], lty = model_ltys[i])
  }
  # create log axis
  ticks <- axis(1, labels = F, tick = F)
  labels <- round(log_grid[pmax(1,ticks)], 1)
  axis(1, at = ticks, labels = labels)
  legend(1, ylim[2], model_names, col = model_colors, lty = model_ltys, lwd = 2)
  dev.off()
}
