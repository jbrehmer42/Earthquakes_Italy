###################################
## Auxiliary functions - General ##

# Poisson scoring function
# Caution: Have to ensure that reports x are never zero
Spois <- function(x, y)  -y * log(x) + x

# Quadratic scoring function
Squad <- function(x,y) (y - x)^2

# Mean Poisson scoring function
Spois2 <- function(x, y) {
  # This is the same as Spois but omitting
  # the if-condition is a speed advantage
  s <- -y * log(x) + x
  s[is.nan(s)] <- 0
  return(mean(s))
}

# Quadratic interval restricted scoring function
Sbin <- function(x, y, a, b) {
  # This function is based on a simplified formula. For any [a,b]
  # the restricted quadratic score for this interval (see Taggart, 2020)
  # is given by S(x,y) = ( k(x) - k(y) ) * ( k(x) + k(y) - 2y ),
  # where k(z) = max(min(z,b),a) is the value of z capped to [a,b].
  x1 <- pmax( pmin(x,b), a)
  y1 <- pmax( pmin(y,b), a)
  return( sum( (x1 - y1) * (x1 + y1 - 2*y) ) )
}

# Elementary scoring function
Sthet <- function(x, y, theta) {
  if ( sum(y) == 0 ) {
    s <- 1/2 * theta * sum(x > theta)
  } else {
    s <- 1/2 * sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x) )
  }
}

# Elementary scoring functions (vectorized)
Sthet_vec <- function(x, y, theta) {
  # Takes a vector of theta values and thus
  # avoids a further loop for these values
  indic <- outer(x, theta, function(x,y) as.numeric(y < x))
  theta_mat <- matrix(theta, nrow = length(x), ncol = length(theta), byrow = T)
  val <- pmax(y - theta_mat, 0) - pmax(x - theta_mat, 0) - (y - x) * indic
  return(1/2 * val)
}



# Compute MCB, DSC, and UNC for a collection of forecasts
bin_decomp <- function(model, obs, scf = NULL, theta = NULL) {
  # ADD SOME FURTHER EXPLANATION
  if (length(dim(model)) == 2) {
    ncol <- dim(model)[2]
  } else {
    ncol <- 1
    model <- as.matrix(model)
    obs <- as.matrix(obs)
  }
  if (missing(theta)) nrow <- 1 else nrow <- length(theta)
  MCB <- DSC <- UNC <- matrix(0, nrow = nrow, ncol = ncol)
  if (missing(scf)) scf <- function(x,y) colMeans( Sthet_vec(x, y, theta) )
  # Loop over all columns (i.e. bins or models)
  for (i in 1:ncol) {
    meanfc <- model[ ,i]
    y <- obs[ ,i]
    rec <- isoreg(meanfc, y)$yf
    s <- scf(meanfc, y)
    s_rc <- scf(rec, y[order(meanfc)])
    s_mg <- scf(rep(mean(y), length(meanfc)), y)
    MCB[ ,i] <- s - s_rc
    DSC[ ,i] <- s_mg - s_rc
    UNC[ ,i] <- s_mg
  }
  return(list(MCB = MCB, DSC = DSC, UNC = UNC))
}

# Print table of mean scores
scores2teX <- function(scores_pois, scores_quad, mnames, filePath) {
  # Input values:
  # scores_pois - Matrix of daily scores
  # scores_quad - Matrix of daily scores
  # mnames      - Names of the forecast models
  # filePath    - File path for the .tex file
  # Take score values and print teX code to
  # file to get values in tabular environment
  # Compute mean scores
  scores_pois <- colMeans(scores_pois)
  scores_quad <- colMeans(scores_quad)
  k <- length(scores_pois)
  begin <- "\\begin{tabular}{l cc}"
  head <- paste("Model", "pois", "quad \\\\", sep = " & ")
  # Write teX code to file
  write(begin, filePath)
  write("\\hline \\hline", filePath, append = T)
  write(head, filePath, append = T)
  write("\\hline", filePath, append = T)
  for (i in 1:k) {
    # Write score values
    row_pois <- sprintf('%.2f', scores_pois[i])
    row_quad <- sprintf('%.4f', scores_quad[i])
    row <- paste(mnames[i], row_pois, row_quad, sep = " & ")
    row <- paste0(row, " \\\\")
    write(row, filePath, append = T)
  }
  write("\\hline", filePath, append = T)
  write("\\end{tabular}", filePath, append = T)
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

# Compute neighborhood matrix for spatial aggregation
# (See Section D in the Supplement)
neigh_mat <- function(cells, k) {
  # Input values:
  # cells - Data frame of grid cells
  # k     - Integer specifying the size of
  #         the neighborhodd for aggregation
  ncells <- dim(cells)[1]
  # Aggregation will usually be done for small values
  # of k so "sparse = T" makes sense in most cases
  mat <- Matrix(0, nrow = ncells, ncol = ncells, sparse = T)
  for (i in 1:ncells) {
    x <- cells$X[i]
    y <- cells$Y[i]
    xind <- (cells$X <= x + k) & (cells$X >= x - k)
    yind <- (cells$Y <= y + k) & (cells$Y >= y - k)
    mat[i, ] <- as.numeric(xind & yind)
  }
  return(mat)
}
