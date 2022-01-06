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
  return(1/2 * val )
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

# Compute the intersection of the testing region
# of the forecast models and climatological model
region_intersect <- function(clima, cells) {
  ## Compute intersection of testing regions
  ccells <- cbind(cells$LON, cells$LAT)
  cclima <- cbind(clima$LON, clima$LAT)
  z <- unique(rbind(ccells, cclima))
  unicells <- c()
  for (i in 1:dim(z)[1]) {
    icells <- any( (ccells[ ,1] == z[i,1]) & (ccells[ ,2] == z[i,2]) )
    iclima <- any( (cclima[ ,1] == z[i,1]) & (cclima[ ,2] == z[i,2]) )
    if (icells && iclima) unicells <- rbind(unicells, z[i, ])
  }
  # compute subsets as logical indices
  model.subs <- apply(ccells, 1, function(x) any( (x[1] == unicells[ ,1]) & (x[2] == unicells[ ,2]) ) )
  clima.subs <- apply(cclima, 1, function(x) any( (x[1] == unicells[ ,1]) & (x[2] == unicells[ ,2]) ) )
  return(list(model = model.subs, clima = clima.subs))
}
