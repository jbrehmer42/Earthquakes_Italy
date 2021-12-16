## Auxilary functions for the analysis of
## Italian earthquake data


Spois <- function(x, y, checkNaN = F) {
  s <- -y * log(x) + x
  if (checkNaN) s[is.nan(s)] <- 0
  return(s)
}

Spois2 <- function(x, y) {
  # This is the same as Spois but omitting
  # the if-condition is a speed advantage
  s <- -y * log(x) + x
  s[is.nan(s)] <- 0
  return(mean(s))
}

Squad <- function(x,y) s <- (y - x)^2

Sbin <- function(x, y, a, b) {
  # This function is based on a simplified formula. Let [a,b] the bin,
  # then the restricted quadratic score for this bin (see Taggart, 2020)
  # is given by S(x,y) = ( k(x) - k(y) ) * ( k(x) + k(y) - 2y ),
  # where k(z) = max(min(z,b),a) is the value of z capped to [a,b].
  x1 <- pmax( pmin(x,b), a)
  y1 <- pmax( pmin(y,b), a)
  return( sum( (x1 - y1) * (x1 + y1 - 2*y) ) )
}

Sthet <- function(x, y, theta) {
  if ( sum(y) == 0 ) {
    s <- 1/2 * theta * sum(x > theta)
  } else {
    s <- 1/2 * sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) * (theta < x) )
  }
}

Sthet_vec <- function(x, y, theta) {
  # This is a vectorized version of Sthet
  # It takes a vector of theta values and thus
  # avoids a further loop for these values
  indic <- outer(x, theta, function(x,y) as.numeric(y < x))
  theta_mat <- matrix(theta, nrow = length(x), ncol = length(theta), byrow = T)
  val <- pmax(y - theta_mat, 0) - pmax(x - theta_mat, 0) - (y - x) * indic
  return(1/2 * val)
}


bin_decomp <- function(model, obs, scf = NULL, theta = NULL) {
  # Compute MCB, DSC, and UNC for a collection
  # of forecasts
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


neigh_mat <- function(bins, k) {
  ##mat <- matrix(0, nrow = nbins, ncol = nbins)
  mat <- Matrix(0, nrow = nbins, ncol = nbins, sparse = T)
  for (i in 1:nbins) {
    x <- bins$X[i]
    y <- bins$Y[i]
    xind <- (bins$X <= x + k) & (bins$X >= x - k)
    yind <- (bins$Y <= y + k) & (bins$Y >= y - k)
    mat[i, ] <- as.numeric(xind & yind)
  }
  return(mat)
}
