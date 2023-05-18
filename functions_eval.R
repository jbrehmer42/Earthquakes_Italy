###################################
## Auxiliary functions - General ##

# Poisson scoring function
# Caution: Have to ensure that reports x are never zero
S_pois <- function(x, y)  -y * log(x) + x

# Quadratic scoring function
S_quad <- function(x,y) (y - x)^2

# Mean quadratic scoring function
S_quad2 <- function(x,y) mean( (y - x)^2 )

# First component of the normalized spatial score
S_spat <- function(x, y) -y * log(x)


S_pois2 <- function(x, y) {
  # Mean Poisson scoring function
  # This is the same as S_pois but omitting
  # the if-condition is a speed advantage
  s <- -y * log(x) + x
  s[is.nan(s)] <- 0
  return(mean(s))
}


S_bin <- function(x, y, a, b) {
  # Quadratic interval restricted scoring function.
  # This function is based on a simplified formula. For any [a,b]
  # the restricted quadratic score for this interval (see Taggart, 2020)
  # is given by S(x,y) = ( k(x) - k(y) ) * ( k(x) + k(y) - 2y ),
  # where k(z) = max(min(z,b),a) is the value of z capped to [a,b].
  x1 <- pmax( pmin(x,b), a)
  y1 <- pmax( pmin(y,b), a)
  return( sum( (x1 - y1) * (x1 + y1 - 2*y) ) )
}


S_theta <- function(x, y, theta) {
  # Elementary scoring function for the mean following Ehm et al. (2016)
  if ( sum(y) == 0 ) {
    s <- theta * sum(x > theta)
  } else {
    s <- sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) *
                      (theta < x) )
  }
  return(s)
}


S_theta_vec <- function(x, y, theta) {
  # Elementary scoring functions for the mean (vectorized.)
  # Takes a vector of theta values and thus avoids a
  # further loop for these values in cell_decomposition
  ind <- outer(x, theta, function(x,y) as.numeric(y < x))
  theta_mat <- matrix(theta, nrow = length(x), ncol = length(theta), byrow = T)
  val <- pmax(y - theta_mat, 0) - pmax(x - theta_mat, 0) - (y - x) * ind
  return(val)
}




cell_decomposition <- function(model, obs, scf = NULL, theta = NULL) {
  # Compute the score decomposition into MCB, DSC, and UNC.
  #
  # Input values:
  # model  -  matrix or vector of forecast values
  # obs    -  observation matrix
  # scf    -  scoring function (optional)
  # theta  -  parameter values for elementary scoring
  #           function for the mean (only if scf missing)
  # Provide either scf or theta! If scf is given, then the
  # decomposition based on scf is calculated for all cells.
  # If theta is given, then the decomposition based on the
  # elementary score with parameter theta is calculated for
  # all cells and all theta. If model is a vector, then we
  # have only one cell.
  # Check input to determine output dimensions
  if (length(dim(model)) == 2) {
    ncol <- dim(model)[2]
  } else {
    ncol <- 1
    model <- as.matrix(model)
    obs <- as.matrix(obs)
  }
  if (missing(theta)) nrow <- 1 else nrow <- length(theta)
  MCB <- DSC <- UNC <- matrix(0, nrow = nrow, ncol = ncol)
  if (missing(scf)) scf <- function(x,y) colMeans( S_theta_vec(x, y, theta) )
  # Loop over all cells
  for (i in 1:ncol) {
    if (sum(obs[ ,i]) == 0) {
      # if obs[ ,i] = 0 everywhere, then there are no events
      # in this bin. This happens quite often, so catching this
      # case is a speed advantage.
      MCB[ ,i] <- scf(model[ ,i], obs[ ,i])
    } else {
      x <- model[ ,i]
      y <- obs[ ,i]
      rec <- isoreg(x, y)$yf
      s <- scf(x, y)
      s_rc <- scf(rec, y[order(x)])
      s_mg <- scf(rep(mean(y), length(x)), y)
      MCB[ ,i] <- s - s_rc
      DSC[ ,i] <- s_mg - s_rc
      UNC[ ,i] <- s_mg
    }
  }
  return(list(MCB = MCB, DSC = DSC, UNC = UNC))
}


score_tex_table <- function(scores_pois, scores_quad, model_names, file_path) {
  # Print lateX table of mean scores
  #
  # Input values:
  # scores_pois  -  Matrix of daily scores
  # scores_quad  -  Matrix of daily scores
  # model_names  -  Names of the forecast models
  # file_path    -  File path for the .tex file
  # Take score values and print teX code to
  # file to get values in tabular environment
  # Compute mean scores
  scores_pois <- colMeans(scores_pois)
  scores_quad <- colMeans(scores_quad)
  k <- length(scores_pois)
  begin <- "\\begin{tabular}{l cc}"
  head <- paste("Model", "pois", "quad \\\\", sep = " & ")
  # Write teX code to file
  write(begin, file_path)
  write("\\hline \\hline", file_path, append = T)
  write(head, file_path, append = T)
  write("\\hline", file_path, append = T)
  for (i in 1:k) {
    # Write score values
    row_pois <- sprintf('%.2f', scores_pois[i])
    row_quad <- sprintf('%.4f', scores_quad[i])
    row <- paste(model_names[i], row_pois, row_quad, sep = " & ")
    row <- paste0(row, " \\\\")
    write(row, file_path, append = T)
  }
  write("\\hline", file_path, append = T)
  write("\\end{tabular}", file_path, append = T)
}



neigh_mat <- function(cells, k) {
  # Compute neighborhood matrix for spatial aggregation
  #
  # Input values:
  # cells - Data frame of grid cells
  # k     - Integer specifying the size of
  #         the neighborhood for aggregation
  # A neighborhood matrix is a binary square matrix where
  # a 1 at position (i,j) indicates that cells i and j 
  # are in each others neighborhoods.
  n_cells <- dim(cells)[1]
  # Aggregation will usually be done for small values
  # of k so "sparse = T" makes sense in most cases
  mat <- Matrix(0, nrow = n_cells, ncol = n_cells, sparse = T)
  for (i in 1:n_cells) {
    x <- cells$X[i]
    y <- cells$Y[i]
    x_index <- (cells$X <= x + k) & (cells$X >= x - k)
    y_index <- (cells$Y <= y + k) & (cells$Y >= y - k)
    mat[i, ] <- as.numeric(x_index & y_index)
  }
  return(mat)
}


region_intersect <- function(clima, cells) {
  # Compute the intersection of the testing region
  # of the forecast models and climatological model
  #
  # Input values:
  # clima  -  data frame with longitudes and latitudes
  #           of the climatological model region
  # clima  -  data frame with longitudes and latitudes
  #           of the four forecast models' region
  coord_cells <- cbind(cells$LON, cells$LAT)
  coord_clima <- cbind(clima$LON, clima$LAT)
  coord_unique <- unique(rbind(coord_cells, coord_clima))
  coord_intersect <- c()
  for (i in 1:dim(coord_unique)[1]) {
    icells <- any( (coord_cells[ ,1] == coord_unique[i,1]) &
                     (coord_cells[ ,2] == coord_unique[i,2]) )
    iclima <- any( (coord_clima[ ,1] == coord_unique[i,1]) &
                     (coord_clima[ ,2] == coord_unique[i,2]) )
    if (icells && iclima) coord_intersect <- rbind(coord_intersect,
                                                   coord_unique[i, ])
  }
  # compute subsets as logical indices
  compare_intersect <- function(x) any( (x[1] == coord_intersect[ ,1]) &
                                          (x[2] == coord_intersect[ ,2]) )
  model_subset <- apply(coord_cells, 1, compare_intersect)
  clima_subset <- apply(coord_clima, 1, compare_intersect)
  return(list(model = model_subset, clima = clima_subset))
}
