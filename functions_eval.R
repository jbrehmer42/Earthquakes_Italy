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

# Mean Poisson scoring function
S_pois2 <- function(x, y) {
  # This is the same as S_pois but omitting
  # the if-condition is a speed advantage
  s <- -y * log(x) + x
  s[is.nan(s)] <- 0
  return(mean(s))
}

# Quadratic interval restricted scoring function
S_bin <- function(x, y, a, b) {
  # This function is based on a simplified formula. For any [a,b]
  # the restricted quadratic score for this interval (see Taggart, 2020)
  # is given by S(x,y) = ( k(x) - k(y) ) * ( k(x) + k(y) - 2y ),
  # where k(z) = max(min(z,b),a) is the value of z capped to [a,b].
  x1 <- pmax( pmin(x,b), a)
  y1 <- pmax( pmin(y,b), a)
  return( sum( (x1 - y1) * (x1 + y1 - 2*y) ) )
}

# Elementary scoring function
S_theta <- function(x, y, theta) {
  if ( sum(y) == 0 ) {
    s <- 1/2 * theta * sum(x > theta)
  } else {
    s <- 1/2 * sum( pmax(y-theta, 0) - pmax(x-theta, 0) - (y - x) *
                      (theta < x) )
  }
}

# Elementary scoring functions (vectorized)
S_theta_vec <- function(x, y, theta) {
  # Takes a vector of theta values and thus
  # avoids a further loop for these values
  ind <- outer(x, theta, function(x,y) as.numeric(y < x))
  theta_mat <- matrix(theta, nrow = length(x), ncol = length(theta), byrow = T)
  val <- pmax(y - theta_mat, 0) - pmax(x - theta_mat, 0) - (y - x) * ind
  return(1/2 * val )
}



# Compute MCB, DSC, and UNC for a collection of forecasts
bin_decomposition <- function(model, obs, scf = NULL, theta = NULL) {
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
  if (missing(scf)) scf <- function(x,y) colMeans( S_theta_vec(x, y, theta) )
  # Loop over all columns (i.e. bins or models)
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

# Print table of mean scores
score_tex_table <- function(scores_pois, scores_quad, model_names, file_path) {
  # Input values:
  # scores_pois - Matrix of daily scores
  # scores_quad - Matrix of daily scores
  # model_names      - Names of the forecast models
  # file_path    - File path for the .tex file
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


# Compute neighborhood matrix for spatial aggregation
neigh_mat <- function(cells, k) {
  # Input values:
  # cells - Data frame of grid cells
  # k     - Integer specifying the size of
  #         the neighborhodd for aggregation
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

# Compute the intersection of the testing region
# of the forecast models and climatological model
region_intersect <- function(clima, cells) {
  ## Compute intersection of testing regions
  # Param
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
