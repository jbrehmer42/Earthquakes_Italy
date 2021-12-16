## First steps for case study
## Data preparation and graphics

path <- '~/EQData'
# Load package for sparse matrices
library(Matrix)


## 1 ##
## Forecast models
file1 <- paste0(path, '/ETAS_LM.txt.xz')
file2 <- paste0(path, '/ETES_FMC.txt.xz')
file3 <- paste0(path, '/STEP_LG.txt.xz')
file4 <- paste0(path, '/Bayesian_corr_27_10.txt.xz')

model1 <- as.matrix( read.table(xzfile(file1)) )
model2 <- as.matrix( read.table(xzfile(file2)) )
model3 <- as.matrix( read.table(xzfile(file3)) )
model4 <- as.matrix( read.table(xzfile(file4)) )

attr(model1, "dimnames") <- NULL
attr(model2, "dimnames") <- NULL
attr(model3, "dimnames") <- NULL
attr(model4, "dimnames") <- NULL


mnames <- c("LM", "FMC", "LG", "SMA")
cols <- c("black", "darkgreen", "blue", "red")

# Load climatological model (constant in time)
cfile <- paste0(path, '/rate_clima.txt')
clima <- read.table(cfile, header = F, col.names = c("LON", "LAT", "RATE"))
#evts_per7 <- 25.95 * 7/365      # See Mail by Warner (08.09.21)
evts_per7 <- 12 * 7/365      # See Mail by Warner (08.09.21)
#evts_per7 <- 16.97 * 7/365      # See n_per_year variable
clima$RATE <- clima$RATE * evts_per7

# Load time stamps corresponding to matrix rows
tfile <- paste(path, '/meta_rows.txt', sep='')
times <- read.table(tfile, col.names = c("DD", "MM", "YY", "H", "M", "S"))

# Filter out all days with multiple model runs and
# before the end of testing period
get_days <- function(times) {
  # identify days which have multiple forecasts
  n <- dim(times)[1]
  ind <- rep(F, n)
  for (i in 2:n) {
    ind[i] <- all(times[i-1, 1:3] == times[i, 1:3])
  }
  return(ind)
}
max.date <- c(20, 5, 2020)

unik <- !get_days(times)
max.ind <- which((times$DD == max.date[1]) & (times$MM == max.date[2]) & (times$YY == max.date[3]))
if ( (length(max.ind) == 1) & (max.ind < length(unik)) ) unik[(max.ind+1):length(unik)] <- FALSE

times2 <- times[unik, ]
model1 <- model1[unik, ]
model2 <- model2[unik, ]
model3 <- model3[unik, ]
model4 <- model4[unik, ]



## 2 ##
## Events
efile <- paste(path, '/meta_catalogo.txt', sep='')
events <- read.table(efile, col.names = c("YY", "MM", "DD", "H", "M", "S", "LAT", "LON", "DEP", "MAG"))

# Using M >= 4 is equivalent to using M >= 3.95 earthquakes
M4ind <- (events$MAG >= 4)
M4events <- events[M4ind, ]
#plot(events$LON[M4ind], events$LAT[M4ind])

# Load bin data
bfile <- paste(path, '/meta_column.csv', sep='')
bins <- read.csv(bfile, header=F, col.names = c("LON", "LAT", "N"))
bins$N <- bins$N + 1

# Calculate x-y-coordinates for bins (relative to testing region)
xvals <- sort(unique(bins$LON))
yvals <- sort(unique(bins$LAT))
bins$X <- sapply(bins$LON, function(x) which(x == xvals))
bins$Y <- sapply(bins$LAT, function(x) which(x == yvals))

# Assign bin numbers to events
size_LON <- 0.1
size_LAT <- 0.1
get_bins <- function(evts, bins) {
  # assign bin number to every event
  # bin numbers start with 1
  n <- dim(evts)[1]
  ind <- rep(0, n)
  for (i in 1:n) {
    bin_ri <- bins$LON + 0.5 * size_LON
    bin_le <- bins$LON - 0.5 * size_LON
    bin_lo <- bins$LAT - 0.5 * size_LAT
    bin_up <- bins$LAT + 0.5 * size_LAT
    isLON <- (bin_le < evts$LON[i]) & (evts$LON[i] <= bin_ri)
    isLAT <- (bin_lo < evts$LAT[i]) & (evts$LAT[i] <= bin_up)
    if (any(isLON & isLAT) ) {
      # add 1 so bin numbers start with 1
      ind[i] <- bins$N[isLON & isLAT]
    } else {
      # set index to -1 outside testing region
      ind[i] <- -1
    }
  }
  return(ind)
}

M4ind2 <- get_bins(M4events, bins)
M4events <- M4events[(M4ind2 > 0), ]
M4events$N <- M4ind2[M4ind2 > 0]

# Assign day number (time index) to events
get_tindex <- function(evts, times) {
  # assign day number/time index to every event
  # day numbers are given by times data.frame
  n <- dim(evts)[1]
  ind <- rep(-1, n)
  for (i in 1:n) {
    DD <- evts$DD[i]
    MM <- evts$MM[i]
    YY <- evts$YY[i]
    gind <- (times$DD == DD) & (times$MM == MM) & (times$YY == YY)
    if (any(gind)) ind[i] <- which(gind)
  }
  return(ind)
}

M4ind3 <- get_tindex(M4events, times2)
M4events <- M4events[(M4ind3 > 0), ]
M4events$TI <- M4ind3[M4ind3 > 0]

# Transform to weekly events
nbins <- dim(model1)[2]
ndays <- dim(model1)[1]
obs <- Matrix(0, ncol = nbins, nrow = ndays, sparse = T)
for (i in 1:ndays) {
  ind <- (M4events$TI >= i) & (M4events$TI < i + 7)
  if (any(ind)) {
    weekbin <- tabulate(M4events$N[ind], nbins = nbins)
    obs[i, ] <- weekbin
  } else next
}

## Clean up
rm(M4ind, M4ind2, M4ind3, unik, xvals, yvals, ind, weekbin, i,
   file1, file2, file3, file4, efile, tfile, bfile)
