## Some auxiliary analysis not relevant
## for the main results


## Reproduce plots by Rob Taggart for the trapezodial
## and interval scoring functions for the mean

y <- seq(0,6, by = 0.001)

a <- 3
b <- 6


Sbin <- function(x,y) ( pmax( pmin(x,b), a) - pmax( pmin(y,b), a) ) * (pmax( pmin(x,b), a) + pmax( pmin(y,b), a) - 2 * y)

plot(y, Sbin(4.5,y), ty="l", ylim = c(0,20))


## test for speed
a <- -1
b <- 1/2
n <- 8000
yval <- c(rbinom(30, 1, 0.5), rep(0, n-30) )

# variante 1
Sbin <- function(x,y) ( pmax( pmin(x,b), a) - pmax( pmin(y,b), a) ) * (pmax( pmin(x,b), a) + pmax( pmin(y,b), a) - 2 * y)

tt <- Sys.time()
for (i in 1:1000) {
  ran <- rnorm(n)
  gg <- Sbin(ran,yval)
}
tt <- Sys.time() - tt

# variante 2 (etwas schneller)
Sbin2 <- function(x,y) {
  x1 <- pmax( pmin(x,b), a)
  y1 <- pmax( pmin(y,b), a)
  val <- (x1 - y1) * (x1 + y1 - 2*y)
}

tt <- Sys.time()
for (i in 1:1000) {
  ran <- rnorm(n)
  gg <- Sbin2(ran,yval)
}
tt <- Sys.time() - tt


##################################################################################################################
##################################################################################################################



## Find a suitable grid for the theta values of the 
## Murphy diagram
## Look at data to find suitable binning
q1 <- quantile(as.numeric(model1), probs = seq(0,1, by = 0.05))
q2 <- quantile(as.numeric(model2), probs = seq(0,1, by = 0.05))
q3 <- quantile(as.numeric(model3), probs = seq(0,1, by = 0.05))
q4 <- quantile(as.numeric(model4), probs = seq(0,1, by = 0.05))

plot(seq(0,1, by = 0.05), log(q1), ylim = c(-25,0.5),
     ylab = "log-quantile", xlab = "probability", pch = 19)
points(seq(0,1, by = 0.05), log(q2), col = "darkgreen", pch = 19)
points(seq(0,1, by = 0.05), log(q3), col = "blue", pch = 19)
points(seq(0,1, by = 0.05), log(q4), col = "red", pch = 19)

## Do the same for fully aggregated forecasts
q1 <- quantile(model1_pav, probs = seq(0,1, by = 0.05))
q2 <- quantile(model2_pav, probs = seq(0,1, by = 0.05))
q3 <- quantile(model3_pav, probs = seq(0,1, by = 0.05))
q4 <- quantile(model4_pav, probs = seq(0,1, by = 0.05))

plot(log(q1), seq(0,1, by = 0.05), xlim = c(-2,0.5),
     xlab = "log-quantile", ylab = "probability", pch = 19)
points(log(q2), seq(0,1, by = 0.05), col = "darkgreen", pch = 19)
points(log(q3), seq(0,1, by = 0.05), col = "blue", pch = 19)
points(log(q4), seq(0,1, by = 0.05), col = "red", pch = 19)



##################################################################################################################
##################################################################################################################

## Blockwise evaluation of the forecasts
## temporal blocks of length e.g. 100 days

# source functions for DM tests
source('~/Documents/_temp/Simulations/PP_DMtest.R')

path <- '~/EQData'
load(paste(path, '/scoresFull.RData', sep = ''))
load(paste(path, '/scores_quadFull.RData', sep = ''))
m <- dim(scoresFull)[1]

# blen <- 182     # gives ca. 30 blocks
blen <- 100
nblock <- (m - (m %% blen)) / blen
nmod <- 4

# block vals for logarithmic
scoresFull <- scoresFull[1:(nblock * blen), 1:nmod]
scoresBlock <- matrix(0, nrow = nblock, ncol = nmod)
for (i in 1:nblock) {
  scoresBlock[i, ] <- colMeans( scoresFull[((i-1) * blen + 1):(i * blen), ] )
}

# block vals for quadratic
scores_quadFull <- scores_quadFull[1:(nblock * blen), 1:nmod]
scoresBlock_quad <- matrix(0, nrow = nblock, ncol = nmod)
for (i in 1:nblock) {
  scoresBlock_quad[i, ] <- colMeans( scores_quadFull[((i-1) * blen + 1):(i * blen), ] )
}

hist(scoresBlock[ ,1], breaks = 10)
hist(scoresBlock[ ,2], breaks = 10)
hist(scoresBlock[ ,3], breaks = 10)
hist(scoresBlock[ ,4], breaks = 10)

cols <- c("black", "darkgreen", "blue", "red")

## logarithmic (blocks)
filePath <- "~/Documents/_temp/Case/Plots/plot_log_blocks.pdf"
pdf(filePath, width = 8, height=6)
plot(1:nblock, scoresBlock[ ,1], ty = "l", ylim = c(0.03,3), xlab = "block",
     ylab = "average score", main = "Blockwise logarithmic score")
for (i in 2:4) {
  lines(1:nblock, scoresBlock[ ,i], col = cols[i])
}
legend(1, 3, mnames, col = cols, lwd = 2)
dev.off()

## quadratic (blocks)
filePath <- "~/Documents/_temp/Case/Plots/plot_quad_blocks.pdf"
pdf(filePath, width = 8, height=6)
plot(1:nblock, scoresBlock_quad[ ,1], ty = "l", ylim = c(0,1.9), xlab = "block",
     ylab = "average score", main = "Blockwise quadratic score")
for (i in 2:4) {
  lines(1:nblock, scoresBlock_quad[ ,i], col = cols[i])
}
legend(1, 1.9, mnames, col = cols, lwd = 2)
dev.off()

## compare model 2 to the rest
plot(1:nblock, scoresBlock[ ,3] - scoresBlock[ ,2], main = "M2 vs M3")
abline(h=mean(scoresBlock[ ,3] - scoresBlock[ ,2]))
plot(1:nblock, scoresBlock[ ,1] - scoresBlock[ ,2], main = "M2 vs M1")
abline(h=mean(scoresBlock[ ,1] - scoresBlock[ ,2]))
plot(1:nblock, scoresBlock[ ,4] - scoresBlock[ ,2], main = "M2 vs M4")
abline(h=mean(scoresBlock[ ,4] - scoresBlock[ ,2]))

# Diebold-Mariano
DMpval(scoresBlock)
DMpval(scoresBlock_quad)

# Wilcoxon signed-rank
WSpval(scoresBlock)
WSpval(scoresBlock_quad)


##################################################################################################################
##################################################################################################################


## Compute information gain per earthquake as proposed
## in Rhoades et al. 2011
evt_bins <- data.frame(TI = M4events$TI, N = M4events$N)
#evt_bins <- unique(evt_bins)
numb_eq <- dim(evt_bins)[1]

info_gains <- -1 * c(sum(model1), sum(model2), sum(model3), sum(model4))
for (i in 1:numb_eq) {
  info_gains[1] <- info_gains[1] + log( model1[evt_bins$TI[i], evt_bins$N[i]] )
  info_gains[2] <- info_gains[2] + log( model2[evt_bins$TI[i], evt_bins$N[i]] )
  info_gains[3] <- info_gains[3] + log( model3[evt_bins$TI[i], evt_bins$N[i]] )
  info_gains[4] <- info_gains[4] + log( model4[evt_bins$TI[i], evt_bins$N[i]] )
}
(info_gains <- info_gains / numb_eq)


##################################################################################################################
##################################################################################################################


## Old function to do k-aggregation for events
## works via data.frame instead of sparse matrix
aggregate_evts <- function(evts, nmat) {
  agg <- data.frame(N = evts$N, TI = evts$TI)
  n <- dim(evts)[1]
  for (i in 1:n) {
    TI <- evts$TI[i]
    N <- evts$N[i]
    neigh <- which(as.logical(nmat[ ,N]))
    neigh <- neigh[!(neigh == N)]
    # delete center bin since each bin is its own neighbour!
    add_evts <- data.frame(N = neigh, TI = rep(TI, length(neigh)))
    agg <- rbind(agg, add_evts)
  }
  return(agg)
}

M4events_agg <- aggregate_evts(M4events, nmat)



##################################################################################################################
##################################################################################################################


## Simple map of climatological model
library(maps)
cfile <- paste0(path, '/rate_clima.txt')
clima <- read.table(cfile, header = F, col.names = c("LON", "LAT", "RATE"))

ncols <- 200
pal <- rev(heat.colors(ncols))
lims <- c(min(clima$RATE), max(clima$RATE))
scl <- (clima$RATE - lims[1]) / (lims[2] - lims[1])
vals <- pal[round(scl * (ncols-1)) + 1]

xlim <- c(min(clima$LON), max(clima$LON))
ylim <- c(min(clima$LAT), max(clima$LAT))
border_col <- rgb(0, 0, 0, alpha = 0.4) 
filePath <- "~/Documents/_temp/Case/Plots/map_climatology.pdf"
pdf(filePath, width = 6, height=7.5)
par(mar = c(2/3,2/3,2,2/3))
plot(1, 1, xlim = xlim, ylim = ylim, col = "white", asp = 1.3, xaxt = "n",
     yaxt = "n", xlab = "", ylab = "", main = "Climatological forecast")
points(clima$LON, clima$LAT, pch = 15, col = vals, cex = 0.5)
map('world', fill = F, add = T, col = border_col)
dev.off()