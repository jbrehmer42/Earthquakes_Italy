## Reduce testing region to intersection of
## forecast models and climatological model

## Compute intersection of testing regions
cbins <- cbind(bins$LON, bins$LAT)
cclima <- cbind(clima$LON, clima$LAT)
z <- unique(rbind(cbins, cclima))
unibins <- c()
for (i in 1:dim(z)[1]) {
  ibins <- any( (cbins[ ,1] == z[i,1]) & (cbins[ ,2] == z[i,2]) )
  iclima <- any( (cclima[ ,1] == z[i,1]) & (cclima[ ,2] == z[i,2]) )
  if (ibins && iclima) unibins <- rbind(unibins, z[i, ])
}

# plot(1, 1, xlim = xlim, ylim = ylim, col = "white", asp = 1.3, xaxt = "n",
#      yaxt = "n", xlab = "", ylab = "", main = "Climatological forecast")
# points(unibins[ ,1], unibins[ ,2], pch = 15, col = "green", cex = 0.3)

## convert to logical and
## compute subsets
bin.subs <- apply(cbins, 1, function(x) any( (x[1] == unibins[ ,1]) & (x[2] == unibins[ ,2]) ) )
clima.subs <- apply(cclima, 1, function(x) any( (x[1] == unibins[ ,1]) & (x[2] == unibins[ ,2]) ) )

bins <- bins[bin.subs, ]
clima <- clima[clima.subs, ]
stopifnot(all(clima$LON == bins$LON))
stopifnot(all(clima$LAT == bins$LAT))
nbins <- dim(bins)[1]

## subset forecast models and observations
obs <- obs[ ,bin.subs]
model1 <- model1[ ,bin.subs]
model2 <- model2[ ,bin.subs]
model3 <- model3[ ,bin.subs]
model4 <- model4[ ,bin.subs]

rm(unibins, cbins, cclima, ibins, iclima, clima.subs, z)


## Additional code for number of
## events in (reduced) testing region
# n4 <- dim(M4events)[1]
# subs4 <- rep(0, n4)
# for (i in 1:n4) {
#   subs4[i] <- any( M4events$N[i] == bins$N)
# }
# n_per_year <- sum(subs4) / (dim(model1)[1] / 365)
