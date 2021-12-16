## Create four maps for the skill score
## use truncation for negative values

# neg_col <- 0.66       # blue # hue in hsv colors
# pos_col <- 0          # red  # hue in hsv colors
pos_col <- 0.35       # green
neg_col <- 0          # red

pdf(filePath, width = 8, height = 7)
layout(matrix(c(1,2,5,3,4,5), 2, 3, byrow = TRUE), widths = c(0.42, 0.42, 0.16))
## start with maps
par(mar = c(2/3,2/3,1.5,1/3))
for (i in 1:4) {
  skills <- pmax(SKL_map[ ,i], skl_min)
  scl_pos <-  skills[skills >= 0]
  scl_neg <-  - (skills[skills < 0] - skl_min) / skl_min
  vals <- rep("", nbins)
  vals[skills >= 0] <- hsv(pos_col, scl_pos, v = 1 - 0.4*scl_pos )
  vals[skills < 0]  <- hsv(neg_col, 1 - scl_neg)
  # main <- paste0("Model ", i)
  main <- mnames[i]
  # plot_map(vals, main = main)
  plot_map(vals, main = main, evts = M4events)
}
## add color bar
par(mar = c(4,1,1,4), mgp = c(3,0,-1))
len_pos <- (ncols/10 - 1)/3
len_neg <- 2 * len_pos
nticks <- 3
plot(1,1, col = "white", xlim = c(0,1), ylim = c(-len_neg, len_pos), asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
labs <- round(seq(0, skl_max, len = nticks+1), 2)
labs <- c(round(seq(skl_min, 0, len = 2*nticks), 2), labs[2:length(labs)])
ats <- seq(0, len_pos, len = nticks+1)
ats <- c(seq(-len_neg, 0, len = 2*nticks), ats[2:length(ats)])
labs[1] <- paste0("< ", skl_min)
axis(4, at = ats, labels = labs, las = 1)
kk <- 50
for (l in 1:(kk-1)) {
  rect(0, (l-1) * len_pos/(kk - 1), 1,  l * len_pos/(kk - 1), 
       col = hsv(pos_col, l/(kk-1), v = 1 - 0.4 * l/(kk-1)), border = NA)
}
kk <- 2*kk
for (l in 1:(kk-1)) {
  rect(0, - (l-1) * len_neg/(kk - 1), 1,  - l * len_neg/(kk - 1), 
       col = hsv(neg_col, l/(kk-1)), border = NA)
}
dev.off()
