## Spatial plots of scores, calibration, and
## discrimination over testing area


# Do Data preparation
source('~/Documents/_temp/Case/data_prep.R')
# Load evaluation and plotting functions
source('~/Documents/_temp/Case/functions_eval.R')
source('~/Documents/_temp/Case/functions_plot.R')

## Compute values for maps
MCB_map <- DSC_map <- matrix(0, nrow = nbins, ncol = 4)
for (i in 1:4) {
  mname <- paste0("model", i)
  decomp <- bin_decomp(get(mname), obs, scf = Spois2)
  MCB_map[ ,i] <- decomp$MCB
  DSC_map[ ,i] <- decomp$DSC
}
UNC_map <- decomp$UNC
SCR_map <- MCB_map - DSC_map + matrix(UNC_map, ncol = 4, nrow = nbins)

## Do spatial plots

# create color vector 
#pal <- paste("gray", round(scl * 100), sep = "")   # white = high score
ncols <- 200
pal <- rev(heat.colors(ncols))
#pal <- gray.colors(ncols)

## Create maps for the scores
lims <- c(min( log(SCR_map) ), max( log(SCR_map) )) + 0.05 * c(-1,1)
file <- "~/Documents/_temp/Case/Plots/map_score_log"
plot_figure(SCR_map, pal, lims, ncols, file, evts = TRUE)

## Create maps for miscalibration
lims[1] <- min( log(MCB_map) ) - 0.05
lims[2] <- max( log(MCB_map) ) + 0.05
file <- "~/Documents/_temp/Case/Plots/map_mcb_log"
plot_figure(MCB_map, pal, lims, ncols, file, evts = TRUE)

## Create maps for discrimination
offs <- 1e-5
lims[1] <- log(offs)
lims[2] <- max( log(DSC_map + offs)) + 0.1
file <- "~/Documents/_temp/Case/Plots/map_dsc_log"
plot_figure(DSC_map, pal, lims, ncols, file, offset = offs)

## Create maps for score differences
rootName <- "~/Documents/_temp/Case/Plots/map_score_diff"
for (i in 1:4) {
  for (j in 1:4) {
    if (i >= j) next
    file <- paste0(rootName, "_", i, j)
    scr_abs <- 1.01 * max( abs(SCR_map[ ,i] - SCR_map[ ,j]) )
    plot_diff_maps(SCR_map, c(i,j), scr_abs, file)
  }
}

## THESE REMAIN TO DO....
## Create maps for skill scores
## Since we use the climatology the testing
## region changes we have to modify SCR_map

# Switch to reduced testing region
source('~/Documents/_temp/Case/clima_subset.R')
# Load plotting functions
source('~/Documents/_temp/Case/plot_functions.R')

# Define climatology
model5 <- matrix(clima$RATE, ncol = nbins, nrow = ndays, byrow = TRUE)
# Compute decomposition
decomp <- bin_decomp(model5, obs, scf = Spois2)
Sclima <- decomp$MCB - decomp$DSC + decomp$UNC
SCR_map <- SCR_map[bin.subs, ]
SKL_map <- ( matrix(Sclima, ncol = 4, nrow = nbins) - SCR_map ) / matrix(Sclima, ncol = 4, nrow = nbins)
# Do plotting
skl_max <- 1
skl_min <- -10
filePath <- "~/Documents/_temp/Case/Plots/map_skill_score.pdf"
source('~/Documents/_temp/Case/map_skill.R')


