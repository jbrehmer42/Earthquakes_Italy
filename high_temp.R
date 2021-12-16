## Spatial plots of scores, calibration, and
## discrimination over testing area


# Do Data preparation
source('~/Documents/_temp/Case/first_steps.R')
# Load evaluation and plotting functions
source('~/Documents/_temp/Case/eval_functions.R')
source('~/Documents/_temp/Case/plot_functions.R')

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
filePath <- "~/Documents/_temp/Case/Plots/map_skill_score_low.pdf"
source('~/Documents/_temp/Case/map_skill.R')


