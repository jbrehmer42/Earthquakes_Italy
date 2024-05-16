# Earthquakes_Italy

This is the repository for an extended case study of earthquake forecasting in
Italy:

Jonas Brehmer, Kristof Kraus, Tilmann Gneiting, Marcus Herrmann, and
Warner Marzocchi (2024+). Comparative evaluation of earthquake forecasting 
models: An application to Italy.

Parts of this case study and the corresponding code have been used in

Jonas Brehmer, Tilmann Gneiting, Marcus Herrmann, Warner Marzocchi, Martin
Schlather, and Kirstin Strokorb (2023). Comparative evaluation of point
process forecasts.

## Differences to previously reported predictive performance

We report slightly different Poisson scores compared to

J. Brehmer, T. Gneiting, M. Herrmann, W. Marzocchi, M.
Schlather, and K. Strokorb (2023). Comparative evaluation of point
process forecasts.

Due to numerical inaccuracies three earthquakes had been unnecessarily excluded 
from analysis (these were the earthquake going with the timestamp 2016-05-30 20:24:20.460,
2016-08-26 04:28:25.890, and 29019-10-25 04:31:38.200). 
Additionally, we adhere to the CSEP binning which includes lower boundaries to a cell and
excludes upper boundaries, whereas Brehmer et al. used the opposite rule. This results in four 
earthquakes being assigned to a neighboring cells.

Likewise, we report different IGPE values compared to

M. Herrmann and W. Marzocchi (2023). Maximizing the forecasting skill of an ensemble
model. Geophysical Journal International.

The reason is again a slightly different binning of earthquakes registered on cell boundaries.
While Herrmann et al. also use the CSEP binning, numerical inaccuracies entailed that five
earthquakes had been assigned to a wrong cell in terms of the CSEP binning.

## Code

The tables and plots can be reproduced with the following files:

- **data_prep.R** Loads all the data (see below) and pre-processes them. Called at
the start of all other main scripts.
- **functions_prep.R** Provides functions to load and preprocess the data
- **precompute.R** pre-compute values for empirical CDFs, Murphy 
diagrams, score component plot and reliability diagram
- **plots_and_tables.R** comprises all the functionality to calculate the values for the
tables and to create the plots of the publication apart from simulation study plots
- **sim_study_calibr.R** manipulate LM forecast in different ways and visualize the effects
on its reliability curve
- **sim_study_tests.R** compare CSEP t-Test and Diebold Mariano test on forecasts derived from the LM, FCM and LG model
- **utils.R** define plot theme, and scoring and test functions


## Data

- **models** (four files): Two-dimensional arrays which contain the forecasts of
the models. The rows represent different model run times. The columns represent
the grid cells into which the testing region is split.
- **time stamps** file: Time stamps belonging to the rows of the model arrays, i.e.
these are the times of the model runs
- **cells** file: Longitude and latitude of the grid cells belonging to the columns
of the model arrays. They specify the centers of a boxes with edge length 0.1
degrees and are numbered consecutively.
- **catalog** file: Contains the details (time stamp, magnitude, etc.) of
earthquakes in the testing region and during the testing period (but time and
space of model forecasts and catalog do not exactly match)
- **climatology** file: Rates for a selection of longitude and latitude values. If 
appropriately scaled, they can be understood as climatological forecast which
are constant in time. The scaling depends on the assumed number of events in a
7-day period, see Mail by Warner Marzocchi (08.09.21)

The code was developed on a machine with 16GB of RAM. This is required if for the analysis all the models are
kept in memory simultaneously.
