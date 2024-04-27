# Earthquakes_Italy

This is the repository for an extended case study of earthquake forecasting in
Italy. Parts of this case study and the corresponding code have been used in

Jonas Brehmer, Tilmann Gneiting, Marcus Herrmann, Warner Marzocchi, Martin
Schlather, and Kirstin Strokorb (2023). Comparative evaluation of point
process forecasts

and

Jonas Brehmer, Kristof Kraus, Tilmann Gneiting, Marcus Herrmann, and
Warner Marzocchi (2024+). Comparative evaluation of earthquake forecasting 
models: An application to Italy.


## Code - Comparative evaluation of point process forecasts

The following files create the diagrams and maps of the paper:

- **data_prep.R** Loads all the data (see below) and pre-processes them. Called at
the start of all other main scripts.
- **temporal_scores_Murphy.R** Calculates overall scores and scores during the 
testing period for the models and their spatially aggregated versions. The
latter allow for Murphy diagrams.
- **mean_calibration.R** Computes the score decomposition into MCB, DSC, and UNC
and plots mean reliability diagrams.
- **spatial_scores.R** Computes scores, skill scores, score differences, and score
decomposition and plots them on a map of the testing region.
- **spatial_scores_aggregated.R** Similar to spatial_scores.R, but now the
quantities (except the skill scores) are calculated and plotted for the
spatially aggregated models.
- **spatial_skill_scores_aggregated.R** Similar to spatial_scores.R, but now the
skill scores are calculated and plotted for the spatially aggregated models.
- **mean_calibration_ex_large.R** Plot mean reliability diagrams similar to
mean_calibration.R, but now exclude very large earthquakes from computation.
Also plot diagrams for each week day separately.
- **test_inter_daily_updates.R** Experimental. Test whether there are changes in
the mean reliability diagrams if we include model updates which happen during
earthquake sequences.
- **functions_plot.R** Provides functions to create the plots.
- **functions_eval.R** Provides functions to evaluate the models.
- **functions_prep.R** Provides functions to load and process the data in
data_prep.R.

Other files:

- **Resin_meanreldiag.R**  Code for mean reliability diagram by Johannes Resin.
Basis for the plot_reliability function in functions_plot.R.

## Code - Comparative evaluation of earthquake forecasting models: An application to Italy

The tables and plots can be reproduced with the following files:

- **data_prep.R** Loads all the data (see below) and pre-processes them. Called at
the start of all other main scripts.
- **functions_prep.R** Provides functions to load and process the data in
- **full-evaluation_heavy-computation.R** pre-compute values for empirical CDFs, Murphy 
diagrams, score component plot and reliability diagram
- **full-evaluation_plots.R** comprises all the functionality to calculate the values for the
tables and to create the plots of the publication regarding the analysis of the earthquake 
forecasting models
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

The code loads all the data to RAM, so on a machine with 8GB or less memory
this can be problematic and crash R. On my current machine with 8GB RAM
everything runs fine if I close other memory-intensive applications. This is
also the reason for many 'rm' and 'gc' calls throughout the code.
