# Earthquakes_Italy

This is the repository for an extended case study of earthquake forecasting in
Italy. Parts of this case study and the corresponding code have been used in

Jonas Brehmer, Tilmann Gneiting, Marcus Herrmann, Warner Marzocchi, Martin
Schlather, and Kirstin Strokorb (2022+). Comparative evaluation of point
process forecasts. 


## Code

The following files reproduce the simulation study:
- *eval1.R*
- ...

Other files:
- *Resin_meanreldiag.R*  Code for mean reliability diagram by Johannes Resin


## Data

- models (four files): Two-dimensional arrays which contain the forecasts of
the models. The rows represent different model run times. The columns represent
the grid cells into which the testing region is split.
- time stamps file: Time stamps belonging to the rows of the model arrays, i.e.
these are the times of the model runs
- cells file: Longitude and latitude of the grid cells belonging to the columns
of the model arrays. They specify the centers of a boxes with edge length 0.1
degrees and are numbered consecutively.
- catalog file: Contains the details (time stamp, magnitude, etc.) of
earthquakes in the testing region and during the testing period (but time and
space of model forecasts and catalog do not exactly match)
- climatology file: Rates for a selection of longitude and latitude values. If 
appropriately scaled, they can be understood as climatological forecast which
are constant in time. The scaling depends on the assumed number of events in a
7-day period, see Mail by Warner Marzocchi (08.09.21)

The code loads all the data to RAM, so on a machine with 8GB or less memory
this can be problematic and crash R. This is also the reason for many 'rm' and
'gc' calls throughout the code.
