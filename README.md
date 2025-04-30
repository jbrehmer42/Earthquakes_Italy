# Earthquakes_Italy

This repository provides the source code (_R_ language) for an extended case study of earthquake forecasting in
Italy:

> Brehmer, J. R., Kraus, K., Gneiting, T., Herrmann, M., and Marzocchi, W. (2025).
> Enhancing the Statistical Evaluation of Earthquake Forecasts—An Application to Italy.
> _Seismological Research Letters, 96_(3).
> doi: [10.1785/0220240209](https://doi.org/10.1785/0220240209).
> arXiv: [2405.10712](https://arxiv.org/abs/2405.10712)

Parts of this case study and the corresponding code have been also used in:

> Brehmer, J. R., Gneiting, T., Herrmann, M., Marzocchi, W., Schlather, M., and Strokorb, K. (2024).
> Comparative evaluation of point process forecasts.
> _Annals of the Institute of Statistical Mathematics 76_(1), 47–71.
> doi: [10.1007/s10463-023-00875-5](https://doi.org/10.1007/s10463-023-00875-5).
> arXiv: [2103.11884](https://arxiv.org/abs/2103.11884)

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

The code was developed on a machine with 16GB of RAM. This is required if for the analysis all the models are
kept in memory simultaneously.

## Data input

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


## Details on slight differences to previously reported scores

Compared to two previous studies, we report slightly different scores due to correcting the
spatial binning of earthquakes that occurred exactly on grid cell boundaries (bin edges).
These binning discrepancies originate from numerical inaccuracies in floating-point arithmetics.
Of the 262 target earthquakes in the catalog, this affects seven events.


**1. Brehmer _et al._ 2024**

Poisson scores in Table 1–3 slightly differ compared to Table 1 in Brehmer _et al._ 2024 (see second reference above).

Three target earthquakes previously had been unnecessarily excluded from analysis:
```
DateTime,                Lat,  Lon,   Depth, Mag
2016-05-30 20:24:20.460, 42.7, 11.976,  7.9, 4.1
2016-08-26 04:28:25.890, 42.6, 13.29,  10.9, 4.8
2019-10-25 04:31:38.200, 39.7, 15.432, 12.1, 4.4
```
(their _Latitude_ component falls exactly on a cell boundary)

Additionally, we now adhere to [_pyCSEP_-style binning](https://docs.cseptesting.org/reference/generated/csep.utils.calc.bin1d_vec.html)
which includes lower boundaries to a cell and excludes upper boundaries,
whereas Brehmer _et al._ 2024 used the opposite rule.
This resulted in four earthquakes being incorrectly assigned to a neighboring cell:
```
DateTime,                Lat,       Lon, Depth, Mag
2006-02-27 04:34:01.830, 38.155,   15.2,   9.2, 4.1
2009-04-06 01:42:49.970, 42.3,   13.429,  10.5, 4.2
2016-12-09 07:21:50.170, 44.33,    10.5,   7.6, 4.0
2016-12-11 12:54:52.860, 42.9,   13.113,   8.3, 4.3
```
(either the _Latitude_ or _Longitude_ component falls exactly on a cell boundary)


**2. Herrmann & Marzocchi 2023**

IGPE values in Table 3 slightly differ compared to Table 1 in:

> Herrmann, M., and W. Marzocchi (2023).
> Maximizing the forecasting skill of an ensemble model.
> _Geophysical Journal International 234_(1), 73–87.
> doi: [10.1093/gji/ggad020](https://doi.org/10.1093/gji/ggad020)

(when transforming reported values to use _ETAS_LM_ as reference instead of _SMA_).

Herrmann & Marzocchi 2023 used [_pyCSEP_'s binning function](https://docs.cseptesting.org/reference/generated/csep.utils.calc.bin1d_vec.html),
but it previously [didn't properly account for floating point precision](https://github.com/SCECcode/pycsep/issues/255),
resulting in five earthquakes being incorrectly assigned to a neighboring cell:
```
DateTime,                Lat,     Lon, Depth, Mag
2009-04-06 01:42:49.970, 42.3, 13.429,  10.5, 4.2
2016-05-30 20:24:20.460, 42.7, 11.976,   7.9, 4.1
2016-08-26 04:28:25.890, 42.6,  13.29,  10.9, 4.8
2016-12-11 12:54:52.860, 42.9, 13.113,   8.3, 4.3
2019-10-25 04:31:38.200, 39.7, 15.432,  12.1, 4.4
```
(only their _Latitude_ component falls exactly on a cell boundary; those five events were already mentioned above)
