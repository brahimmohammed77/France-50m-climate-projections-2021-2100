This repository contains the R and Python scripts used to generate high-resolution (50 m) monthly climate projections for Metropolitan France (2021–2100) from CMIP6 models.
The workflow statistically downscales CHELSA V2.1 (~1 km) climate anomalies using the DIGITALIS V3 (50 m) topoclimatic baseline.

These scripts accompany the data descriptor paper:

Brahim, Piedallu & Serra-Diaz (2025)
Downscaling climate change projections to topoclimates across metropolitan France (2021–2100)

Contents

Model evaluation (R): comparison of CMIP6 models to Météo-France observations
CHELSA pre-processing (R): extraction, mosaicking, reprojection
Anomaly calculation (R): monthly deviations relative to 1981–2010
Tiles downscaling (R): Gaussian Kernel interpolation to 50 m
Mosaicking & climatology integration (Python): assembling final national rasters
The scripts reproduce the entire downscaling chain that produced the 23,040 GeoTIFF files included in the final dataset.

Requirements
R (≥ 4.3)

Required packages:

sf
terra
raster
sp
dplyr
tidyverse
stringr
parallel
lubridate
Rcpp

Python (≥ 3.10)

Required packages:

rasterio
numpy
pyproj
multiprocessing
glob
pathlib

Output Format

Final rasters follow the structure:

<variable>_<year>_<month>_<model>_<scenario>.tif


Projection: EPSG:2154 (Lambert-93)

Resolution: 50 m

Data type: int32

Values stored ×100 (divide by 100 to retrieve °C or mm)

Variables:

 tas = mean temperature
 tasmax = maximum temperature
 tasmin = minimum temperature
 pr = precipitation

Citation

If you use these scripts, please cite:

Brahim, M., Piedallu, C., Serra-Diaz, J.M. (2025).
Downscaling climate change projections to topoclimates across metropolitan France (2021–2100).

Author Contributions (CRediT)

 Mohammed Brahim – Conceptualization, Methodology, Software, Data Curation, Analysis, Visualization, Writing – Original Draft
 Christian Piedallu – Supervision, Methodology, Resources, Review & Editing
 Josep M. Serra-Diaz – Conceptualization, Funding Acquisition, Supervision, Review & Editing
