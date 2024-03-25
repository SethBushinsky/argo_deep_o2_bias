# float_bgc_synthesis_products


[![DOI](https://zenodo.org/badge/656396759.svg)](https://zenodo.org/doi/10.5281/zenodo.10866941)


## Download float data:
### float_download_sprof_meta.py
- Saves:
     - Sprof, meta files
## Load glodap, float data, perform derived calculations, interpolate, find and save crossovers
### float_bgc_bias_correction_parallel.ipynb
1. Calls:
  - carbon_utils.py
  - float_data_processing.py
  - argo_interp_and_crossover.py
  - outlier_filter_ESD_test.py
  - gdap_crossover_intermediate_script.py
2. Saves:
  - Derived files
  - Interpolated files
  - Glodap crossovers
## Calculate mean and trimmed (outliers removed) glodap crossovers for different pressure levels. Can run multiple files (for instance if I want to use different crossover criteria for different pressure ranges)
### glodap_crossover_approach_comparison_parallel_attempt.ipynb
1. Calls:
  - pressure_level_glodap_mean.py
  - Plot_offsets_crossover_plot_only.py
3. Saves:
  - Individual netcdf files for each pressure level
  - Glodap crossover plots for each file / pressure level (need to add two more loops)
  - Plots of offsets vs. pressure with different criteria 
## Figures for the manuscript
### deep_o2_analysis_plotting.ipynb
