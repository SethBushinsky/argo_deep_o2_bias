# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Compare bgc float data crossovers (float and GLODAP) and apply bias correction

# 1. Download and process float + GLODAP data
# 2. apply float bias corrections
# 3. calculate derivative variables (pH, TALK)
# 4. do float/glodap crossover matchups
# 5. plot results

import numpy as np
import glob, os
from pathlib import Path
from datetime import datetime, date, time
import pandas as pd
import xarray as xr
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import colorbar, colors
import cartopy
import cartopy.crs as ccrs
import PyCO2SYS as pyco2
import gsw
import float_data_processing as fl
import carbon_utils

# Create data directories

# +
# Set the paths
output_dir = 'output/'
data_dir = 'data/'

#check directories exist
if not os.path.isdir('output'):
    os.mkdir('output')
if not os.path.isdir('data'):
    os.mkdir('data')
# -

# Download datasets and read in as xarray

gdap = fl.get_glodap(data_dir, year = 2021)

# GLODAP data processing

gdap.G2longitude[gdap.G2longitude < 0.] = gdap.G2longitude[gdap.G2longitude < 0.] + 360.


# GLODAP derived variables: density, MLD and pH

# +
#set flagged data to NaN (is this needed? or masked array better?)
flagvars = ['G2salinity','G2oxygen','G2nitrate','G2tco2','G2talk','G2phts25p0']

for v in flagvars:
    flag = v+'f'
    naninds = gdap[flag]!=2
    gdap[v][naninds] = np.nan
# -

#need to group by cruise and station
#get unique stations
gdap_s = gdap.groupby(['G2cruise','G2station','G2cast'])

# +
for n, group in gdap_s:
    #need to only do if not all -9999!!!
    zmin = np.absolute(group.G2depth-10.).argmin()
    zminind = group.G2depth.index[zmin]
    mld_sigma = group.G2sigma0[zminind]+0.03 #density threshold value
    #look only below 10m
    mldind = group[zmin:].index[group.G2sigma0[zmin:]>=mld_sigma] #identify MLD 
    if len(mldind):
        group['MLD_sigma0'] = group.G2depth[mldind[0]]
        
#put groups back to original dataframe. pd.join?

# +
#pH from LIPHR --> need wrapper
# calculate LIPHR pH at Glodap points below 1480 m and above 2020m
Coordinates = [gdap.G2longitude gdap.G2latitude, gdap.G2pressure]
Measurements = [gdap.G2salinity, gdap.G2temperature, gdap.G2nitrate, gdap.G2oxygen]
MeasIDVec = [1, 7, 3, 6]

pHEstimates,UncertaintyEstimates,MinUncertaintyEquation = carbon_utils.LIPHR_matlab(Coordinates,
                                                                                    Measurements,
                                                                                    MeasIDVec, 
                                                                                    OAAdjustTF = False)                                  

gdap['pH_in_situ_total'] = pHEstimates
gdap.pH_in_situ_total[np.isnan(gdap.G2phts25p0)] = np.nan

# -

# gdap pH 25C - NEED TO UPDATE INPUTS HERE for SOCCOM version-> see Nancy's note in Slack
results = pyco2.sys(
    par1=2300., 
    par2=gdap.PH_IN_SITU_TOTAL_ADJUSTED,
    par1_type=1,
    par2_type=3,
    temperature=gdap.TEMP_ADJUSTED, 
    pressure=gdap.PRES_ADJUSTED, 
    salinity=gdap.PSAL_ADJUSTED, 
    temperature_out=25., #fixed 25C temperature
    pressure_out=gdap.PRES_ADJUSTED,
    opt_k_carbonic=10,
    opt_k_bisulfate=1, # this is the default
    opt_k_fluoride=1, # this is the default
    buffers_mode='auto',
)
#not sure what's happening here with constant inputs?
#[DATA,~,~]=CO2SYSSOCCOM_smb(2300.*ones(length(gdap_SO.TEMP_ADJUSTED),1), gdap_SO.PH_IN_SITU_TOTAL_ADJUSTED, ...
#    1,3, gdap_SO.PSAL_ADJUSTED, gdap_SO.TEMP_ADJUSTED, 25.*ones(length(gdap_SO.TEMP_ADJUSTED),1),...
#    gdap_SO.PRES_ADJUSTED, gdap_SO.PRES_ADJUSTED, zeros(length(gdap_SO.TEMP_ADJUSTED),1), zeros(length(gdap_SO.TEMP_ADJUSTED),1),1,10,3);
#gdap_SO.pH_25C_TOTAL = DATA(:,37);
#
#set pH to nan where there was no original pH data from GLODAP
#gdap_SO.pH_25C_TOTAL(isnan(gdap_SO.G2phts25p0))=nan;

# +
# change Glodap names to Argo names
#name_convert = {'G2longitude' 'LONGITUDE'; 'G2latitude', 'LATITUDE'; 'G2pressure', 'PRES_ADJUSTED'; 'G2temperature', 'TEMP_ADJUSTED'; 'G2salinity', 'PSAL_ADJUSTED';...
#   'G2oxygen' 'DOXY_ADJUSTED'; 'G2nitrate' 'NITRATE_ADJUSTED';'G2tco2' 'DIC' ; 'G2talk' 'TALK_LIAR' ; 'G2MLD' 'MLD'; 'G2o2sat' 'o2sat' ; 'G2PTMP' 'PTMP';'pH_in_situ_total' 'PH_IN_SITU_TOTAL_ADJUSTED'};
# -


# Apply float bias corrections - this is the meaty part



# Calculate float TALK



# Calculate float pH and apply pH bias correction



# Compare float/GLODAP crossovers




