# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# Plotting script for deep O2 bias manuscript
#
# 2022_07_21 Seth Bushinsky - will try to piggy-back off code being writting by Veronica Tamsitt & the rest of the Argo BGC bias adjustment project (Bushinsky, Nachod, Fassbender, Williams)
#

import numpy as np
import pandas as pd
import xarray as xr
import glob, os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, date, time
import time


# Function for converting date to year fraction (taken from stack overflow: https://stackoverflow.com/questions/6451655/how-to-convert-python-datetime-dates-to-decimal-float-years)

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction



# +
# read in a user-created text file to point to local directories to avoid having to change this every time 
# we update code
lines=[]
with open('path_file.txt') as f:
    lines = f.readlines()
    
count = 0
for line in lines:
    count += 1
    index = line.find("=")
    #print(f'line {count}: {line}')
    #print(index)
    #print(line[0:index])
    line = line.rstrip()
    if line[0:index].find("argo")>=0:
        argo_path=line[index+1:]
    elif line[0:index].find("liar")>=0:
        liar_dir=line[index+1:]
    elif line[0:index].find("matlab")>=0:
        matlab_dir=line[index+1:]


# +
output_dir = 'figures_o2_bias/'

#check directories exist
if not os.path.isdir('figures_o2_bias'):
    os.mkdir('figures_o2_bias')
# -

# 1. Read in all float Sprof files from argo path

# +
argolist = []
for file in os.listdir(argo_path):
    if file.endswith('Sprof.nc'):
        argolist.append(file)

last = 0
# -

# 2. For each file, find all profiles with valid O2 data, save out lat/lon/date into an array

#wmo_list = list()
for n in range(last, len(argolist)):
    print(f' {n}' ' File: ' + argolist[n]) 

    #n = 0
    argo_n = xr.load_dataset(argo_path+argolist[n])
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))
    
    var_list = list(argo_n.data_vars)
    if "DOXY_ADJUSTED" not in var_list:
        continue
            
    doxy_trimmed = argo_n.DOXY_ADJUSTED.where(~np.isnan(argo_n.DOXY_ADJUSTED), drop=True)
    wmo_n = argo_n.PLATFORM_NUMBER.values.astype(int)[0]

    prof_loc = xr.Dataset()
    prof_loc['wmo']=(['N_PROF'],np.repeat(wmo_n,len(doxy_trimmed)))

    prof_loc['LATITUDE'] = (['N_PROF'], doxy_trimmed.LATITUDE.data)
    prof_loc['LONGITUDE'] = (['N_PROF'], doxy_trimmed.LONGITUDE.data)
    prof_loc['juld'] = (['N_PROF'],doxy_trimmed.JULD.data)
    #prof_loc['profile'] = (['N_PROF'],doxy_trimmed.CYCLE_NUMBER.data)
    #argo_n
    # append all files into one long xarray
    if n == 0:
        argo_all = prof_loc
    else:
        argo_all = xr.concat([argo_all,prof_loc], 'N_PROF')

    last = n + 1
argo_all

# Plotting

# plotting all data
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
plt.plot(argo_all.LONGITUDE, argo_all.LATITUDE, linestyle='none', marker='.', markersize=1)
plt.show()

# +
### Note that this currently screws up Longitude averaging for floats that cross -180/180. Need to deal with that
argo_wmo = argo_all.groupby('wmo').mean()

avg_date = argo_all.juld.astype('int64').groupby(argo_all['wmo']).mean().astype('datetime64[ns]')

# calculate the year fraction for plotting
avg_date_dec_year = []
for i in avg_date.values:
    x = pd.Timestamp(i)
    y = toYearFraction(x)
    avg_date_dec_year.append(y)


# +
# Figure 1

fig = plt.figure(figsize=(20,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
sct = plt.scatter(x=argo_wmo.LONGITUDE, 
            y=argo_wmo.LATITUDE, 
            c=avg_date_dec_year,cmap='turbo',s=10)
cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
cbar.set_label('Year', labelpad=15, fontsize=14)
plt.savefig(output_dir+ 'Fig_1_O2_float_map_date.png')
plt.show()

