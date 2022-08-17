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

# Plotting script for analysing offsets and applying bias adjustments

import numpy as np
import pandas as pd
import xarray as xr
import glob, os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, date, time
import time

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
        
# Set the paths
output_dir = 'output/'
data_dir = 'data/'
# -

# load glodap and float offsets and plot

#load saved glodap offsets
glodap_offsets = xr.load_dataset(output_dir+'glodap_offsets.nc')

plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.DOXY_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
plt.xlabel('DOXY Offset')
plt.savefig(output_dir + 'Glodap_offsets_doxy_plus_minus_400.png')


#check location of O2 offsets > 100
offset_100 = glodap_offsets.DOXY_ADJUSTED_offset.where(glodap_offsets.DOXY_ADJUSTED_offset<-100)


fig = plt.figure(figsize=(20,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
sct = plt.scatter(glodap_offsets.main_float_longitude, 
                glodap_offsets.main_float_latitude, 
                c=offset_100,cmap='viridis',s=10,vmin=-300,vmax=-100)
cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
cbar.set_label('O2 offset', labelpad=15, fontsize=14)
plt.scatter(glodap_offsets.glodap_longitude,glodap_offsets.glodap_longitude,s=4)
plt.savefig(output_dir+ 'map_o2_offsets_100.png')
plt.show()

plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.NITRATE_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
plt.xlabel('NITRATE Offset')
plt.savefig(output_dir + 'Glodap_offsets_nitrate_plus_minus_400.png')

