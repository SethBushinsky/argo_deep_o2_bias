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
float_offsets = xr.load_dataset(output_dir+'float_offsets.nc')

# Plot histograms of all global glodap offsets combined

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

#nitrate all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.NITRATE_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
plt.xlabel('NITRATE Offset')
plt.savefig(output_dir + 'Glodap_offsets_nitrate_plus_minus_400.png')


# Now plot histograms of offsets for each main float with crossovers

# +
#load saved argo_interp data 
argo_interp = xr.open_dataset(data_dir+'argo_interp_temp.nc')

#group by float wmo
argo_wmo = argo_interp.groupby('wmo')


#loop through each float
fig = plt.figure(figsize=(12,10))

for wmo, group in argo_wmo:
    #get index all glodap crossovers 
    g_inds = np.flatnonzero(gdap_offsets[len(var_list_plot)*2] == wmo)
    #get index of all float crossovers
    fl_inds = np.flatnonzero(float_offsets[len(var_list_plot)*2] == wmo)
    
    if len(g_inds)==0 and len(fl_inds)==0:
        continue
    elif len(g_inds)==0:
        #float crossover only
        #loop through each variable
        for idx, var in enumerate(var_list_plot):
            if len(float_offsets[idx*2]):
                f_plot = np.array(float_offsets[idx*2])[fl_inds]
                if not np.all(np.isnan(f_plot)):
                    axn = plt.subplot(3,3,idx+1)
                    axn.hist(f_plot,color='r',alpha=0.5)
                axn.set_title(var)
        
        #add crossover location map
        #main float positions
        fl_lon = np.array(float_offsets[17])[fl_inds]
        fl_lat = np.array(float_offsets[19])[fl_inds]
        g_lon = np.array(gdap_offsets[14])[g_inds]
        g_lat = np.array(gdap_offsets[16])[g_inds]
        
        axn = plt.subplot(3,2,6)
        axn.plot(group.LONGITUDE,group.LATITUDE,'bo',markersize=10,label='Current float')
        #test float positions
        axn.plot(fl_lon,fl_lat,'go',label = 'test floats',markersize=10)

        plt.grid(linestyle=':')
        plt.legend()
        plt.savefig(output_dir+str(wmo)+'_v_float.png')
        plt.clf()
        
    elif len(fl_inds)==0:
        #glodap crossover only
        #loop through each variable
        for idx, var in enumerate(var_list_plot):
            if len(gdap_offsets[idx*2]):
                axn = plt.subplot(3,3,idx+1)
                g_plot = np.array(gdap_offsets[idx*2])[g_inds]
                if not np.all(np.isnan(g_plot)):
                    axn.hist(g_plot,color='b',alpha=0.5)
                axn.set_title(var)
        
        #add crossover location map
        #main float positions
        fl_lon = np.array(float_offsets[17])[fl_inds]
        fl_lat = np.array(float_offsets[19])[fl_inds]
        g_lon = np.array(gdap_offsets[14])[g_inds]
        g_lat = np.array(gdap_offsets[16])[g_inds]
        
        axn = plt.subplot(3,2,6)
        axn.plot(group.LONGITUDE,group.LATITUDE,'bo',markersize=10,label='Current float')
        #glodap
        axn.plot(g_lon,g_lat,'mv',label = 'Glodap',markersize=10)
        
        plt.grid(linestyle=':')
        plt.legend()
        plt.savefig(output_dir+str(wmo)+'_v_glodap.png')
        plt.clf()
        
    else:
        #loop through each variable
        for idx, var in enumerate(var_list_plot):
            if len(float_offsets[idx*2]) and len(gdap_offsets[idx*2]):
                axn = plt.subplot(3,3,idx+1)
                #print(np.array(float_offsets[idx]))
                f_plot = np.array(float_offsets[idx*2])[fl_inds]
                g_plot = np.array(gdap_offsets[idx*2])[g_inds]
                if not np.all(np.isnan(g_plot)):
                    axn.hist(g_plot,color='b',alpha=0.5)
                if not np.all(np.isnan(f_plot)):
                    axn.hist(f_plot,color='r',alpha=0.5)
                axn.set_title(var)
        
        #add crossover location map
        #main float positions
        fl_lon = np.array(float_offsets[17])[fl_inds]
        fl_lat = np.array(float_offsets[19])[fl_inds]
        g_lon = np.array(gdap_offsets[14])[g_inds]
        g_lat = np.array(gdap_offsets[16])[g_inds]
        
        axn = plt.subplot(3,2,6)
        axn.plot(group.LONGITUDE,group.LATITUDE,'bo',markersize=10,label='Current float')
        #test float positions
        axn.plot(fl_lon,fl_lat,'go',label = 'test floats',markersize=10)
        #glodap
        axn.plot(g_lon,g_lat,'mv',label = 'Glodap',markersize=10)
        
        plt.grid(linestyle=':')
        plt.legend()
        plt.savefig(output_dir+str(wmo)+'_v_float_and_glodap.png')
        plt.clf()
