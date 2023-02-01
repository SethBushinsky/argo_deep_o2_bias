# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
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

# Check for a glodap_offsets_plots directory, create if it does not exist
offset_dir = output_dir + 'glodap_offset_plots/'
if not os.path.isdir(offset_dir):
    os.mkdir(offset_dir)
# -

# load glodap and float offsets and plot

#load saved glodap offsets 
glodap_offsets = xr.load_dataset(output_dir+'glodap_offsets.nc')

# plot histograms of offsets for each main float with crossovers

#group by main float wmo
offsets_g = glodap_offsets.groupby(glodap_offsets.main_float_wmo)

# +
#loop through offsets from each float
fig = plt.figure(figsize=(16,16))

for n,g in offsets_g:
    ncross = len(g.DOXY_ADJUSTED_offset)
    
    #if doxy_adjusted_offset is nan, skip
    if np.all(np.isnan(g.DOXY_ADJUSTED_offset.values)):
        continue
        
    #add crossover location map  
    axn = plt.subplot(4,4,1)
    axn.plot(g.main_float_longitude,g.main_float_latitude,'g+',markersize=10,label='Current float')
    #glodap
    glodap_lon = xr.where(g.glodap_longitude>180,g.glodap_longitude-360.,g.glodap_longitude)
    axn.plot(glodap_lon,g.glodap_latitude,'ro',label = 'Glodap',markersize=10)
    axn.legend()
    axn.set_title('WMO: %d, N crossovers: %d' % (g.main_float_wmo.values[0],ncross))
    
    #time
    axn = plt.subplot(4,4,(2,4))
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset,'rx',label='float')
    axn.plot(g.glodap_datetime,g.DOXY_ADJUSTED_offset,'bx',label='GDAP')
    axn.axhline(y=0, color='k', linestyle='--')
    axn.legend()
    axn.set_title('O2 Offsets vs date')
    
    #plot offset histograms
    axn = plt.subplot(4,5,6)
    g_plot = g.TEMP_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('TEMP %.2f +/- %.2f' % (g_mean,g_std))
    
    axn = plt.subplot(4,5,7)
    g_plot = g.PSAL_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('PSAL %.2f +/- %.2f' % (g_mean,g_std))
    
    axn = plt.subplot(4,5,8)
    g_plot = g.DOXY_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('DOXY %.2f +/- %.2f' % (g_mean,g_std))
    
    axn = plt.subplot(4,5,9)
    g_plot = g.NITRATE_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('NO3 %.2f +/- %.2f' % (g_mean,g_std))
 
    axn = plt.subplot(4,5,10)
    g_plot = g.NITRATE_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('PH %.2f +/- %.2f' % (g_mean,g_std))
    

    #O2 vs pressure 
    axn = plt.subplot(4,4,9)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_float,'bx')
    axn.set_title('Float pres vs DOXY offset')
    plt.gca().invert_yaxis()
    axn.set_ylabel('pres')
    
    axn = plt.subplot(4,4,10)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_glodap,'bx')
    axn.set_title('GDAP pres vs DOXY offset')
    plt.gca().invert_yaxis()
    
    #vs density 
    axn = plt.subplot(4,4,11)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PDENS_float,'bx')
    axn.set_title('Float Pdens vs DOXY offset')
    axn.set_ylabel('dens')
    
    axn = plt.subplot(4,4,12)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PDENS_glodap,'bx')
    axn.set_title('GDAP Pdens vs DOXY offset')
    
    #O2 offset vs o2 concentration
    axn = plt.subplot(4,4,13)
    axn.plot(g.DOXY_ADJUSTED_offset,g.DOXY_ADJUSTED_float,'bx')
    axn.set_title('Float O2 vs DOXY offset')
    axn.set_ylabel('O2 concentration')
    
    #T-S coloured by O2 offset
    axn = plt.subplot(4,4,14)
    axn.scatter(g.PSAL_ADJUSTED_float,g.TEMP_ADJUSTED_float,c=g.DOXY_ADJUSTED_offset)
    axn.set_title('float T-S')
    
    axn = plt.subplot(4,4,15)
    cn = axn.scatter(g.PSAL_ADJUSTED_glodap,g.TEMP_ADJUSTED_glodap,c=g.DOXY_ADJUSTED_offset)
    axn.set_title('GDAP T-S')
    
    cx = fig.add_axes([0.72,0.12,0.02,0.17])
    plt.colorbar(cn,cax=cx,label='DOXY OFFSET')

    
    plt.savefig(offset_dir+ 'individual_floats/' + str(g.main_float_wmo.values[0])+'_v_glodap.png')
    plt.clf()
  
# -

# load and process float - glodap offset metadata

# +
#list of DOXY_ADJUSTED SCIENTIFIC_CALIB_COMMENT substrings to group together for bias analysis

#no calibration bad data
bad_cal_list = ['Sensor issue','out of order','Bad data','Biofouling','unadjustable']

#no calibration, reason unspecified
no_cal_list = ['no adjustment','No QC','none','not applicable']

#air cal
air_cal_list = ['DOXY_ADJUSTED corrected using co','SVU Foil','Partial pressure','Bittig',
                'Adjusted with SAGE02 using co','Adjusted with SAGE02 with in-air',
                'PPOX converted from DOXY','G determined from float measure']
#no air cal
noair_cal_surf_list = ['DOXY_ADJUSTED is computed from','DOXY_ADJUSTED is estimated from',
                       'Adjustment done on PPOX_DOXY;Tem','Polynomial calibration','Adjustment via comparison of',
                      'optode multi calibration','RT adjustment by comparison','adjusted with median',
                      'Adjusted with annual WOA','Adjusted on WOA monthly','Adjusted PPOX at zero for',
                      'G determined by surface']

noair_cal_subsurf_list = ['1-point multiplicative corr']

noair_cal_funcofdoxy_list = ['Percent saturation corrected as','DOXY_ADJUSTED computed using Ste',
                        'Oxygen concentration corrected']

noair_cal_unspec_list = ['DOXY_ADJUSTED corrected based','Adjust with WOA monthly','GAIN determined from WOA2013',
                        'Adjusted with WOA climatology','Adjusted with SAGEO2 based on WO',
                         'Adjusted with SAGEO2 on WOA','Adjusted with WOA 2018','Takeshita and all, 2013']

noair_cal_withdrift_list = ['Adjustment done on PPOX_DOXY;Tem'] #this is incomplete

noair_cal_combined_list = noair_cal_surf_list+noair_cal_subsurf_list+noair_cal_funcofdoxy_list+noair_cal_unspec_list

# +
#initialise metadata DataArrays
glodap_offsets['o2_calib_comment'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['o2_calib_equation'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['o2_calib_coeff'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['o2_calib_group'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['o2_calib_drift'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['o2_calib_date'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['project_name'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['plat_type'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)
glodap_offsets['data_centre'] = xr.DataArray(dims=['N_CROSSOVERS'],coords={'N_CROSSOVERS':range(glodap_offsets.main_float_wmo.shape[0])}).astype(str)

num_crossovers = glodap_offsets.main_float_wmo.size

# I think this reloads the Sprof file for each individual crossover, this seems inefficient
for i,g in enumerate(glodap_offsets.main_float_wmo):
    
    #find full float file matching offset
    fn = argo_path + str(g.values) + '_Sprof.nc'
    float_og = xr.open_dataset(fn)
    
    #retrieve calibration info
    #get single profile
    p = glodap_offsets.main_float_profile[i].values
    p_index = np.flatnonzero(float_og.CYCLE_NUMBER==p)
    prof_og = float_og.isel(N_PROF=p_index)
    
    #retrieve calibration information
    #need to add option for multiple calibrations- defaul|t to use most recent only for now?
    #O2 calibration
    #if len(float_og.N_CALIB) >1 :
        #print(str(g.values)+' has more than 1 calibration')
    
    no_cal = []
    if 'DOXY_ADJUSTED' in float_og.keys():
        cal_str = prof_og.STATION_PARAMETERS.values.astype(str)[0]
        o2_ind = [idx for idx, s in enumerate(cal_str) if 'DOXY' in s]
        if len(o2_ind):
            o2_cal = prof_og.SCIENTIFIC_CALIB_COMMENT.values[0][-1][o2_ind[0]].decode("utf-8")
            o2_eq = prof_og.SCIENTIFIC_CALIB_EQUATION.values[0][-1][o2_ind[0]].decode("utf-8")
            o2_co = prof_og.SCIENTIFIC_CALIB_COEFFICIENT.values[0][-1][o2_ind[0]].decode("utf-8")
            if not prof_og.SCIENTIFIC_CALIB_DATE.values[0][-1][o2_ind[0]].decode("utf-8"):
                #print('no calibration info')
                o2_calib_date = 'nan'
                no_cal.append(glodap_offsets.main_float_wmo.values)
            else:
                o2_calib_date = prof_og.SCIENTIFIC_CALIB_DATE.values[0][-1][o2_ind[0]].decode("utf-8")
        else:
            o2_cal = 'empty O2 calib comment'
            o2_eq = 'none'
            o2_co = 'none'
            o2_calib_date = 'none'
    else:
        o2_cal = 'no O2'
        o2_eq = 'none'
        o2_co = 'none'
        o2_calib_date = 'none'
        
    
    #group together similar calib comments into one larger group
    if any(substring in o2_cal for substring in bad_cal_list):
        o2_group = 'bad'
    elif any(substring in o2_cal for substring in no_cal_list):
        o2_group = 'no cal'
    elif any(substring in o2_cal for substring in air_cal_list):
        o2_group = 'air cal'
    elif any(substring in o2_cal for substring in noair_cal_surf_list):
        o2_group = 'no air surface cal'
    elif any(substring in o2_cal for substring in noair_cal_subsurf_list):
        o2_group = 'no air subsurface cal'
    elif any(substring in o2_cal for substring in noair_cal_funcofdoxy_list):
        o2_group = 'no air function of DOXY cal'
    elif any(substring in o2_cal for substring in noair_cal_unspec_list):
        o2_group = 'no air unspecified'
    else:
        o2_group = 'ungrouped'
    
        
    #separate variable for drift/no drift
    if any(substring in o2_cal for substring in noair_cal_withdrift_list):
        o2_group_drift = 'drift corrected'
    else:
        o2_group_drift = 'no drift corrected'
    
    #add metadata to glodap_offset Dataset
    glodap_offsets.o2_calib_comment[i] = o2_cal
    glodap_offsets.o2_calib_equation[i] = o2_eq
    glodap_offsets.o2_calib_coeff[i] = o2_co
    glodap_offsets.o2_calib_group[i] = o2_group
    glodap_offsets.o2_calib_drift[i] = o2_group_drift
    glodap_offsets.o2_calib_date[i] = o2_calib_date
    glodap_offsets.project_name[i] = prof_og.PROJECT_NAME.values.astype(str)[0]
    glodap_offsets.plat_type[i] = prof_og.PLATFORM_TYPE.values.astype(str)[0]
    glodap_offsets.data_centre[i] = prof_og.DATA_CENTRE.values.astype(str)[0]

    print(str(i) + ' out of ' + str(num_crossovers))
#save offsets with cal info
glodap_offsets.to_netcdf(output_dir+'glodap_offsets_withcalibration.nc')
    
# -

# Group and plot offsets

#load saved offsets
glodap_offsets = xr.load_dataset(output_dir+'glodap_offsets_withcalibration.nc')

glodap_offsets.main_float_wmo

#create meta groups based on calibration groups (air cal, no air cal, no cal)
g = glodap_offsets.o2_calib_group.copy(deep=True)
g = g.where(glodap_offsets.o2_calib_group == 'air cal','not air')
glodap_offsets['o2_calib_air'] = xr.where(glodap_offsets.o2_calib_group=='bad','no cal',g)

# +
#extract gain from calibration coefficients
# -

# Group offsets by chosen float metadata and plot histograms

# +
glodap_offsets['o2_calib_air_key'] = xr.where(glodap_offsets.o2_calib_air=='air cal',1,0)

glodap_offsets['o2_not_air_key'] = xr.where(glodap_offsets.o2_calib_air=='not air',-1,0)
glodap_offsets['o2_cal_key'] = glodap_offsets['o2_calib_air_key'] + glodap_offsets['o2_not_air_key']
# -

plt.hist(glodap_offsets.o2_cal_key)

glodap_offsets.groupby("main_float_wmo")[3901668].o2_calib_comment

wmo_group = glodap_offsets.groupby('main_float_wmo')
list(glodap_offsets.groupby('main_float_wmo'))

glodap_offsets.groupby("main_float_wmo")[1900722]

wmo_means = glodap_offsets.groupby('main_float_wmo').mean(...)

# +
a = (wmo_means.where(wmo_means.o2_cal_key>0.6))
print(a.DOXY_ADJUSTED_offset.median())
print(a.DOXY_ADJUSTED_offset.mean())

#np.isnan(a.DOXY_ADJUSTED_offset).sum(1)
print(a.DOXY_ADJUSTED_offset.count())
plt.hist(a.DOXY_ADJUSTED_offset, bins=np.linspace(-20, 20, 41),)
plt.title('Air Cal Offsets, first averaged by float')
plt.grid()

plt.savefig(offset_dir + 'Glodap_offsets_doxy_air_cal_averaged_by_float.png')

# -

fig = plt.figure(figsize=(20,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
sct = plt.scatter(a.main_float_longitude, 
                a.main_float_latitude, 
                c=a.DOXY_ADJUSTED_offset,cmap='RdBu_r',s=40,vmin=-6,vmax=6)
cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
cbar.set_label('O2 offset', labelpad=15, fontsize=14)
#plt.scatter(glodap_offsets.glodap_longitude,glodap_offsets.glodap_longitude,s=4)
plt.savefig(offset_dir+ 'map_o2_offsets_air_only.png')
plt.show()

plt.hist(wmo_means['DOXY_ADJUSTED_offset'], bins=np.linspace(-20, 20, 61),label=str(n))
print(np.around(wmo_means['DOXY_ADJUSTED_offset'].median().values, decimals=2))
plt.title('All Offsets, first averaged by float')

# +
# which metadata variable to group by
group_variable = 'o2_calib_air'


#iterate through groups to plot offsets by group
offsets_g = glodap_offsets.groupby(glodap_offsets[group_variable])

plt.figure(figsize=(16,10))
for n,group in offsets_g:
    print(n)
    print(group['DOXY_ADJUSTED_offset'].shape[0])
    #print equations associated with comment
    print(np.unique(group['main_float_wmo'].values))
    #plot histogram
    plt.hist(group['DOXY_ADJUSTED_offset'], bins=np.linspace(-60, 60, 121),label=str(n))
    grp_mean = group['DOXY_ADJUSTED_offset'].mean()
    grp_std = group['DOXY_ADJUSTED_offset'].std()
    grp_med = group['DOXY_ADJUSTED_offset'].median()

    #calc mean/median values
    plt.grid()
    plt.xlabel('DOXY Offset')
    plt.title('mean: ' + str(grp_mean.values) + ' \pm ' + str(grp_std.values) + '; median: ' + str(grp_med.values))
    plt.savefig(offset_dir + 'Glodap_offsets_doxy_'+group_variable+'_'+str(n)+'.png')
    plt.clf()
    


# +
#plot histogram on same figure
plt.figure(figsize=(16,10))
for n,group in offsets_g:
    print(n)
    
    #calc mean values
    nmean = np.around(group['DOXY_ADJUSTED_offset'].mean().values, decimals=2)
    nmedian = np.around(group['DOXY_ADJUSTED_offset'].median().values, decimals=2)
    print(nmean)
    print(nmedian)
    
    #plot histogram
    plt.hist(group['DOXY_ADJUSTED_offset'], bins=np.linspace(-60, 60, 121),
             alpha=0.3,label=str(n)+', mean='+str(nmean))


plt.xlabel('DOXY Offset')
plt.legend()
plt.savefig(offset_dir + 'Glodap_offsets_doxy_all_'+group_variable+'.png')
# -

# Plot histograms of all global glodap offsets combined

plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.DOXY_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
plt.xlabel('DOXY Offset')
plt.savefig(offset_dir + 'Glodap_offsets_doxy_plus_minus_400.png')


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
plt.savefig(offset_dir+ 'map_o2_offsets_100.png')
plt.show()

# +
#nitrate all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.NITRATE_ADJUSTED_offset, bins=np.linspace(-4, 4, 30))
nmean = np.around(glodap_offsets['NITRATE_ADJUSTED_offset'].mean().values, decimals=2)
nmedian = np.around(glodap_offsets['NITRATE_ADJUSTED_offset'].median().values, decimals=2)
plt.grid()
    
plt.title('Mean: ' + str(nmean))
plt.xlabel('NITRATE Offset')
plt.savefig(offset_dir + 'Glodap_offsets_nitrate_plus_minus_400.png')


# +
#nitrate all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.DIC_offset, bins=np.linspace(-40, 40, 30))
nmean = np.around(glodap_offsets['DIC_offset'].mean().values, decimals=2)
nmedian = np.around(glodap_offsets['DIC_offset'].median().values, decimals=2)
plt.grid()
    
plt.title('Mean: ' + str(nmean))
plt.xlabel('DIC Offset')
plt.savefig(offset_dir + 'Glodap_offsets_DIC.png')
# -

glodap_offsets

# plot histograms of offsets for each main float with crossovers

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
        plt.savefig(offset_dir+str(wmo)+'_v_float.png')
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
        plt.savefig(offset_dir+str(wmo)+'_v_glodap.png')
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
        plt.savefig(offset_dir+str(wmo)+'_v_float_and_glodap.png')
        plt.clf()
