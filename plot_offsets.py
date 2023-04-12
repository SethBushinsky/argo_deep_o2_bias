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

# ## Plotting script for analysing offsets and applying bias adjustments

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
    os.mkdir(offset_dir+'individual_floats/')
# -

# ### 1. Load offsets and look at offset summary plots for each float
# This step helps to 1) identify patterns in float offsets in particular floats so that we can go back and refine offset criteria and 2) make final decisions about whether a float should be adjusted or not.
#
# First load saved glodap offsets

glodap_offsets = xr.load_dataset(output_dir+'glodap_offsets.nc')

# plot histograms of offsets for each main float with all crossovers

#group by main float wmo
offsets_g = glodap_offsets.groupby(glodap_offsets.main_float_wmo)

offsets_g

# +
#loop through each float
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

# ### 2. Combine metadata (calibration, instrument etc) with offsets so offsets can be sorted and plotted using different criteria
#
# For the following analyis we only need the mean offset for each float, so can take the mean by each float wmo first.

glodap_offsets_mean = offsets_g.mean(skipna='True')

# We also need to create broad groupings of calibration information provided in SCIENTIFIC_CALIB_COMMENT to enable easy sorting and plotting

# +
#list of DOXY_ADJUSTED SCIENTIFIC_CALIB_COMMENT substrings to group together for bias analysis

#no calibration bad data
bad_cal_list = ['Sensor issue','out of order','Bad data; not adjustable','Biofouling','unadjustable']

#no calibration, reason unspecified
no_cal_list = ['no adjustment','No QC','none','not applicable']

#blank cal
blank_cal = ['                             ']

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
# -

# Now we can loop through each WMO and combine the calibration + sensor information with the mean offset for each float (note that in some cases individual profiles from a float may have different calibration comments and coeffients, but each float should be able to be grouped by O2 calibration)

# +
#initialise metadata DataArrays
varnames = ['o2_calib_comment','o2_calib_equation','o2_calib_coeff','o2_calib_group','o2_calib_air_group',
           'o2_calib_drift','project_name','plat_type','data_centre']
coord= range(glodap_offsets_mean.main_float_wmo.shape[0])
for v in varnames:
    glodap_offsets_mean[v] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

#also create groups by sensor 
glodap_offsets_mean['pH_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
glodap_offsets_mean['pH_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
glodap_offsets_mean['nitrate_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
glodap_offsets_mean['nitrate_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
glodap_offsets_mean['DOXY_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
glodap_offsets_mean['DOXY_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

#group for floats that have some under ice profiles with some air calibrated some not
glodap_offsets_mean['ice_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

num_crossovers = glodap_offsets_mean.main_float_wmo.size

# loop through each wmo in offsets
for n,w in enumerate(glodap_offsets_mean.main_float_wmo):
    
    #find full float file matching offset
    wmo_n = w.values
    fn = argo_path + str(wmo_n) + '_Sprof.nc'
    float_og = xr.open_dataset(fn)
    
    #also load meta file  for same float
    fn1 = argo_path + str(wmo_n) + '_meta.nc'
    meta_og = xr.open_dataset(fn1)

    #option for multiple calibrations at different times- default to use most recent only for now?
    #O2 calibration
    #if len(float_og.N_CALIB) >1 :
        #print(str(g.values)+' has more than 1 calibration')
    
    no_cal = []
    if 'DOXY_ADJUSTED' in float_og.keys():
        is_doxy = 'DOXY'
        
        #get O2 sensor type
        #check sensor index
        sensors = meta_og.SENSOR.values
        sind = [idx for idx, s in enumerate(sensors) if 'DOXY' in s.decode("utf-8")] 
        #grab sensor model for index
        if len(sind):
            o2_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
        else:
            o2_sensor = 'not listed'

        cal_str = float_og.STATION_PARAMETERS.values.astype(str)[0]
        o2_ind = [idx for idx, s in enumerate(cal_str) if 'DOXY' in s]
        if len(o2_ind):
 
            
            #get calibration info
            o2_cal_full = float_og.SCIENTIFIC_CALIB_COMMENT.values[:,-1,o2_ind[0]]
            o2_cal_date = float_og.SCIENTIFIC_CALIB_DATE.values[:,-1,o2_ind[0]]
            o2_cal_unique = np.unique(o2_cal_full)
            
            #check if any unique cal comments are bad/no cal
            isbad = []
            for i in range(len(o2_cal_unique)):
                o2_cal_i = o2_cal_unique[i].decode("utf-8")
                if any(substring in o2_cal_i for substring in bad_cal_list) or \
                   any(substring in o2_cal_i for substring in no_cal_list) or \
                   any(substring in o2_cal_i for substring in blank_cal):
                    isbad.append(i)
            #remove bad indices
            o2_cal_unique = np.delete(o2_cal_unique,isbad)

            #check if more than one unique calibration comment remaining after excluding bad data
            if len(o2_cal_unique)>1:
                print('multiple cal comments')  
                print(wmo_n)
                #how to integrate multiple comments? Do any contain different groupings?
                o2_cal = float_og.SCIENTIFIC_CALIB_COMMENT.values[0,-1,o2_ind[0]].decode("utf-8")
                if np.any(float_og.POSITION_QC.values==8):
                    under_ice = 'some air cal, some under ice'
                    print('under ice')
                else:
                    under_ice = 'no ice'
                    
            else:
                o2_cal = float_og.SCIENTIFIC_CALIB_COMMENT.values[0,-1,o2_ind[0]].decode("utf-8")
            o2_eq = float_og.SCIENTIFIC_CALIB_EQUATION.values[0,-1,o2_ind[0]].decode("utf-8")
            o2_co = float_og.SCIENTIFIC_CALIB_COEFFICIENT.values[0,-1,o2_ind[0]].decode("utf-8")
            under_ice = 'no ice'
        else:
            o2_cal = 'empty O2 calib comment'
            o2_eq = 'none'
            o2_co = 'none'

    else:
        is_doxy = 'no DOXY'
        o2_cal = 'no O2'
        o2_eq = 'none'
        o2_co = 'none'
        o2_sensor = 'none'
    

    ######## group calibration comments into categories
    
    #o2 calib comment groups
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
    
    
    #group O2 air cal and no air cal meta groups
    if any(substring in o2_cal for substring in air_cal_list):
        o2_air_group = 'air cal'
    elif any(substring in o2_cal for substring in noair_cal_combined_list):
        o2_air_group = 'no air cal'
    else:
        o2_air_group = 'no cal/bad'
        
    #separate groups for drift/no drift
    if any(substring in o2_cal for substring in noair_cal_withdrift_list):
        o2_group_drift = 'drift corrected'
    else:
        o2_group_drift = 'no drift corrected'
    
    if 'PH_IN_SITU_TOTAL_ADJUSTED' in float_og.keys():
        is_pH = 'pH'
        #get pH sensor type
        #check sensor index
        sensors = meta_og.SENSOR.values
        sind = [idx for idx, s in enumerate(sensors) if 'PH' in s.decode("utf-8")] 
        #grab sensor model for index
        pH_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
    else:
        is_pH = 'no pH'
        pH_sensor = 'none'
    
    if 'NITRATE_ADJUSTED' in float_og.keys():
        is_nitrate = 'nitrate'
        #get nitrate sensor type
        #check sensor index
        sensors = meta_og.SENSOR.values
        sind = [idx for idx, s in enumerate(sensors) if 'NITRATE' in s.decode("utf-8")] 
        #grab sensor model for index
        nitrate_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
    else:
        is_nitrate = 'no nitrate'
        nitrate_sensor = 'none'
    

    ###### get sensor information from meta files

    
    
    ##################################################
    #add metadata to glodap_offset Dataset
    glodap_offsets_mean.o2_calib_comment[n] = o2_cal
    glodap_offsets_mean.o2_calib_equation[n] = o2_eq
    glodap_offsets_mean.o2_calib_coeff[n] = o2_co
    glodap_offsets_mean.o2_calib_group[n] = o2_group
    glodap_offsets_mean.o2_calib_air_group[n] = o2_air_group
    glodap_offsets_mean.o2_calib_drift[n] = o2_group_drift
    glodap_offsets_mean.project_name[n] = float_og.PROJECT_NAME.values.astype(str)[0]
    glodap_offsets_mean.plat_type[n] = float_og.PLATFORM_TYPE.values.astype(str)[0]
    glodap_offsets_mean.data_centre[n] = float_og.DATA_CENTRE.values.astype(str)[0]
    glodap_offsets_mean.DOXY_group[n] = is_doxy
    glodap_offsets_mean.DOXY_sensor[n] = o2_sensor
    glodap_offsets_mean.pH_group[n] = is_pH
    glodap_offsets_mean.pH_sensor[n] = pH_sensor
    glodap_offsets_mean.nitrate_group[n] = is_nitrate
    glodap_offsets_mean.nitrate_sensor[n] = nitrate_sensor
    glodap_offsets_mean.ice_group[n] = under_ice
    
    print(str(n) + ' out of ' + str(num_crossovers))
#save offsets with cal info
glodap_offsets_mean.to_netcdf(output_dir+'glodap_offsets_floatmean_withcalibration.nc')
    
# -

# ### 3. Example: group offsets by chosen float variables/metadata and plot histograms
#
# In this example, float offsets are grouped by a) floats with/without pH and b) DOXY air calibrated/not air calibrated. Parameters a) and b) can be substituted by any other dataarray in glodap_offsets_mean
#
# First define parameters a) and b) to group offsets

glodap_offsets_mean


parameter_a = 'pH_group'
parameter_b = 'o2_calib_air_group'
#can pandas group by three different things?

# Now group mean offsets by parameters a and b and print statistics

#convert to pandas
glodap_offsets_p = glodap_offsets_mean.to_dataframe()
#do groupby paramater a and b
offsets_g = glodap_offsets_p.groupby([parameter_a,parameter_b])

# Plot pdfs for each group

glodap_offsets_mean

n

# +
#loop through groups and plot histograms and mean/median values of offsets
plt.figure(figsize=(16,10))
for n,group in offsets_g:
    #if n[1] == 'no cal/bad':
    #    continue
    print(n)
    
    #calc mean values
    nmean = np.around(group['DOXY_ADJUSTED_offset'].mean(), decimals=2)
    nmedian = np.around(group['DOXY_ADJUSTED_offset'].median(), decimals=2)
    ncount = group['DOXY_ADJUSTED_offset'].count()
    print(nmean)
    print(nmedian)
    print(ncount)

    #plot histogram
    plt.hist(group['DOXY_ADJUSTED_offset'], bins=np.linspace(-60, 60, 121),
             alpha=0.3,label=str(n)+', mean='+str(nmean) + ', n='+str(ncount))



plt.xlabel('DOXY Offset')
plt.legend()
plt.savefig(offset_dir + 'Glodap_offsets_doxy_all_'+parameter_a+parameter_b+'_v3.png')
# -

# ### 4. Example: map offsets for single float sub-group
#
# In this example plot histogram of float-mean DOXY_ADJUSTED offsets for air-calibrated floats only

parameter = 'o2_calib_air_group'
group_name = 'air cal'

# group offsets by single parameter

offsets_g1 = glodap_offsets_mean.groupby(parameter)

# now map offsets for single group, in this example floats with air calibration only - this should be profile by profile offsets not float by float??

for n, group in offsets_g1:
    if n == group_name:
        fig = plt.figure(figsize=(20,12))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_global()
        sct = plt.scatter(group.main_float_longitude, 
                group.main_float_latitude, 
                c=group.DOXY_ADJUSTED_offset,cmap='RdBu_r',s=40,vmin=-10,vmax=10)
        cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
        cbar.set_label('O2 offset', labelpad=15, fontsize=14)
        #plt.scatter(glodap_offsets.glodap_longitude,glodap_offsets.glodap_longitude,s=4)
        plt.savefig(offset_dir+ 'map_o2_offsets_air_only_v3.png')
        plt.show()

# ### 5. Example: histogram of all global DOXY offsets, mean of each float

plt.figure()
plt.hist(glodap_offsets_mean['DOXY_ADJUSTED_offset'], bins=np.linspace(-20, 20, 61),label=str(n))
print(np.around(glodap_offsets_mean['DOXY_ADJUSTED_offset'].median().values, decimals=2))
plt.title('All Offsets, first averaged by float')

# ### Other code

plt.figure(figsize=(20,12))
plt.hist(glodap_offsets_mean.DOXY_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
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
#DIC all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.DIC_offset, bins=np.linspace(-40, 40, 30))
nmean = np.around(glodap_offsets['DIC_offset'].mean().values, decimals=2)
nmedian = np.around(glodap_offsets['DIC_offset'].median().values, decimals=2)
plt.grid()
    
plt.title('Mean: ' + str(nmean))
plt.xlabel('DIC Offset')
plt.savefig(offset_dir + 'Glodap_offsets_DIC.png')
