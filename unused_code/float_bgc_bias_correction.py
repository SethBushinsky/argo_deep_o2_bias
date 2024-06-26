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

# ## Apply BGC float bias correction and compare float/GLODAP crossovers

# Authors: Veronica Tamsitt (USF) et al...
#
# Adapted from MATLAB code written by Seth Bushinsky (UH)
#
# 1. Download and process GLODAP data
# 3. apply float bias corrections and calculate derivative variables (pH, TALK)
# 3. do float - glodap crossover comparison
# 4. do float - float crossover comparison
#
# Link to MATLAB LIAR/LIPHR code: https://github.com/BRCScienceProducts/LIRs
#

# import modules

import numpy as np
import glob, os
from pathlib import Path
from datetime import datetime, date, time
import pandas as pd
import xarray as xr
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import colorbar, colors
import cartopy
import cartopy.crs as ccrs
import PyCO2SYS as pyco2
import gsw
import float_data_processing as fl
import carbon_utils
import re
import time

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

# Check for a glodap_offsets_plots directory, create if it does not exist
offset_dir = output_dir + 'glodap_offset_plots/'
if not os.path.isdir(offset_dir):
    os.mkdir(offset_dir)

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

#add derived float file directory within argo_path
argo_path_derived = argo_path+'../derived/'
if not os.path.isdir(argo_path_derived):
    os.mkdir(argo_path_derived)
# -

# ### User inputs: 

# +
#pressure limits for interpolation
p_interp_min = 1450 #minimum pressure for float crossover comparison
p_interp_max = 2000 #maximum pressure for float crossover comparison
#pressure levels to interpolate to, every 1db
p_interp = np.arange(p_interp_min,p_interp_max+1)

#pressure limits for crossover comparison
p_compare_min = 1400
p_compare_max = 2100

#max density difference to store crossover
delta_dens = 0.005
# delta_dens = 0.05

#max spice difference to store crossover
delta_spice = 0.005
# delta_spice = 0.05

# max pressure difference to store crossover
delta_press = 100

#crossover distance range
dist = 100

#variables to do crossovers
var_list_plot = ['PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED','DOXY_ADJUSTED','NITRATE_ADJUSTED',
                 'DIC','pH_25C_TOTAL_ADJUSTED','PH_IN_SITU_TOTAL_ADJUSTED','PDENS']
# -

# ## 1. Download and process GLODAP data

# +
gdap = fl.get_glodap(data_dir, year = 2022)
gdap.G2longitude[gdap.G2longitude < 0.] = gdap.G2longitude[gdap.G2longitude < 0.] + 360.
#set flagged data to NaN (is this needed? or masked array better?)
flagvars = ['G2salinity','G2oxygen','G2nitrate','G2tco2','G2talk','G2phts25p0']

for v in flagvars:
    flag = v+'f'
    naninds = gdap[flag]!=2
    gdap[v][naninds] = np.nan

# GLODAP derived variables: density, MLD and pH

#calc potential density
gdap['sigma0_calculated'] = carbon_utils.sigma0(gdap.G2salinity.values,gdap.G2temperature.values,
                                  gdap.G2longitude.values,gdap.G2latitude.values,gdap.G2pressure.values)
#calculate spice
gdap['spice'] = carbon_utils.spiciness0(gdap.G2salinity.values,gdap.G2temperature.values,
                                  gdap.G2longitude.values,gdap.G2latitude.values,gdap.G2pressure.values)

#pH from LIPHR
# calculate LIPHR pH at Glodap points below 1480 m and above 2020m (V: where does the depth restriction come in?)
LIPHR_path = liar_dir
Coordinates = np.stack((gdap.G2longitude.values.flatten(), 
                        gdap.G2latitude.values.flatten(), 
                        gdap.G2pressure.values.flatten()),
                        axis=1)
Measurements = np.stack((gdap.G2salinity.values.flatten(), 
                         gdap.G2temperature.values.flatten(), 
                         gdap.G2oxygen.values.flatten()),
                         axis=1)
MeasIDVec = [1, 7, 6]
                                
results = carbon_utils.LIPHR_matlab(LIPHR_path,
                                    Coordinates.tolist(),
                                    Measurements.tolist(),
                                    MeasIDVec, 
                                    OAAdjustTF = False)                                  

gdap['pH_in_situ_total'] = results
gdap.pH_in_situ_total[np.isnan(gdap.G2phts25p0)] = np.nan
# gdap pH 25C 
gdap['pH_25C_TOTAL_ADJUSTED'] = carbon_utils.co2sys_pH25C(2300.,gdap.pH_in_situ_total,gdap.G2temperature,
                                                         gdap.G2salinity,gdap.G2pressure)
#set pH to nan where there was no original pH data from GLODAP
gdap.pH_25C_TOTAL_ADJUSTED[np.isnan(gdap.G2phts25p0)]=np.nan

#rename GLODAP comparison variables to match argo
gdap = gdap.rename(columns={'G2longitude':'LONGITUDE', 'G2latitude':'LATITUDE', 'G2pressure':'PRES_ADJUSTED',
                            'G2temperature':'TEMP_ADJUSTED','G2salinity':'PSAL_ADJUSTED', 
                            'G2oxygen':'DOXY_ADJUSTED','G2nitrate':'NITRATE_ADJUSTED', 'G2tco2':'DIC', 
                            'G2talk':'TALK_LIAR', 'G2MLD':'MLD','G2o2sat':'o2sat', 'G2PTMP':'PTMP', 
                            'pH_in_situ_total':'PH_IN_SITU_TOTAL_ADJUSTED','sigma0_calculated':'PDENS'})

gdap['obs_index']=gdap.reset_index().index
# -


# ## 2. Apply float bias corrections 

# +
start_time = time.perf_counter()

# 0: overwrites and runs all floats in the argo_path directory 
# 1: reads in and adds to argo_interp_temp.nc rather than overwriting and running all floats
# 2: runs specific floats listed below
append_data = 0

# when making major changes, list version number here
ver_n = '6a' 
# v2 - moving interpolated spice and density calculation to post-PSAL and TEMP interpolation
# v3 - fixed PH_25C calculation for float data, fixed in situ pH comparison (I think)
# v4 - added back in SI and NO3 to DIC calculation - makes a difference apparently (also changes which points have valid data)
# v5 - trying to do near-surface comparisons as well 
# v6 - working on full depth comparison that I will then separate by depth 

### 
if 'argo_interp' in locals():
    argo_interp.close()
    
argolist = []

for file in os.listdir(argo_path):
    if file.endswith('Sprof.nc'):
        argolist.append(file)
argolist.sort()

LIAR_path = liar_dir

#float QC data fields
qc_data_fields = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 
                  'PRES_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED']

bgc_data_fields = ['DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED']

# variables to do interpolation on:
interpolation_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED',
            'pH_25C_TOTAL_ADJUSTED', 'PRES_ADJUSTED', 'DIC','TALK_LIAR']
#variables to save to derived file
derived_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED',
            'pH_25C_TOTAL_ADJUSTED', 'PDENS', 'spice', 'PRES_ADJUSTED', 'DIC','TALK_LIAR']
# #variables to do crossover calculation
# crossover_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED',
#             'pH_25C_TOTAL_ADJUSTED', 'PDENS', 'spice', 'PRES_ADJUSTED', 'DIC','TALK_LIAR']

#if append data is set to 1, reads in argo_interp_temp.nc which contains prior argo_interp array, 
#compare wmo numbers between argolist and the wmo numbers in argo_interp, and continues on processing 
#float files. otherwise, start from the beginning
if append_data==1 and os.path.exists(data_dir+'argo_interp_temp.nc'):
    #load previously saved argo_interp
    argo_interp = xr.load_dataset(data_dir+'argo_interp_temp.nc')

    #extract wmo #s as integers from argolist
    s = [int(s) for s in re.findall("[0-9]+", str(argolist))]
    
    #find indices of argolist where argo_interp has no existing matching wmos
    indices = [s.index(i) for i in s if i not in argo_interp.wmo]
    argolist_run = [argolist[i] for i in indices]
elif append_data==0:
    argolist_run=argolist
    if 'argo_interp' in locals():
        del argo_interp # deletes argo_interp in case this code is being run multiple times. 
else:
#     argolist_run = ['1902385_Sprof.nc']
# argolist_run = ['7900566_Sprof.nc',
#                      '7900560_Sprof.nc',
#                      '6904115_Sprof.nc', 
#                      '6903557_Sprof.nc']
# #     argolist_run = ['1900722_Sprof.nc']
#     argolist_run = ['1901154_Sprof.nc',
#                      '1902332_Sprof.nc',
#                      '2900961_Sprof.nc', 
#                      '4900871_Sprof.nc', 
#                      '4901135_Sprof.nc', 
#                      '5901447_Sprof.nc', 
#                      '5901451_Sprof.nc',
#                      '4901135_Sprof.nc',
#                      '5906312_Sprof.nc']
    argolist_run = ['5906547_Sprof.nc',
                     '5906548_Sprof.nc',
                     '5906549_Sprof.nc', 
                     '5906550_Sprof.nc', 
                     '5906551_Sprof.nc', 
                     '5906552_Sprof.nc', 
                     '5906553_Sprof.nc',
                     '5906562_Sprof.nc',
                     '5906554_Sprof.nc',
                     '5906561_Sprof.nc', 
                     '5906556_Sprof.nc',
                     '5906558_Sprof.nc',
                     '5906559_Sprof.nc',
                     '5906557_Sprof.nc']
    if 'argo_interp' in locals():
        del argo_interp # deletes argo_interp in case this code is being run multiple times. 
  
       # argolist_run = ['4901135_Sprof.nc']

        
#####
#iterate through each float file 
#calculates derived carbonate parameters (currently TALK, DIC), bias corrected pH and stores all in argo_n
#interpolates all variables in "interpolation_list" to 1 m resolution and stores in argo_interp_n
#saves out individual float netcdf files with variables to be adjusted/used for crossovers
#appends 1m interpolated dataset in argo_interp for comparison to glodap (not saved currently)

argo_interp_filename = 'argo_interp_temp_' + str(dist) + 'km_' \
    + str(p_compare_min) + '_to_' + str(p_compare_max) + '_' + str(delta_press) + 'm_' + \
    str(delta_dens) + 'dens_' + str(delta_spice) + 'spice' + '_' + ver_n + '.nc'

wmo_list= list()
for n in range(len(argolist_run)):
    print(str(n)+' Processing float file '+ argolist_run[n])
    argo_n = xr.load_dataset(argo_path+argolist_run[n])
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))

    wmo_n = argo_n.PLATFORM_NUMBER.values.astype(int)[0]
    wmo_list.append(wmo_n)
    
    nprof_n = argo_n.dims['N_PROF']
    
    #   set bad data and possibly bad data to NaN 
    for q in qc_data_fields:      
        if q in argo_n.keys():
            qc_val = argo_n[q+'_QC'].values.astype('float')
            
            # for some reason the .where statement was not filtering out bad values. 
            #This code is now changing QC values of 3 or 4 to nans, not sure if it is the best approach
            #argo_n[q].where(np.logical_and(qc_val<3.,qc_val>4.))
            argo_n[q].values[np.logical_or(qc_val==4,qc_val==3)]=np.nan
            
            #check for any Inf values not included in QC flag and set to NaN
            argo_n[q].values[np.isinf(argo_n[q]).values] = np.nan
      
    # check for interpolated profile positions (under ice) and set all BGC data to nan
    qc_val = argo_n['POSITION_QC'].values.astype('float')
    for b in bgc_data_fields:
        if b in argo_n.keys() and np.any(qc_val==8):
            naninds = np.argwhere(qc_val==8)[:,0]
            argo_n[b][naninds,:] = np.nan
    
    # we are currently processing floats that have no valid biogeochemical data. 
    #Should check to see if data in key 
    #original bgc parameters (O2, NO3, pH) is valid and skip the rest if not
    bgc_valid = 0
    for b in bgc_data_fields:
        if b in argo_n.keys() and np.any(~np.isnan(argo_n[b])):
            bgc_valid = bgc_valid+1
    if bgc_valid >=1:
        print('float has valid BGC data')
    else:
        print('float has no valid BGC data')
        continue
    
    argo_n['PDENS'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
    argo_n.PDENS[:] = np.nan
    argo_n['spice'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
    argo_n.spice[:] = np.nan
    
    #initialise interpolated dataset for float
    nan_interp = np.empty((nprof_n,p_interp.shape[0]))
    nan_interp[:] = np.nan
    argo_interp_n = xr.Dataset()
    argo_interp_n['wmo'] = (['N_PROF'],np.repeat(wmo_n,nprof_n))
    argo_interp_n['profile'] = (['N_PROF'],argo_n.CYCLE_NUMBER.data) # added .data 
    argo_interp_n['juld'] = (['N_PROF'],argo_n.JULD_LOCATION.data)
    #add lat -lons to Dataset
    argo_interp_n['LATITUDE']  = (['N_PROF'],argo_n.LATITUDE.data)
    argo_interp_n['LONGITUDE']  = (['N_PROF'],argo_n.LONGITUDE.data)
    argo_interp_n['num_var'] = (['N_PROF'],np.zeros((nprof_n))) # changed from np.empty to np.zeros to avoid filling array with random large numbers
    for v in derived_list: # all the variables that will be saved out in the derived and interpolated files
        argo_interp_n[v] = (['N_PROF','N_LEVELS'],np.copy(nan_interp))

    #check first if PH_IN_SITU_TOTAL_ADJUSTED exists
    if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
        
        print('Calculating TALK, DIC and pH 25C correction for float '+str(wmo_n))
        
        #initialise pH 25c and DIC variables - could do this only if float has pH
        argo_n['TALK_LIAR'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.TALK_LIAR[:] = np.nan
        argo_n['pH_25C_TOTAL_ADJUSTED'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.pH_25C_TOTAL_ADJUSTED[:] = np.nan
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.DIC[:] = np.nan
#         argo_n['pH_insitu_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
#         argo_n.pH_insitu_corr[:] = np.nan
#         argo_n['pH_25C_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
#         argo_n.pH_25C_corr[:] = np.nan
#         argo_n['bias_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof
#         argo_n.bias_corr[:] = np.nan

    
        ##### Calc float TALK       
        #repeat lats, lons to match pressure shape
        lons_rep = np.tile(argo_n.LONGITUDE.values,(argo_n.PRES_ADJUSTED.shape[1],1)).T
        lats_rep = np.tile(argo_n.LATITUDE.values,(argo_n.PRES_ADJUSTED.shape[1],1)).T

        #set Si and PO4 inputs
        #if nitrate, then use redfield for Si and PO4?, otherwise set to 0    
        if 'NITRATE_ADJUSTED' in argo_n.keys():
            SI = argo_n.NITRATE_ADJUSTED*2.5
            SI.where(~np.isnan(SI), 0)
            PO4 = argo_n.NITRATE_ADJUSTED/16
            PO4.where(~np.isnan(PO4),0)
            Coordinates = np.stack((lons_rep.flatten(), 
                            lats_rep.flatten(), 
                            argo_n.PRES_ADJUSTED.values.flatten()),
                            axis=1)
            Measurements = np.stack((argo_n.PSAL_ADJUSTED.values.flatten(), 
                             argo_n.TEMP_ADJUSTED.values.flatten(), 
                             argo_n.NITRATE_ADJUSTED.values.flatten(), 
                             argo_n.DOXY_ADJUSTED.values.flatten()),
                             axis=1)
            MeasIDVec = [1, 7, 3, 6]
   
        else:
            SI = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
            PO4 = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
            Coordinates = np.stack((lons_rep.flatten(), 
                            lats_rep.flatten(), 
                            argo_n.PRES_ADJUSTED.values.flatten()),
                            axis=1)
            Measurements = np.stack((argo_n.PSAL_ADJUSTED.values.flatten(), 
                             argo_n.TEMP_ADJUSTED.values.flatten(),
                             argo_n.DOXY_ADJUSTED.values.flatten()),
                             axis=1)
            MeasIDVec = [1, 7, 6]                            


        results = carbon_utils.LIAR_matlab(LIAR_path,
                                               Coordinates.tolist(),
                                               Measurements.tolist(),
                                               MeasIDVec,
                                               VerboseTF=False)                                  

        argo_n['TALK_LIAR'] = (['N_PROF','N_LEVELS'],
                                np.reshape(np.asarray(results),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
  
    
        # Keep DIC bc I might want it for crossover comparison
        ##### Calculate float pH at 25C, DIC and apply bias corr
        results = pyco2.sys(
                par1=argo_n.TALK_LIAR, 
                par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
                par1_type=1,
                par2_type=3,
                temperature=argo_n.TEMP_ADJUSTED, 
                pressure=argo_n.PRES_ADJUSTED, 
                salinity=argo_n.PSAL_ADJUSTED, 
                temperature_out=25.,#*np.ones(argo_n.PRES_ADJUSTED.shape), #fixed 25C temperature
                pressure_out=argo_n.PRES_ADJUSTED,
                total_silicate=SI,
                total_phosphate=PO4,
                opt_pH_scale = 1, #total
                opt_k_carbonic=10, #Lueker et al. 2000
                opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
                opt_total_borate=2, # Lee et al. 2010
                opt_k_fluoride=2, # Perez and Fraga 1987
                opt_buffers_mode=1,
        )
        
#         argo_n['pH_25C_TOTAL_ADJUSTED'] = (['N_PROF','N_LEVELS'],results['pH_total_out'])
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],results['dic'])  
        
        argo_n['pH_25C_TOTAL_ADJUSTED'] = (['N_PROF','N_LEVELS'],carbon_utils.co2sys_pH25C(argo_n.TALK_LIAR,
                                                    argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
                                                    argo_n.TEMP_ADJUSTED,
                                                    argo_n.PSAL_ADJUSTED,
                                                    argo_n.PRES_ADJUSTED))
# #         gdap['pH_25C_TOTAL_ADJUSTED'] = carbon_utils.co2sys_pH25C(2300.,gdap.pH_in_situ_total,gdap.G2temperature,
#                                                          gdap.G2salinity,gdap.G2pressure)
      
#         for p in range(nprof_n):
            
#             # skip a profile if pH is above 10.  There seem to be pH's above 10 that causing 
#             if any(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]>10) or all(np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:])):
#                 'print pH out of range'
#                 continue
    
#             #apply pH bias correction   
#             #find if there are valid values between fixed p levels 
#             pH_p_min = 1480
#             pH_p_max = 1520
    
#             #if there are valid pressure levels between 1480-1520 db, 
#             #calc bias correction only in this depth band, if not, calc correction between 970 and 1520
#             if any((argo_n.PRES_ADJUSTED[p,:]>1480) & (argo_n.PRES_ADJUSTED[p,:]<1520)):
                
#                 inds = (argo_n.PRES_ADJUSTED[p,:]>1480) & (argo_n.PRES_ADJUSTED[p,:]<1520)
#                 correction = -0.034529*argo_n.pH_25C_TOTAL_ADJUSTED[p,inds]+0.26709
                              
#             else:
#                 inds = (argo_n.PRES_ADJUSTED[p,:]>970) & (argo_n.PRES_ADJUSTED[p,:]<1520)
#                 correction = -0.034529*argo_n.pH_25C_TOTAL_ADJUSTED[p,inds]+0.26709

#             if len(correction):
#                 argo_n.bias_corr[p] = np.nanmean(correction)
#                 argo_n.pH_insitu_corr[p,:] = argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]+argo_n.bias_corr[p]
#                 argo_n.pH_25C_corr[p,:] = argo_n.pH_25C_TOTAL_ADJUSTED[p,:]+argo_n.bias_corr[p]
    
#         #call CO2sys again to get pCO2 with corrected PH- do we need to include this here?
#         results = pyco2.sys(
#                 par1=argo_n.TALK_LIAR, 
#                 par2=argo_n.pH_insitu_corr,
#                 par1_type=1,
#                 par2_type=3,
#                 temperature=argo_n.TEMP_ADJUSTED, 
#                 pressure=argo_n.PRES_ADJUSTED, 
#                 salinity=argo_n.PSAL_ADJUSTED, 
#                 temperature_out=25.* np.ones(argo_n.PRES_ADJUSTED.shape), #fixed 25C temperature
#                 pressure_out=argo_n.PRES_ADJUSTED,
#                 total_silicate=SI,
#                 total_phosphate=PO4,
#                 opt_pH_scale = 1, #total
#                 opt_k_carbonic=10, #Lueker et al. 2000
#                 opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
#                 opt_total_borate=2, # Lee et al. 2010
#                 opt_k_fluoride=2, # Perez and Fraga 1987
#                 opt_buffers_mode=1,
#                 )
 
#         argo_n['pCO2_pH_corr'] = (['N_PROF','N_LEVELS'],results['pCO2'])  
    
    ##### now calc potential density, save, and interpolate data for comparison
    for p in range(nprof_n):
        #pressure for profile
        p_prof = argo_n.PRES_ADJUSTED[p,:]
        
        # For interpolated data, shouldn't calculate pdens and spice and then interpolate - 
        # should interpolate psal and temp and then calculate spice and pdens
        # Do both so that you are able to have PDENS and spice in the derived files too (do I need them?)
        argo_n['PDENS'][p,:] = carbon_utils.sigma0(argo_n.PSAL_ADJUSTED[p,:].values,
                                                   argo_n.TEMP_ADJUSTED[p,:].values,
                                                  argo_n.LONGITUDE[p].values,
                                                  argo_n.LATITUDE[p].values,
                                                  argo_n.PRES_ADJUSTED[p,:].values)
        argo_n['spice'][p,:] = carbon_utils.spiciness0(argo_n.PSAL_ADJUSTED[p,:].values,
                                                   argo_n.TEMP_ADJUSTED[p,:].values,
                                                  argo_n.LONGITUDE[p].values,
                                                  argo_n.LATITUDE[p].values,
                                                  argo_n.PRES_ADJUSTED[p,:].values)

        #for each profile get pressure values > p_interp_min db
        p100 = p_prof[p_prof>p_interp_min].values
            
        #if only 1 value of pressure or if there is not valid profile data down to p-max, continue loop
        if (len(p100) <= 1) or (np.nanmax(p100)<p_interp_min):
            continue
        
        #find which crossover variables exist in main float file
        var_list_n = []
        for vname in interpolation_list:
            if (vname in argo_n.keys()) and (np.any(~np.isnan(argo_n[vname]))):
                var_list_n.append(vname)
                
        argo_interp_n['num_var'][p] = len(var_list_n) 
        
        for var in var_list_n:
            var100 = argo_n[var][p,p_prof>p_interp_min]

            #if there are non-unique pressure values, 
            #then grab only unique pressure values and matching data points
            if len(p100)>len(np.unique(p100)):
                p100u,unique_inds = np.unique(p100, return_index=True)
                var100u = var100[unique_inds]
            else:
                p100u = p100
                var100u = var100
                
            #interpolate 1d profile data onto p_interp levels 
            # use valid var data from p_interp_min to p_interp_max OR maximum valid pressure 
            #(greater than minimum comparison pressure)

            if len(p100u[~np.isnan(var100u.values)])>1 and \
                (np.nanmax(p100u[~np.isnan(var100u.values)])>p_interp_min) and \
                (np.nanmin(p100u[~np.isnan(var100u.values)])<p_interp_max):
                
                #interpolation function
                f = interpolate.interp1d(p100u[~np.isnan(var100u.values)],var100u[~np.isnan(var100u.values)])
                
                #check if non-NaN data does not extend down to p_interp_max
                if np.logical_and((p100u[~np.isnan(var100u.values)][-1]<p_interp_max),
                                  (p100u[~np.isnan(var100u.values)][0]>p_interp_min)):
                    pmin_ind = np.argwhere(p_interp>p100u[~np.isnan(var100u.values)][0])[0][0]
                    pmax_ind = np.argwhere(p_interp>p100u[~np.isnan(var100u.values)][-1])[0][0]
                    #if  p100u[~np.isnan(var100u)][0]>p_interp_min:                   
                    var_interp_p = f(p_interp[pmin_ind:pmax_ind])
                    #assign interpolated variables to array 
                    argo_interp_n[var][p,pmin_ind:pmax_ind] = var_interp_p
                    
                elif p100u[~np.isnan(var100u.values)][-1]<p_interp_max:
                    pmax_ind = np.argwhere(p_interp>p100u[~np.isnan(var100u.values)][-1])[0][0]
                    var_interp_p = f(p_interp[:pmax_ind])
                    #assign interpolated variables to array 
                    argo_interp_n[var][p,:pmax_ind] = var_interp_p
                        
                elif p100u[~np.isnan(var100u.values)][0]>p_interp_min:
                    pmin_ind = np.argwhere(p_interp>p100u[~np.isnan(var100u.values)][0])[0][0]
                    var_interp_p = f(p_interp[pmin_ind:])
                    #assign interpolated variables to array 
                    argo_interp_n[var][p,pmin_ind:] = var_interp_p
                    
                else:
                    var_interp_p = f(p_interp)
                    #assign interpolated variables to array 
                    argo_interp_n[var][p,:] = var_interp_p
            
#             else: 
                # print('profile data not deep enough to interpolate ' + str(p) + ' ' +  var)
                #                       str(np.nanmax(p100u[~np.isnan(var100u.values)])))
                # print('values greater than min ' + str(var100u[p100u>p_interp_min].values))

    # loop through profiles again to calculate PDENS and spice for interpolated dataset
    for p in range(nprof_n):
        #pressure for profile
        p_prof = argo_interp_n.PRES_ADJUSTED[p,:]

        # For interpolated data, shouldn't calculate pdens and spice and then interpolate - 
        # should interpolate psal and temp and then calculate spice and pdens
        # Do both so that you are able to have PDENS and spice in the derived files too (do I need them?)
        argo_interp_n['PDENS'][p,:] = carbon_utils.sigma0(argo_interp_n.PSAL_ADJUSTED[p,:].values,
                                                   argo_interp_n.TEMP_ADJUSTED[p,:].values,
                                                  argo_interp_n.LONGITUDE[p].values,
                                                  argo_interp_n.LATITUDE[p].values,
                                                  argo_interp_n.PRES_ADJUSTED[p,:].values)
        argo_interp_n['spice'][p,:] = carbon_utils.spiciness0(argo_interp_n.PSAL_ADJUSTED[p,:].values,
                                                   argo_interp_n.TEMP_ADJUSTED[p,:].values,
                                                  argo_interp_n.LONGITUDE[p].values,
                                                      argo_interp_n.LATITUDE[p].values,
                                                      argo_interp_n.PRES_ADJUSTED[p,:].values)
        
    #create new dataset with relevant crossover variables only
    argo_n_derived = xr.Dataset()
    argo_n_derived['wmo'] = wmo_n
    argo_n_derived['CYCLE_NUMBER'] = (['N_PROF'],argo_n.CYCLE_NUMBER.values)
    argo_n_derived['LONGITUDE'] = (['N_PROF'],argo_n.LONGITUDE.values)
    argo_n_derived['LATITUDE'] = (['N_PROF'],argo_n.LATITUDE.values)
    argo_n_derived['JULD_LOCATION'] = (['N_PROF'],argo_n.JULD_LOCATION.values)
    for var in derived_list:
        if var in argo_n.keys():
            argo_n_derived[var] = (['N_PROF','N_LEVELS'],argo_n[var].values)
    argo_n_derived.to_netcdf(argo_path_derived+str(wmo_n)+'_derived.nc')
    
    #Also save interpolated dataset as one mega Dataset for doing crossovers
    if n == 0 or 'argo_interp' not in locals(): # modified to deal w/ situation where n==0 skipped defining argo_interp
        argo_interp = argo_interp_n
    else:
        argo_interp = xr.concat([argo_interp,argo_interp_n],'N_PROF')
                    
    #save out argo_interp periodically:
    if n/20==round(n/20):
        argo_interp.to_netcdf(data_dir+argo_interp_filename)

argo_interp.to_netcdf(data_dir+argo_interp_filename) # need to save one final time so that final floats are not missed
argo_interp.close() # close dataset to avoid keeping in memory by accident 

## 3. Compare float - GLODAP crossovers

if 'argo_wmo' in locals():
       del argo_wmo # deletes argo_interp in case this code is being run multiple times. 

glodap_offsets_filename = 'glodap_offsets_' + str(dist) + 'km_' \
    + str(p_compare_min) + '_to_' + str(p_compare_max) + '_' + str(delta_press) + 'm_' + \
    str(delta_dens) + 'dens_' + str(delta_spice) + 'spice' + '_' + ver_n + '.nc'

        
#toggle to plot offsets profile by profile
plot_profile = 0
    
#restrict glodap data to comparison pressure range
gdap_p = gdap[(gdap.PRES_ADJUSTED.values>p_compare_min) & (gdap.PRES_ADJUSTED.values<p_compare_max)]

#load saved argo_interp data if needed
argo_interp = xr.open_dataset(data_dir+argo_interp_filename)
#group by float wmo
argo_wmo = argo_interp.groupby('wmo')

#initiate offset list
#number of additional rows for saving metadata items
num_meta_items = 11
gdap_offsets =  [[] for _ in range(3*len(var_list_plot)+num_meta_items)]

#iterate over each float & profile
float_count = 0

for wmo, group in argo_wmo:
        
    #number of profiles
    nprof = group.LATITUDE.shape[0]
    
     #check sufficient non-NaN data
    if np.sum(group.num_var>3)==0: #changed to look at all profiles and to only select floats 
        # with more than 3 variables present (i.e. more than T, S, Press)
        print('No non-NAN bgc adjusted data for: '+str(wmo))
        continue
    
    float_count = float_count+1

    # if float_count==20:
    #    break
        
    print('Float ' + str(float_count) + ' ' + str(wmo))
    
    #  Find lat-lon limits within dist of float location
    lat_tol = dist/ 111.6 
    lon_tol = dist/ (111.6*np.cos(np.deg2rad(np.nanmean(group.LATITUDE))))  
    #set lat/lon crossover limits
    lat_min = group.LATITUDE.values-lat_tol
    lat_max = group.LATITUDE.values+lat_tol
    lon_min = group.LONGITUDE.values-lon_tol
    lon_max = group.LONGITUDE.values+lon_tol

    #find all data in lat-lon limits
    for n in range(nprof): #range(4,5):# range(3, 4): #   #
        #print(group.profile[n].values)

        #index of all gdap profiles within distance range
        if lon_min[n] < 0 and lon_max[n]>0:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n]+360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]))
            
        elif lon_min[n] < 0 and lon_max[n]<0:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_and(gdap_p.LONGITUDE.values>lon_min[n]+360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]+360))  
            
        elif lon_max[n] > 360 and lon_min[n] < 360:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n],
                                                  gdap_p.LONGITUDE.values<lon_max[n]-360))

        elif lon_max[n] > 360 and lon_min[n] > 360:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_and(gdap_p.LONGITUDE.values>lon_min[n]-360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]-360))

        else:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_and(gdap_p.LONGITUDE.values>lon_min[n],
                                                  gdap_p.LONGITUDE.values<lon_max[n]))

        match = np.squeeze(match)
        match_inds = np.argwhere(match)
        

        #get matched glodap data subset
        if len(match_inds)==0:
            continue

        
        #get glodap data points that match
        gdap_match = gdap_p[match]

        
        #find index of interp profile density that most closely matches each glodap match point density
        dens_inds = np.empty(match_inds.shape[0],dtype=int)
        o2_offsets = np.empty(match_inds.shape[0])
        o2_offsets[:] = np.nan
        dens_offsets = np.empty(match_inds.shape[0])
        o2_offsets[:] = np.nan
        

        for m in range(len(match_inds)):
        
            #if no interpolated density (not deep enough) then move on to next profile
            # probably should move the float profile check earlier in this cell, no need to run so much extra code 
            if np.all(np.isnan(group.PDENS.values[n,:])) or np.isnan(gdap_match.PDENS.values[m]):
#                 print('no interpolated density in depth range')
                dens_inds[m]=-1
                continue
                
            #calculate density difference and nearest density index
            dens_diff = np.absolute(group.PDENS.values[n,:]-gdap_match.PDENS.values[m])
            dens_ind = np.nanargmin(dens_diff)

            dens_inds[m] = dens_ind #save indices for plotting 

            
            #don't use matches if the density difference to the closest matching density value exceeds delta_dens
            if dens_diff[dens_ind] > delta_dens:
                continue
                
            #don't use matches at top/bottom of interpolated vector
            if (dens_ind == 0) or (dens_ind == len(p_interp) -1):
                continue
                
            ##optional: add spice criteria
            spice_diff = np.absolute(group.spice.values[n,dens_ind]-gdap_match.spice.values[m])
            
            #don't use matches if spice difference exceeds delta_spice
            if spice_diff > delta_spice:
                continue
                
            ## add depth criteria
            press_diff = np.absolute(group.PRES_ADJUSTED.values[n,dens_ind]-gdap_match.PRES_ADJUSTED.values[m])
            if press_diff> delta_press:
                continue
            #calc offset at the interp profile index for each match for each var to calculate crossover and or plot 
            #check var is present in both float and glodap
            for idx, var in enumerate(var_list_plot):
                if var in group.keys():
                    gdap_offset_ind = gdap_match[var].index[m]
                    offset = group[var][n,dens_ind] - gdap_match[var][gdap_offset_ind]
                    
                    #save offsets for optional plotting
                    if var=='DOXY_ADJUSTED':
                        o2_offsets[m] = offset.values
                    elif var=='PDENS':
                        dens_offsets[m] = offset.values
                        
                    #append to full offset list
                    gdap_offsets[idx*3].append(offset.values)
                    
                    #also save absolute glodap value at crossover
                    gdap_offsets[idx*3+1].append(gdap_match[var][gdap_offset_ind])
                    #save absolute float value at crossover
                    gdap_offsets[idx*3+2].append(group[var][n,dens_ind])
#                     print(var)
                #append nan if variable is not there so lists all remain same length?
                else:
                    gdap_offsets[idx*3].append(np.nan) 
                    gdap_offsets[idx*3+1].append(np.nan) 
                    gdap_offsets[idx*3+2].append(np.nan) 
            
            #append metadata to offset list
            gdap_offsets[len(var_list_plot)*3].append(wmo)
            gdap_offsets[len(var_list_plot)*3+1].append(group.profile[n].values)
            gdap_offsets[len(var_list_plot)*3+2].append(group.juld[n].values)
            gdap_offsets[len(var_list_plot)*3+3].append(gdap_match.datetime.values[m])
            gdap_offsets[len(var_list_plot)*3+4].append(group.LONGITUDE[n].values)
            gdap_offsets[len(var_list_plot)*3+5].append(gdap_match.LONGITUDE.values[m])
            gdap_offsets[len(var_list_plot)*3+6].append(group.LATITUDE[n].values)
            gdap_offsets[len(var_list_plot)*3+7].append(gdap_match.LATITUDE.values[m])
            gdap_offsets[len(var_list_plot)*3+8].append(gdap_match.G2cruise.values[m])
            gdap_offsets[len(var_list_plot)*3+9].append(gdap_match.G2station.values[m])
            gdap_offsets[len(var_list_plot)*3+10].append(gdap_match.obs_index.values[m])

            #can add additional float metadata variable to list here
     

        ################################
        #optional: plot individual profiles and offsets to check density match up
        if plot_profile == 1:
            # if there are no non-nan values of o2_offsets for this profile, move on
            if np.all(np.isnan(o2_offsets)):
                continue
            
            #if no interpolated density (not deep enough) then move on to next profile
            if np.all(np.isnan(group.PDENS.values[n,:])) or np.isnan(gdap_match.PDENS.values[m]):
                continue
                
            # create a folder for profile plots if one does not exist:
            if not os.path.isdir(offset_dir + 'individual_floats'):
                os.mkdir(offset_dir + 'individual_floats')
                
            if not os.path.isdir(offset_dir + 'individual_floats/' + str(wmo)):
                os.mkdir(offset_dir + 'individual_floats/' + str(wmo))
    
            print('Plotting float '+str(group.wmo[0].values) + ' profile' + str(group.profile[n].values))
            fig = plt.figure(figsize=(10,16))
            
            #DOXY_ADJUSTED
            #dens_inds = dens_inds[dens_inds>0] #only keep valid indices
            
            max_o2 = 0
            min_o2 = 300
            max_dens = 0
            min_dens = 300
            min_pres = 1600
            max_pres = 1700
            
            axn = plt.subplot(3,3,1)
            #plot float profile and matchups
            axn.scatter(group.DOXY_ADJUSTED[n,:],group.PRES_ADJUSTED[n,:],s=4,c='orangered',marker='o')
            # commenting this out, because dens_inds is specific to an individual glodap crossover   
            # if group.DOXY_ADJUSTED[n, dens_inds].size>0:
            #     axn.scatter(group.DOXY_ADJUSTED[n,dens_inds],group.PRES_ADJUSTED[n,dens_inds],s=50,c='k',marker='s')
            #plot glodap profiles and matchups- (split glodap into profiles by cruise and station?)

            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match.DOXY_ADJUSTED.index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    #print(m)
                    axn.scatter(group.DOXY_ADJUSTED[n,dens_inds[m]],group.PRES_ADJUSTED[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                    #plot all glodap data that match these cruise and stations
                    axn.plot(gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'DOXY_ADJUSTED'].values,
                         gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PRES_ADJUSTED'].values,
                        'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.DOXY_ADJUSTED.values[m],
                            gdap_match.PRES_ADJUSTED.values[m],s=50,c='b',marker='D')
                    
                                    
                    if np.nanmin(group.DOXY_ADJUSTED[n,dens_inds[m]])<min_o2:
                        min_o2 = np.nanmin(group.DOXY_ADJUSTED[n,dens_inds[m]])
                    
                    if np.nanmin(gdap_match.DOXY_ADJUSTED.values[m])<min_o2:
                        min_o2 = np.nanmin(gdap_match.DOXY_ADJUSTED.values[m])

                    if np.nanmax(gdap_match.DOXY_ADJUSTED.values[m])>max_o2:
                        max_o2 = np.nanmax(gdap_match.DOXY_ADJUSTED.values[m])
                    
                    if np.nanmax(group.DOXY_ADJUSTED[n,dens_inds[m]])>max_o2:
                        max_o2 = np.nanmax(group.DOXY_ADJUSTED[n,dens_inds[m]])
                                        
                    if np.nanmin(group.PDENS[n,dens_inds[m]])<min_dens:
                        min_dens = np.nanmin(group.PDENS[n,dens_inds[m]])
                    
                    if np.nanmax(group.PDENS[n,dens_inds[m]])>max_dens:
                        max_dens = np.nanmax(group.PDENS[n,dens_inds[m]])
                        
                    if np.nanmin(group.PRES_ADJUSTED[n,dens_inds[m]])<min_pres:
                        min_pres = np.nanmin(group.PRES_ADJUSTED[n,dens_inds[m]])
                    
                    if np.nanmax(group.PRES_ADJUSTED[n,dens_inds[m]])>max_pres:
                        max_pres = np.nanmax(group.PRES_ADJUSTED[n,dens_inds[m]])

            plt.gca().invert_yaxis()
            plt.ylim([2000,1200])
            plt.ylabel('pressure')
            plt.xlabel('DOXY_ADJUSTED')
            plt.title(str(group.wmo[0].values))
            plt.xlim([min_o2 - 20, 
                      max_o2 + 20])            
             
             #PDENS
            axn = plt.subplot(3,3,2)
            #plot float profile and matchups
            axn.scatter(group.PDENS[n,:],group.PRES_ADJUSTED[n,:],s=4,c='orangered',marker='o')
            
            #plot glodap profiles and matchups- (split glodap into profiles by cruise and station?)
            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match[var].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
    
                    axn.scatter(group.PDENS[n,dens_inds[m]],group.PRES_ADJUSTED[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                                #plot all glodap data that match these cruise and stations
                    axn.plot(gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PDENS'].values,
                         gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PRES_ADJUSTED'].values,
                        'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.PDENS.values[m],
                            gdap_match.PRES_ADJUSTED.values[m],s=50,c='b',marker='D')
            plt.ylim([min_pres- 100, max_pres+100])

            plt.gca().invert_yaxis()
#             plt.ylim([2000,1200])
#             plt.xlim([27.5,28.0])
            plt.ylabel('pressure')
            plt.xlabel('PDENS')

            plt.xlim([min_dens-.01, max_dens+.01])

            
           #DOXY vs PDENS
            axn = plt.subplot(3,3,4)
            #plot float profile and matchups
            axn.scatter(group.DOXY_ADJUSTED[n,:],group.PDENS[n,:],s=4,c='orangered',marker='o')

            
            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match[var].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    axn.scatter(group.DOXY_ADJUSTED[n,dens_inds[m]],group.PDENS[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                    #plot all glodap data that match these cruise and stations
                    axn.plot(gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'DOXY_ADJUSTED'].values,
                             gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PDENS'].values,
                            'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.DOXY_ADJUSTED.values[m],
                                gdap_match.PDENS.values[m],s=50,c='b',marker='D')
    
            plt.gca().invert_yaxis()
            plt.ylabel('PDENS')
            plt.xlabel('DOXY_ADJUSTED')
            plt.xlim([min_o2 - 10, 
                      max_o2 + 10])
            plt.ylim([min_dens-.05, max_dens+.05])

            #print cruise names and dates
            axn = plt.subplot(3,3,3)
            axn.text(0.05,0.95,'G2cruise G2station G2date',fontsize=12)
            yn=0.9
            for m in range(len(match_inds)):
                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover

                    off_ind = gdap_match['G2cruise'].index[m]
                
                    #if m == 0:
                    axn.text(0.05,yn,str(gdap_match.G2cruise[off_ind])+' '
                             +str(gdap_match.G2station[off_ind])+' '
                             +str(gdap_match.datetime[off_ind])+' '
                             +str(round(o2_offsets[m],2)))
                    print()
                    yn=yn-0.05
                    #elif gdap_match.G2cruise[off_ind] != cc:
                    #    axn.text(0.05,yn,str(gdap_match.G2cruise[off_ind])+' '+str(gdap_match.G2station[off_ind])+' '+str(gdap_match.datetime[off_ind]))
                    #    yn=yn-0.05    
                    #cc = gdap_match.G2cruise[off_ind] #current cruise to compare with next 
            
            # map of offsets:
         #   ax = plt.subplot(3,3,5,subplot_kw={'projection': ccrs.PlateCarree()})
            ax = plt.subplot(3,3,5,projection=ccrs.PlateCarree())

            #ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines()
            ax.set_global()
            # find any longitude values over 180, then subtract
            temp_lon = gdap_match.LONGITUDE.values
            temp_lon[temp_lon>180] = temp_lon[temp_lon>180]-360

            max_lon = np.nanmax(temp_lon[~np.isnan(o2_offsets)])
            min_lon = np.nanmin(temp_lon[~np.isnan(o2_offsets)])
            max_lat = np.nanmax(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])
            min_lat = np.nanmin(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])

            # only plot gdap locations where offsets exist
            #plt.plot(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)],marker='x')
            scat = plt.scatter(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)], c=o2_offsets[~np.isnan(o2_offsets)])
            plt.xlim([min_lon - .05, max_lon+.05])
            plt.ylim([min_lat - .05, max_lat+.05])
            plt.colorbar(scat, ax=ax)
            #plt.show()
            
            # plot concentration vs. offset
            axn = plt.subplot(3,3,6)
#             plt.plot(o2_offsets, group.DOXY_ADJUSTED[n,:], 'bx')
            
            for m in range(len(match_inds)):
#                 gdap_offset_ind = gdap_match[var].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    scat = plt.scatter(o2_offsets[m], group.DOXY_ADJUSTED[n,dens_inds[m]],s=50,
                                c=group.PRES_ADJUSTED[n,dens_inds[m]], marker='s')

#                     cruise = gdap_match.G2cruise[gdap_offset_ind]
#                     stat = gdap_match.G2station[gdap_offset_ind]
                    #plot all glodap data that match these cruise and stations
#                     axn.plot(gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'DOXY_ADJUSTED'].values,
#                              gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PDENS'].values,
#                             'm+-',linewidth=0.3)
#                     axn.scatter(gdap_match.DOXY_ADJUSTED.values[m],
#                                 gdap_match.PDENS.values[m],s=50,c='b',marker='D')
    
            plt.colorbar()

            
            #histogram of offsets
            axn = plt.subplot(3,3,7)
            o2_offsets[o2_offsets==np.inf] = np.nan
            if not np.all(np.isnan(o2_offsets)):
                axn.hist(o2_offsets[~np.isnan(o2_offsets)],color='b')
                plt.xlabel('DOXY_OFFSETS')
                
            axn = plt.subplot(3,3,8)
            dens_offsets[dens_offsets==np.inf] = np.nan
            if not np.all(np.isnan(dens_offsets)):
                axn.hist(dens_offsets[~np.isnan(o2_offsets)],color='b')
                plt.xlabel('DENS_OFFSETS')
                plt.xticks(rotation = 45)
            plt.tight_layout()

            plt.savefig(offset_dir+ 'individual_floats/' + str(wmo) + '/' + str(group.wmo[0].values) + ' profile' + str(int(group.profile[n].values)))
            plt.clf()


#convert GLODAP offset lists to xarray and save to netcdf
glodap_offsets = xr.Dataset(coords={'N_CROSSOVERS':(['N_CROSSOVERS'],np.arange(len(gdap_offsets[0])))})

glodap_offsets['p_compare_min'] = p_compare_min
glodap_offsets['p_compare_max'] = p_compare_max
glodap_offsets['delta_dens'] = delta_dens
glodap_offsets['delta_spice'] = delta_spice
glodap_offsets['delta_press'] = delta_press
glodap_offsets['dist'] = dist

for idx, var in enumerate(var_list_plot):
    glodap_offsets[var+'_offset'] = (['N_CROSSOVERS'], gdap_offsets[idx*3])
    glodap_offsets[var+'_glodap'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+1])
    glodap_offsets[var+'_float'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+2])
    
glodap_offsets['main_float_wmo'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3])
glodap_offsets['main_float_profile'] = (['N_CROSSOVERS'], gdap_offsets[len(var_list_plot)*3+1])
glodap_offsets['main_float_juld'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+2])
glodap_offsets['glodap_datetime'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+3])
glodap_offsets['main_float_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+4])
glodap_offsets['glodap_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+5])
glodap_offsets['main_float_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+6])
glodap_offsets['glodap_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+7])
glodap_offsets['glodap_cruise'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+8])
glodap_offsets['glodap_station'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+9])
glodap_offsets['glodap_obs_index'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+10])

glodap_offsets.to_netcdf(output_dir+glodap_offsets_filename)

print('Total number of glodap crossovers: ' + str(len(gdap_offsets[len(var_list_plot)*3])))
finish_time = time.perf_counter()
print("Program finished in {} seconds".format(finish_time - start_time))
print("---")

# +
print(np.nanmean(glodap_offsets['pH_25C_TOTAL_ADJUSTED_glodap']))
print(np.nanmean(glodap_offsets['pH_25C_TOTAL_ADJUSTED_float']))
print(np.nanmean(glodap_offsets['pH_25C_TOTAL_ADJUSTED_offset']))

print(np.nanmean(glodap_offsets['PH_IN_SITU_TOTAL_ADJUSTED_glodap']))
print(np.nanmean(glodap_offsets['PH_IN_SITU_TOTAL_ADJUSTED_float']))
print(np.nanmean(glodap_offsets['PH_IN_SITU_TOTAL_ADJUSTED_offset']))
print(var_list_plot)

# +
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
# find any longitude values over 180, then subtract
temp_lon = gdap_match.LONGITUDE.values
temp_lon[temp_lon>180] = temp_lon[temp_lon>180]-360

max_lon = np.nanmax(temp_lon[~np.isnan(o2_offsets)])
min_lon = np.nanmin(temp_lon[~np.isnan(o2_offsets)])
max_lat = np.nanmax(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])
min_lat = np.nanmin(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])

# only plot gdap locations where offsets exist
#plt.plot(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)],marker='x')
scat = plt.scatter(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)], c=o2_offsets[~np.isnan(o2_offsets)])
plt.xlim([min_lon - .05, max_lon+.05])
plt.ylim([min_lat - .05, max_lat+.05])
plt.colorbar(scat, ax=ax)
plt.show()

# +
for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match[var].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
    
                   # axn.scatter(group.PDENS[n,dens_inds[m]],group.PRES_ADJUSTED[n,dens_inds[m]],s=50,c='k',marker='s')

                   # cruise = gdap_match.G2cruise[gdap_offset_ind]
                  #  stat = gdap_match.G2station[gdap_offset_ind]
                                #plot all glodap data that match these cruise and stations
                  #  axn.plot(gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PDENS'].values,
                  #       gdap.loc[(gdap['G2cruise']==cruise) & (gdap['G2station']==stat),'PRES_ADJUSTED'].values,
                  #      'm+-',linewidth=0.3)
                   # axn.scatter(gdap_match.PDENS.values[m],
                  #          gdap_match.PRES_ADJUSTED.values[m],s=50,c='b',marker='D')
                    plt.plot(gdap_match.LONGITUDE.values[m], gdap_match.LATITUDE.values[m], linestyle='none', marker='.', markersize=10, color='m')
                    print(gdap_match.LONGITUDE.values[m])
                    print(gdap_match.LATITUDE.values[m])

                    #plt.plot(argo_all.LONGITUDE, argo_all.LATITUDE, linestyle='none', marker='.', markersize=1)

plt.show()
# -

match_inds

# +
#convert GLODAP offset lists to xarray and save to netcdf
glodap_offsets = xr.Dataset(coords={'N_CROSSOVERS':(['N_CROSSOVERS'],np.arange(len(gdap_offsets[0])))})
for idx, var in enumerate(var_list_plot):
    glodap_offsets[var+'_offset'] = (['N_CROSSOVERS'], gdap_offsets[idx*3])
    glodap_offsets[var+'_glodap'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+1])
    glodap_offsets[var+'_float'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+2])
    
glodap_offsets['main_float_wmo'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3])
glodap_offsets['main_float_profile'] = (['N_CROSSOVERS'], gdap_offsets[len(var_list_plot)*3+1])
glodap_offsets['main_float_juld'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+2])
glodap_offsets['glodap_datetime'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+3])
glodap_offsets['main_float_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+4])
glodap_offsets['glodap_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+5])
glodap_offsets['main_float_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+6])
glodap_offsets['glodap_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+7])

glodap_offsets.to_netcdf(output_dir+'glodap_offsets.nc')

print('Total number of glodap crossovers: ' + str(len(gdap_offsets[len(var_list_plot)*3])))
# -
# ## 4. Compare float - float crossovers

# +
#float- float crossover 
#variables to do crossover plot 
var_list_plot = ['TEMP_ADJUSTED','PSAL_ADJUSTED','DOXY_ADJUSTED','NITRATE_ADJUSTED',
                 'DIC','pH_25C_TOTAL_ADJUSTED','PDENS']


#initiate offset list
#number of additional rows for saving metadata items
num_meta_items = 10
float_offsets = [[] for _ in range(2*len(var_list_plot)+num_meta_items)]

#iterate over each interpolated float & profile
for wmo, group in argo_wmo:
    
    #number of profiles
    nprof = group.LATITUDE.shape[0]
    
     #check sufficient non-NaN data
    if group.num_var[0]<4:
        print('No non-NAN bgc adjusted data for: '+str(wmo))
        continue
    
    #  Find lat-lon limits within dist of float location
    lat_tol = dist/ 111.6 
    lon_tol = dist/ (111.6*np.cos(np.deg2rad(np.nanmean(group.LATITUDE))))  
    #set lat/lon crossover limits
    lat_min = group.LATITUDE.values-lat_tol
    lat_max = group.LATITUDE.values+lat_tol
    lon_min = group.LONGITUDE.values-lon_tol
    lon_max = group.LONGITUDE.values+lon_tol

    #find all data in lat-lon limits
    for n in range(nprof):
        
        #exclude profiles from same float
        testinds = argo_interp.wmo != wmo
        argo_wmo_test = argo_interp.wmo[testinds]
        argo_prof_test = argo_interp.profile[testinds]

        #index of all profiles within distance range
        match_inds = np.argwhere(np.logical_and(np.logical_and(argo_interp.LATITUDE[testinds].values>lat_min[n],
                                                            argo_interp.LATITUDE[testinds].values<lat_max[n]),
                            np.logical_and(argo_interp.LONGITUDE[testinds].values>lon_min[n],
                                           argo_interp.LONGITUDE[testinds].values<lon_max[n])))
        
        if match_inds.size==0:
            print('no float crossovers for float ' +str(wmo)+ ', profile ' +str(group.profile[n].values))
            continue

        #now load and loop through raw float data floats/profiles for matched test floats
        for m in np.arange(match_inds.shape[0]):
            w = argo_wmo_test[match_inds[m]].values[0]
            p = argo_prof_test[match_inds[m]].values[0]
            print('Comparing float '+str(w)+', profile '+str(p))
                
            #get single test profile data
            fn = argo_path+str(w)+'_adjusted.nc'
            test_float = xr.open_dataset(fn)
            prof_ind = np.argwhere(test_float.CYCLE_NUMBER.values==p)[0][0]
            

            #find index of interp profile density that most closely matches each test profile point density
            #only loop through valid data points
            p_valid = ~np.isnan(test_float.PDENS[prof_ind,:])
                                  
            test_pdens = test_float.PDENS[prof_ind,p_valid]
            
            if len(test_pdens):
 
                for i,pd in enumerate(test_pdens.values):
                    dens_ind = np.argmin(np.absolute(group.PDENS.values[n,:]-pd))
            
                    #don't use matches at top/bottom of interpolated vector
                    if (dens_ind == 0) or (dens_ind == len(p_interp) -1):
                        continue
                    #calc offset at the interp profile index for each match for each var to plot 
                    #check var is present in both float and glodap
                    for idx, var in enumerate(var_list_plot):
                        if var in group.keys() and var in test_float.keys():
                            #offset_ind = test_float[var].index[i]
                            offset = group[var][n,dens_ind] - test_float[var][prof_ind,i]
                                                                                                          
                            #append to offset list
                            float_offsets[idx*2].append(offset.values)
                    
                            #also save absoute test float profile value at crossover
                            float_offsets[idx*2+1].append(test_float[var][prof_ind,i].values)
                        
                        #append nan if variable is not there so lists all remain same length?
                        else:
                            float_offsets[idx*2].append(np.nan)
                            float_offsets[idx*2+1].append(np.nan)
                            
                #append metadata to offset list (from main and test float)
                float_offsets[len(var_list_plot)*2].append(wmo)
                float_offsets[len(var_list_plot)*2+1].append(w)
                float_offsets[len(var_list_plot)*2+2].append(group.profile[n].values)
                float_offsets[len(var_list_plot)*2+3].append(p)
                float_offsets[len(var_list_plot)*2+4].append(group.juld[n].values)
                float_offsets[len(var_list_plot)*2+5].append(test_float.JULD_LOCATION[prof_ind].values)
                float_offsets[len(var_list_plot)*2+6].append(group.LONGITUDE[n].values)
                float_offsets[len(var_list_plot)*2+7].append(test_float.LONGITUDE[prof_ind].values)
                float_offsets[len(var_list_plot)*2+8].append(group.LATITUDE[n].values)
                float_offsets[len(var_list_plot)*2+9].append(test_float.LATITUDE[prof_ind].values)
                #can add additional float metadata variable to list here


# +
#save float-float offsets to file
argo_offsets = xr.Dataset()
for idx, var in enumerate(var_list_plot):
    argo_offsets[var+'_offset'] = (float_offsets[idx*2])
    argo_offsets[var] = (float_offsets[idx*2+1])
    
argo_offsets['main_float_wmo'] = (float_offsets[len(var_list_plot)*2])
argo_offsets['test_float_wmo'] = (float_offsets[len(var_list_plot)*2+1])
argo_offsets['main_float_profile'] = (float_offsets[len(var_list_plot)*2+2])
argo_offsets['test_float_profile'] = (float_offsets[len(var_list_plot)*2+3])
argo_offsets['main_float_juld'] = (float_offsets[len(var_list_plot)*2+4])
argo_offsets['test_float_juld'] = (float_offsets[len(var_list_plot)*2+5])
argo_offsets['main_float_longitude'] = (float_offsets[len(var_list_plot)*2+6])
argo_offsets['test_float_longitude'] = (float_offsets[len(var_list_plot)*2+7])
argo_offsets['test_float_latitude'] = (float_offsets[len(var_list_plot)*2+8])
argo_offsets['test_float_latitude'] = (float_offsets[len(var_list_plot)*2+9])

argo_offsets.to_netcdf(output_dir+'float_offsets.nc')

print('Total number of float crossovers: ' + str(len(float_offsets[len(var_list_plot)*2])))
