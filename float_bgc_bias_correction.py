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

# ## Apply BGC float bias correction and compare float/GLODAP crossovers

# Authors: Veronica Tamsitt (USF) ...
#
# Adapted from MALTAB code written by Seth Bushinsky (UH)
#
# 1. Download GLODAP data 
# 2. process data
# 3. apply float bias corrections and calculate derivative variables (pH, TALK)
# 4. do float/glodap crossover matchups

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

# ### User inputs: 

matlab_dir  = '/Users/veronicatamsitt/Documents/MATLAB/' #set paths for MATLAB LIAR/LIPHR (on local computer!)
liar_dir = matlab_dir + 'LIRs-master/'
p_interp_min = 1200 #minimum pressure for float crossover comparison
p_interp_max = 2000 #maximum pressure for float crossover comparison

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

# ## 1. Download GLODAP data

gdap = fl.get_glodap(data_dir, year = 2021)

# ## 2. Data processing

gdap.G2longitude[gdap.G2longitude < 0.] = gdap.G2longitude[gdap.G2longitude < 0.] + 360.


# GLODAP derived variables: density, MLD and pH

# +
#set flagged data to NaN (is this needed? or masked array better?)
flagvars = ['G2salinity','G2oxygen','G2nitrate','G2tco2','G2talk','G2phts25p0']

for v in flagvars:
    flag = v+'f'
    naninds = gdap[flag]!=2
    gdap[v][naninds] = np.nan

# +
###make into separate function?
#iterate over grouped data to calc MLD
#need to group by cruise and station first
gdap_s = gdap.groupby(['G2cruise','G2station','G2cast'])
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
#pH from LIPHR
# calculate LIPHR pH at Glodap points below 1480 m and above 2020m (V: where does the depth restriction come in?)
LIPHR_path = liar_dir
Coordinates = np.stack((gdap.G2longitude.values.flatten(), 
                        gdap.G2latitude.values.flatten(), 
                        gdap.G2pressure.values.flatten()),
                        axis=1)
Measurements = np.stack((gdap.G2salinity.values.flatten(), 
                         gdap.G2temperature.values.flatten(), 
                         gdap.G2nitrate.values.flatten(), 
                         gdap.G2oxygen.values.flatten()),
                         axis=1)
MeasIDVec = [1, 7, 3, 6]
                                
results = carbon_utils.LIPHR_matlab(LIPHR_path,
                                    Coordinates.tolist(),
                                    Measurements.tolist(),
                                    MeasIDVec, 
                                    OAAdjustTF = False)                                  

gdap['pH_in_situ_total'] = results
gdap.pH_in_situ_total[np.isnan(gdap.G2phts25p0)] = np.nan

# +
# gdap pH 25C -> Nancy or Seth to double check this matches the MATLAB SOCCOM version
#   *SOCCOM* version modified by Nancy Williams on 10/15/15 according to
#    Dickson in 9/7/15 e-mail and in Dickson et al. 2007 
#    changed KF to Perez and Fraga 1987
#    Last three inputs should be ... 1,10,3)
results = pyco2.sys(
    par1=2300., 
    par2=gdap.pH_in_situ_total,
    par1_type=1,
    par2_type=3,
    temperature=gdap.G2temperature, 
    pressure=gdap.G2pressure, 
    salinity=gdap.G2salinity, 
    temperature_out=25., #fixed 25C temperature
    pressure_out=gdap.G2pressure,
    opt_pH_scale = 1, #total
    opt_k_carbonic=10, #Lueker et al. 2000
    opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
    opt_total_borate=2, # Lee et al. 2010
    opt_k_fluoride=2, # Perez and Fraga 1987
    buffers_mode='auto',
)


gdap.pH_25C_TOTAL = results['pH_total_out']

#set pH to nan where there was no original pH data from GLODAP
gdap.pH_25C_TOTAL[np.isnan(gdap.G2phts25p0)]=np.nan
# -

# ## 3. Apply float bias corrections 

# +
##for now load fixed argo snapshot pre-downloaded Sprof files to make sure no updated QC 
#will add optionality to re-download updated/more recent Argo Sprof
#!!!!! NOTE assumes user has Sprof file downloaded from shared Dropbox
argo_path = data_dir+'Sprof/' #USER LOCAL ARGO PATH
argolist = os.listdir(argo_path)
LIAR_path = liar_dir

qc_data_fields = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED',  'PRES_ADJUSTED']

#iterate through each float file 
#-> is this the best way? can't easily combine floats into 1 dataset unless interpolated onto same p levels, 
#so do initial processing in loop and then combine into a single dataset once interpolated onto p levels?
wmo_list= list()
for n in argolist:
    argo_n = xr.open_dataset(argo_path+n)
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))

    wmo_n = argo_n.PLATFORM_NUMBER.values.astype(int)[0]
    wmo_list.append(wmo_n)
    
    #check first if PH_IN_SITU_TOTAL_ADJUSTED exists
    if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
        
        print('doing TALK, DIC and pH bias correction for float '+str(wmo_n))
        
        #initialise pH 25c and DIC variables
        argo_n['TALK_LIAR'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.TALK_LIAR[:] = np.nan
        argo_n['pH_25C_TOTAL'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.pH_25C_TOTAL[:] = np.nan
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.DIC[:] = np.nan
        argo_n['pH_insitu_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.pH_insitu_corr[:] = np.nan
        argo_n['bias_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof
        argo_n.bias_corr[:] = np.nan
    
        #   set bad data and possibly bad data to NaN 
        for q in qc_data_fields:
            qc_val = argo_n[q+'_QC'].values.astype('float')
            argo_n[q].where(np.logical_and(qc_val<3.,qc_val>4.))

        
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
            SI = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape()))
            PO4 = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape()))
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

        argo_n['TALK_LIAR'] =   (['N_PROF','N_LEVELS'],np.reshape(np.asarray(results),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
  
    
        ##### Calculate float pH at 25C, DIC and apply bias corr
        results = pyco2.sys(
                par1=argo_n.TALK_LIAR, 
                par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
                par1_type=1,
                par2_type=3,
                temperature=argo_n.TEMP_ADJUSTED, 
                pressure=argo_n.PRES_ADJUSTED, 
                salinity=argo_n.PSAL_ADJUSTED, 
                temperature_out=25., #fixed 25C temperature
                pressure_out=argo_n.PRES_ADJUSTED,
                total_silicate=SI,
                total_phosphate=PO4,
                opt_pH_scale = 1, #total
                opt_k_carbonic=10, #Lueker et al. 2000
                opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
                opt_total_borate=2, # Lee et al. 2010
                opt_k_fluoride=2, # Perez and Fraga 1987
                buffers_mode='auto',
        )
           
        argo_n['pH_25C_TOTAL'] = (['N_PROF','N_LEVELS'],results['pH_total_out'])
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],results['dic'])  
        
        #is it necessary to loop through each float profile for co2sys or can we use apply_ufunc instead?
        for p in range(argo_n.PRES_ADJUSTED.shape[0]):
            # skip a profile if pH is above 10.  There seem to be pH's above 10 that causing 
            #if any(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]>10) or all(np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:])):
            #    continue
    
            #apply pH bias correction   
            #find if there are valid values between fixed p levels 
            pH_p_min = 1480
            pH_p_max = 1520
    
            #V- not 100% confident I understood this correctly- if there are valid pressure levels between 1480-1520 db, 
            #calc bias correction only in this depth band, if not, calc correction between 970 and 1520
            if any(argo_n.PRES_ADJUSTED[p,:]>1480 & argo_n.PRES_ADJUSTED[p,:]<1520):
        
                inds = argo_n.PRES_ADJUSTED[p,:]>1480 & argo_n.PRES_ADJUSTED[p,:]<1520
                correction = -0.034529*argo_n.pH_25C_TOTAL[inds]+0.26709
                              
            else:
                inds = argo_n.PRES_ADJUSTED[p,:]>970 & argo_n.PRES_ADJUSTED[p,:]<1520
                correction = -0.034529*argo_n.pH_25C_TOTAL[inds]+0.26709
                              
            if len(correction):
                argo_n.bias_corr[p] = np.nanmean(correction)
                argo_n.pH_insitu_corr[p,:] = argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]+argo_n.bias_corr[p]
    
    #now interpolate to 1db pressure levels for crossover comparison
    #pressure levels to interpolate to, every 1db
    p_interp = np.arange(p_interp_min,p_interp_max+1) 
    

    
    #for each profile get pressure values(indices) >100db
    #if only 1 value of pressure, continue loop
    #for comparison vars make sure at least one datapoint for p>100db
    #if there are non-unique pressure values, then grab only unique pressure values and matching data points
    #interpolate 1d onto p_interp levels (non-NaN data as input only)
    #calc potential density on interpolated p and T (gsw_pdens)
    
    #at this point, can merge float interpolated data to one Dataset all with same p levels and exit mega loop?
# -

# ## 4. Compare float/GLODAP crossovers

# +
#make this into a separate script and call data processing from within? Can input all user variables at the top?
#variables to do crossover plots
var_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 
            'pH_25C_TOTAL', 'PDENS', 'PRES_ADJUSTED', 'DIC']

#crossover distance range
dist = 50.





#loop through floats to find crossovers

#find which crossover variables exist in main float file
for var in var_list:
    if var in argo_n.keys() and any(~np.isnan(argo_n[var])):
        var_list_n.append(var)
    
if len(var_list_n)<4:
    print('No non-NAN bgc adjusted data for: '+wmo_n)
    continue

#fl.float_float_crossovers(floatn,varlist,dist,p_interp)
# -




