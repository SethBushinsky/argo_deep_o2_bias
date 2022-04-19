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

# User input: set paths for MATLAB LIAR/LIPHR (on local computer!)

matlab_dir  = '/Users/veronicatamsitt/Documents/MATLAB/'
liar_dir = matlab_dir + 'LIRs-master/'

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
# -

#need to group by cruise and station first
gdap_s = gdap.groupby(['G2cruise','G2station','G2cast'])

# +
###make into separate function?
#iterate over grouped data to calc MLD
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
                                
results = carbon_utils.LIPHR_matlab(LIPHR_path,Coordinates.tolist(),Measurements.tolist(),MeasIDVec, OAAdjustTF = False)                                  

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
for n in argolist[:1]:
    argo_n = xr.open_dataset(argo_path+n)
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))

    wmo_n = argo_n.PLATFORM_NUMBER
    wmo_list.append(wmo_n)
    
    
    #   set bad data and possibly bad data to NaN or mask? Let's try mask
    for q in qc_data_fields:
        qc_val = argo_n[q+'_QC'].values
        print(qc_val)
        argo_n[q].where(qc_val.decode()<3. and qc_val.decode()>4.)
        #argo_n[q][badind] = np.nan
        #badind = argo_n[q+'_QC'] == 4
        #argo_n[q][badind] = np.nan

    ##### Calc float TALK   
    #set Si and PO4 inputs
    #if nitrate, then use redfield for Si and PO4?, otherwise set to 0    
    if 'NITRATE_ADJUSTED' in argo_n.keys():
        SI = argo_n.NITRATE_ADJUSTED*2.5
        SI[np.isnan(SI)] = 0
        PO4 = argo_n.NITRATE_ADJUSTED/16
        PO4[isnan(PO4)] = 0
        Coordinates = np.stack((argo_n.LONGITUDE.values.flatten(), 
                        argo_n.LATITUDE.values.flatten(), 
                        argo_n.PRESSURE.values.flatten()),
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
        Coordinates = np.stack((argo_n.LONGITUDE.values.flatten(), 
                        argo_n.LATITUDE.values.flatten(), 
                        argo_n.PRESSURE.values.flatten()),
                        axis=1)
        Measurements = np.stack((argo_n.PSAL_ADJUSTED.values.flatten(), 
                         argo_n.TEMP_ADJUSTED.values.flatten(),
                         argo_n.DOXY_ADJUSTED.values.flatten()),
                         axis=1)
        MeasIDVec = [1, 7, 6]                            


    results = carbon_utils.LIAR_matlab(LIAR_path,Coordinates.tolist(),Measurements.tolist(),MeasIDVec,'VeroboseTF',False)                                  
    print(results)
    
# argo_n.TALK_LIAR =   reshape(AlkalinityEstimates,size(Argo.(SNs{f}).PRES_ADJUSTED,1),size(Argo.(SNs{f}).PRES_ADJUSTED,2));
  
    
    ##### Calculate float pH at 25C, DIC and apply bias corr
    
    #need to initialise pH 25c and DIC variables empty?
    #is it necessary to loop through each float profile for co2sys or can we use apply_ufunc?
    for p in range(argo_n.PRES_ADJUSTED.shape[0]):
        # skip a profile if pH is above 10.  There seem to be pH's above 10 that causing 
        # CO2SYS to hang up and probably not even worth considering otherwise
        if any(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]>10) or all(np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:])):
            continue

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
            total_phosphae=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            buffers_mode='auto',
        )
           
        argo_n.pH_25C_TOTAL[p,:] = results[:,37]
        argo_n.DIC[p,:] = results[:,2]
    
        #apply pH bias correction   
        #find if there are valid values between fixed p levels 
        pH_p_min = 1480
        pH_p_max = 1520
    
        #V- not confident I understood this correctly- if there are valid pressure levels between 1480-1520, 
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
# -

#



# ## 4. Compare float/GLODAP crossovers

# +
#make this into a separate script and call data processing from within? Can input all user variables at the top?
#variables to do crossover plots
var_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 
            'pH_25C_TOTAL', 'PDENS', 'PRES_ADJUSTED', 'DIC']

#crossover distance range
dist = 50.

#pressure levels to interpolate to, every 1db
p_interp = np.arange(1200,2001) 



#loop through floats to find crossovers
#fl.float_float_crossovers(floatn,varlist,dist,p_interp)
# -




