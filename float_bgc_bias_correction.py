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

# -

# ### User inputs: 

# +
#User local directories
#argo_path = data_dir+'Sprof/' #USER LOCAL ARGO PATH !!!!! NOTE assumes user has Sprof files downloaded from Dropbox
#for now load fixed argo snapshot pre-downloaded Sprof files to make sure no updated QC 
#will add optionality to re-download updated/more recent Argo Sprof
#matlab_dir  = '/Users/veronicatamsitt/Documents/MATLAB/' #set paths for MATLAB LIAR/LIPHR (on local computer!)
#liar_dir = matlab_dir + 'LIRs-master/'

#pressure limits for interpolation
p_interp_min = 1200 #minimum pressure for float crossover comparison
p_interp_max = 2000 #maximum pressure for float crossover comparison
#pressure levels to interpolate to, every 1db
p_interp = np.arange(p_interp_min,p_interp_max+1)

#pressure limits for crossover comparison
p_compare_min = 1480
p_compare_max = 2000

#crossover distance range
dist = 50.

#variables to do crossover plot 
var_list_plot = ['TEMP_ADJUSTED','PSAL_ADJUSTED','DOXY_ADJUSTED','NITRATE_ADJUSTED',
                 'DIC','pH_25C_TOTAL_ADJUSTED','PDENS']
# -

# ## 1. Download and process GLODAP data

gdap = fl.get_glodap(data_dir, year = 2021)
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
    opt_buffers_mode=1, # used to be "buffers_mode='auto'" but seems to have changed in versions of pyco2?
)


gdap.pH_25C_TOTAL = results['pH_total_out']

#set pH to nan where there was no original pH data from GLODAP
gdap.pH_25C_TOTAL[np.isnan(gdap.G2phts25p0)]=np.nan
# -

#rename GLODAP comparison variables to match argo
gdap = gdap.rename(columns={'G2longitude':'LONGITUDE', 'G2latitude':'LATITUDE', 'G2pressure':'PRES_ADJUSTED',
                            'G2temperature':'TEMP_ADJUSTED','G2salinity':'PSAL_ADJUSTED', 
                            'G2oxygen':'DOXY_ADJUSTED','G2nitrate':'NITRATE_ADJUSTED', 'G2tco2':'DIC', 
                            'G2talk':'TALK_LIAR', 'G2MLD':'MLD','G2o2sat':'o2sat', 'G2PTMP':'PTMP', 
                            'pH_in_situ_total':'PH_IN_SITU_TOTAL_ADJUSTED','G2sigma0':'PDENS'})

# ## 2. Apply float bias corrections 

# +
append_data = 1 #reads in and adds to argo_interp_temp.nc rather than overwriting and running all floats
argolist = []
for file in os.listdir(argo_path):
    if file.endswith('Sprof.nc'):
        argolist.append(file)
LIAR_path = liar_dir

#float QC data fields
qc_data_fields = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED',  'PRES_ADJUSTED']

#variables to do crossover calculation
var_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 
            'pH_25C_TOTAL', 'PDENS', 'PRES_ADJUSTED', 'DIC']

#if append data is set to 1, reads in argo_interp_temp.nc which contains prior argo_interp array, 
#compare wmo numbers between argolist and the wmo numbers in argo_interp, and continues on processing 
#float files. otherwise, start from the beginning
if append_data==1 and os.path.exists(data_dir+'argo_interp_temp.nc'):
    #load previously saved argo_interp
    argo_interp = xr.load_dataset(data_dir+'argo_interp_temp.nc')

    #extract wmo #s as integers from argolist
    s = [int(s) for s in re.findall("[0-9]+", str(argolist_orig))]
    
    #find indices of argolist where argo_interp has matching wmos - don't need to run these again
    indices = [s.index(i) for i in s if i not in argo_interp.wmo]
    
    #set start_index to first argolist index not found in argo_interp
    start_index = indices[0]

else:
    #argolist=argolist_orig
    start_index=0            
                    
#####
#iterate through each float file 
#calculates derived carbonate parameters (currently TALK, DIC), bias corrected pH and stores all in argo_n
#interpolates all variables in "var_list" to 1 m resolution and stores in argo_interp_n
#saves out individual float netcdf files with variables to be adjusted/used for crossovers
#appends 1m interpolated dataset in argo_interp for comparison to glodap (not saved currently)
wmo_list= list()
for n in range(start_index,len(argolist)):
    print(str(n)+' Processing float file '+ argolist[n])
    argo_n = xr.load_dataset(argo_path+argolist[n])
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))

    wmo_n = argo_n.PLATFORM_NUMBER.values.astype(int)[0]
    wmo_list.append(wmo_n)
    
    nprof_n = argo_n.dims['N_PROF']
    
    #   set bad data and possibly bad data to NaN 
    for q in qc_data_fields:
        if q in argo_n.keys():
            qc_val = argo_n[q+'_QC'].values.astype('float')
            argo_n[q].where(np.logical_and(qc_val<3.,qc_val>4.))
            
    argo_n['PDENS'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
    argo_n.PDENS[:] = np.nan
    
    #initialise interpolated dataset
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
    for v in var_list:
        argo_interp_n[v] = (['N_PROF','N_LEVELS'],np.copy(nan_interp))
    
    #check first if PH_IN_SITU_TOTAL_ADJUSTED exists
    if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
        
        print('doing TALK, DIC and pH bias correction for float '+str(wmo_n))
        
        #initialise pH 25c and DIC variables - could do this only if float has pH
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
                opt_buffers_mode=1,
        )
           
        argo_n['pH_25C_TOTAL'] = (['N_PROF','N_LEVELS'],results['pH_total_out'])
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],results['dic'])  
        
        #is it necessary to loop through each float profile for co2sys or can we use apply_ufunc instead?
        for p in range(nprof_n):
            # skip a profile if pH is above 10.  There seem to be pH's above 10 that causing 
            #if any(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]>10) or all(np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:])):
            #    continue
    
            #apply pH bias correction   
            #find if there are valid values between fixed p levels 
            pH_p_min = 1480
            pH_p_max = 1520
    
            #if there are valid pressure levels between 1480-1520 db, 
            #calc bias correction only in this depth band, if not, calc correction between 970 and 1520
            if any((argo_n.PRES_ADJUSTED[p,:]>1480) & (argo_n.PRES_ADJUSTED[p,:]<1520)):
        
                inds = (argo_n.PRES_ADJUSTED[p,:]>1480) & (argo_n.PRES_ADJUSTED[p,:]<1520)
                correction = -0.034529*argo_n.pH_25C_TOTAL[p,inds]+0.26709
                              
            else:
                inds = (argo_n.PRES_ADJUSTED[p,:]>970) & (argo_n.PRES_ADJUSTED[p,:]<1520)
                correction = -0.034529*argo_n.pH_25C_TOTAL[p,inds]+0.26709
                              
            if len(correction):
                argo_n.bias_corr[p] = np.nanmean(correction)
                argo_n.pH_insitu_corr[p,:] = argo_n.PH_IN_SITU_TOTAL_ADJUSTED[p,:]+argo_n.bias_corr[p]
    

    
    ##### now calc potential density, save, and interpolate data for comparison
    for p in range(nprof_n):
        #pressure for profile
        p_prof = argo_n.PRES_ADJUSTED[p,:]
        
        #calc potential density 
        SA = gsw.SA_from_SP(argo_n.PSAL_ADJUSTED[p,:].values,
                            argo_n.PRES_ADJUSTED[p,:].values,
                            argo_n.LONGITUDE[p].values,
                            argo_n.LATITUDE[p].values)

        CT = gsw.CT_from_t(SA,
                           argo_n.TEMP_ADJUSTED[p,:].values,
                           argo_n.PRES_ADJUSTED[p,:].values)

        argo_n['PDENS'][p,:] = gsw.sigma0(SA,CT)
        

        #for each profile get pressure values >100db
        p100 = p_prof[p_prof>100.].values
            
        #if only 1 value of pressure or if there is not valid profile data down to p-max, continue loop
        if (len(p100) <= 1) or (np.nanmax(p100)<p_interp_min):
            continue
        
        #find which crossover variables exist in main float file
        var_list_n = []
        for vname in var_list:
            if (vname in argo_n.keys()) and (np.any(~np.isnan(argo_n[vname]))):
                var_list_n.append(vname)
                
        argo_interp_n['num_var'][p] = len(var_list_n) # changed to argo_interp_n from argo_interp
        
        for var in var_list_n:
            var100 = argo_n[var][p,p_prof>100.]

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
                (np.nanmax(p100u[~np.isnan(var100u.values)])>p_compare_min) and \
                (np.nanmin(p100u[~np.isnan(var100u.values)])<p_compare_max):
                
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
            
            else: 
                print('profile data not deep enough to interpolate')
    
    #save adjusted/processed, non-interpolated data to use for crossovers
    #could instead save each comparison var as a list of lists that can support different sizes r
    #rather than write to file?
    #or just output a dataset with only the relevant variables for crossovers?
    
    ###NOTE: currently getting a ValueError in to_netcdf for 1 float (2903176)
    #clean_dataset(argo_n)
    #argo_n.to_netcdf(argo_path+str(wmo_n)+'_adjusted.nc')
    #workaround: create new dataset with relevant crossover variables only
    argo_n_adjusted = xr.Dataset()
    argo_n_adjusted['wmo'] = wmo_n
    argo_n_adjusted['CYCLE_NUMBER'] = (['N_PROF'],argo_n.CYCLE_NUMBER.values)
    argo_n_adjusted['LONGITUDE'] = (['N_PROF'],argo_n.LONGITUDE.values)
    argo_n_adjusted['LATITUDE'] = (['N_PROF'],argo_n.LATITUDE.values)
    argo_n_adjusted['JULD_LOCATION'] = (['N_PROF'],argo_n.JULD_LOCATION.values)
    for var in var_list_plot:
        if var in argo_n.keys():
            argo_n_adjusted[var] = (['N_PROF','N_LEVELS'],argo_n[var].values)
    argo_n_adjusted.to_netcdf(argo_path+str(wmo_n)+'_adjusted.nc')
    
    #Also save interpolated dataset as one mega Dataset for doing crossovers
    if n == 0:
        argo_interp = argo_interp_n
    else:
        argo_interp = xr.concat([argo_interp,argo_interp_n],'N_PROF')
                    
    #save out argo_interp periodically:
    if n/20==round(n/20):
        argo_interp.to_netcdf(data_dir+'argo_interp_temp.nc')
     


# -

# ## 3. Compare float - GLODAP crossovers

# +
#float- GLODAP crossover 

#restrict glodap data to comparison pressure range
gdap_p = gdap[(gdap.PRES_ADJUSTED.values>p_compare_min) & (gdap.PRES_ADJUSTED.values<p_compare_max)]

#group by float wmo
argo_wmo = argo_interp.groupby('wmo')

#initiate offset list
#number of additional rows for saving metadata items
num_meta_items = 7
gdap_offsets =  [[] for _ in range(2*len(var_list_plot)+num_meta_items)]

#iterate over each float & profile
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
        
        #index of all gdap profiles within distance range
        if lon_min[n] < 0:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n]+360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]))
        elif lon_max[n] > 360:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n],
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
        for m in range(len(match_inds)):
            #if no interpolated density (not deep enough) then move on to next profile
            if np.all(np.isnan(group.PDENS.values[n,:])):
                continue
            dens_ind = np.nanargmin(np.absolute(group.PDENS.values[n,:]-gdap_match.PDENS.values[m]))
            
            #don't use matches at top/bottom of interpolated vector
            if (dens_ind == 0) or (dens_ind == len(p_interp) -1):
                continue
            #calc offset at the interp profile index for each match for each var to plot 
            #check var is present in both float and glodap
            for idx, var in enumerate(var_list_plot):
                if var in group.keys():
                    gdap_offset_ind = gdap_match[var].index[m]
                    offset = group[var][n,dens_ind] - gdap_match[var][gdap_offset_ind]
                    
                    #append to offset list
                    gdap_offsets[idx*2].append(offset.values)
                    
                    #also save absoute glodap value at crossover
                    gdap_offsets[idx*2+1].append(gdap_match[var][gdap_offset_ind])
                #append nan if variable is not there so lists all remain same length?
                else:
                    gdap_offsets[idx*2].append(np.nan) # changed to gdap_offsets from float_offsets
                    gdap_offsets[idx*2+1].append(np.nan) # changed to gdap_offsets from float_offsets
            
            #append metadata to offset list
            gdap_offsets[len(var_list_plot)*2].append(wmo)
            gdap_offsets[len(var_list_plot)*2+1].append(group.profile[n].values)
            gdap_offsets[len(var_list_plot)*2+2].append(group.juld[n].values)
            gdap_offsets[len(var_list_plot)*2+3].append(group.LONGITUDE[n].values)
            gdap_offsets[len(var_list_plot)*2+4].append(gdap_match.LONGITUDE[gdap_match.LONGITUDE.index[0]])
            gdap_offsets[len(var_list_plot)*2+5].append(group.LATITUDE[n].values)
            gdap_offsets[len(var_list_plot)*2+6].append(gdap_match.LATITUDE[gdap_match.LATITUDE.index[0]])
            #can add additional float metadata variable to list here


# +
#convert GLODAP offset lists to xarray Dataset and save to netcdf file
glodap_offsets = xr.Dataset()
for idx, var in enumerate(var_list_plot):
    glodap_offsets[var+'_offset'] = (gdap_offsets[idx*2])
    glodap_offsets[var] = (gdap_offsets[idx*2+1])
    
glodap_offsets['main_float_wmo'] = (gdap_offsets[len(var_list_plot)*2])
glodap_offsets['main_float_profile'] = ( gdap_offsets[len(var_list_plot)*2+1])
glodap_offsets['main_float_juld'] = (gdap_offsets[len(var_list_plot)*2+2])
glodap_offsets['main_float_longitude'] = (gdap_offsets[len(var_list_plot)*2+3])
glodap_offsets['glodap_longitude'] = (gdap_offsets[len(var_list_plot)*2+4])
glodap_offsets['main_float_latitude'] = (gdap_offsets[len(var_list_plot)*2+5])
glodap_offsets['glodap_latitude'] = (gdap_offsets[len(var_list_plot)*2+6])

glodap_offsets.to_netcdf(output_dir+'glodap_offsets.nc')

print('Total number of glodap crossovers: ' + str(len(gdap_offsets[len(var_list_plot)*2])))
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
# -

# ## Plot crossovers


# +
#now plot histograms of offsets for each main float with crossovers
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
# -
# compare to SOCCOM floats

