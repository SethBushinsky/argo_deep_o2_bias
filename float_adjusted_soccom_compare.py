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

# ### Compare adjusted BGC float pH, DIC, pCO2 to SOCCOM published data

import numpy as np
import glob, os
from pathlib import Path
from datetime import datetime, date, time
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import colorbar, colors

# +
save_dir = './output/SOCCOM_figs'

#check directories exist
if not os.path.isdir(save_dir):
    os.mkdir(save_dir)
# -

# 1. SOCCOM data snapshot

soccom_dir = './data/SOCCOM/SOCCOM_HiResQC_LIAR_19May2022_netcdf/'
soccom_list = glob.glob(soccom_dir+'*.nc')
soccom_list.sort()

# 2. Our adjusted float files

float_dir = './data/argo_Sprof/'
f_list = glob.glob(float_dir+'*_adjusted.nc')

# load SOCCOM floats one by one and compare to our adjusted float data

# +
qc_data_fields = ['Nitrate','Oxygen','DIC_LIAR']

for s in soccom_list:
    print(s)
    #get wmo
    s_wmo  = s[51:58]

    #see if matching wmo in our adjsted files
    res = [i for i in f_list if s_wmo in i]

    if len(res)==0:
        print('no adjusted DIC/pCO2 data')
        continue
        
    #load both data files
    s_float = xr.open_dataset(s)
    a_float = xr.open_dataset(res[0])
    
    #apply soccom qc
    #   set bad data and possibly bad data to NaN 
    for q in qc_data_fields:
        if q in s_float.keys():
            qc_val = s_float[q+'_QF'].values.astype('float')
            s_float[q].where(np.logical_and(qc_val<3.,qc_val>4.))
    
    for prof in range(s_float.JULD.shape[0]):
            
        #sp_ind = np.flatnonzero(~np.isnan(s_p[prof,:]))[0]
        #ap_ind = np.flatnonzero(~np.isnan(a_float.PRES_ADJUSTED[prof,:]))[-1]+1
            
        s_pressure = s_float.Pressure[prof,::-1]
        a_pressure = a_float.PRES_ADJUSTED[prof,:]
            
        fig,axs = plt.subplots(1,3,figsize=(12,4))
        
        #compare nitrate
        if 'NITRATE_ADJUSTED' in a_float.keys():
            s_nitrate = s_float.Nitrate[prof,::-1]
            a_nitrate = a_float.NITRATE_ADJUSTED[prof,:]
            axs[0].scatter(s_nitrate,s_pressure,s=5,c='orangered',label='SOCCOM')
            axs[0].scatter(a_nitrate,a_pressure,s=5,c='blue',label='our data')
            plt.gca().invert_yaxis()
            axs[0].set_xlabel('nitrate (umol/kg)')
            axs[0].set_ylabel('pressure (db)')
            axs[0].legend()
            axs[0].set_title('nitrate float ' + s_wmo +', profile '+str(prof))
            #nitrate_diff = s_nitrate[sp_ind:] - a_nitrate[:ap_ind]
            
        #compare oxygen
        if 'DOXY_ADJUSTED' in a_float.keys():
            s_doxy = s_float.Oxygen[prof,::-1]
            a_doxy = a_float.DOXY_ADJUSTED[prof,:]
            axs[1].scatter(s_doxy,s_pressure,s=5,c='orangered',label='SOCCOM')
            axs[1].scatter(a_doxy,a_pressure,s=5,c='blue',label='our data')
            plt.gca().invert_yaxis()
            axs[1].set_xlabel('oxygen (umol/kg)')
            axs[1].set_ylabel('pressure (db)')
            axs[1].legend()
            axs[1].set_title('O2 float ' + s_wmo +', profile '+str(prof))
            #doxy_diff = s_float.Oxygen - a_float.DOXY_ADJUSTED
    
        
        #compare DIC
        if 'DIC' in a_float.keys():
            s_dic = s_float.DIC_LIAR[prof,::-1]
            a_dic = a_float.DIC[prof,:]
            axs[2].scatter(s_dic,s_pressure,s=5,c='orangered',label='SOCCOM')
            axs[2].scatter(a_dic,a_pressure,s=5,c='blue',label='our data')
            plt.gca().invert_yaxis()
            axs[2].set_xlabel('DIC (umol/kg)')
            axs[2].set_ylabel('pressure (db)')
            axs[2].legend()
            axs[2].set_title('DIC float ' + s_wmo +', profile '+str(prof))
        
        plt.savefig(save_dir+'/'+s_wmo+'_profile'+str(prof)+'_vsSOCCOM.png')

            #dic_diff = s_float.DIC_LIAR - a_float.DIC
# -


