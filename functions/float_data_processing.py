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

# ## Download and process Argo float data

import numpy as np
import glob, os
import pandas as pd
import xarray as xr
import gsw
import ftplib
import getpass


def get_glodap(save_dir, year):
    if year==2021:
        url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0237935/GLODAPv2.'+str(year)+'_Merged_Master_File.csv'
    elif year==2022:
        url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0257247/GLODAPv2.2022_Merged_Master_File.csv'
    elif year==2023:
        url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0283442/GLODAPv2.2023_Merged_Master_File.csv'
        
    print(url)
    if not os.path.exists(save_dir+'GLODAPv2.'+str(year)+'_Merged_Master_File.csv'):
        os.system('wget %s -P %s' % (url,save_dir))
    gdap = pd.read_csv(save_dir+'GLODAPv2.'+str(year)+'_Merged_Master_File.csv')
    
    #add datetime
    df = pd.DataFrame({'year': gdap['G2year'],
                   'month': gdap['G2month'],
                   'day': gdap['G2day'],
                   'hour': gdap['G2hour'],
                   'minute': gdap['G2minute']})
    gdap['datetime'] = pd.to_datetime(df)
    
    return gdap


#argo ftp download
def get_argo(save_dir, wmo_list):
    IP = 'ftp.ifremer.fr'
    filepath = './ifremer/argo/dac/'

    ftp = ftplib.FTP(IP,user='anonymous')
    ftp.cwd(filepath)
    daclist = ftp.nlst()
    for dac in daclist:
        ftp.cwd(dac+'/')
        wmos = ftp.nlst() #get list of wmos
        print(wmos)
        #for w in wmos:
        #    if 
        ftp.cwd('..')

    ftp.quit()


# +
def float_float_crossovers(floatn,test_floats,varlist,dist,p_interp):

    #find any other float profile that crosses this float and save the float in the
    #match list, along with the profiles that fall within the lat and lon test
    
    #ignore any profiles with <4 data points

    #  Find lat-lon limits within dist of float location
    lat_tol = dist/ 111.6 
    lon_tol = dist/ (111.6*np.cosd(np.nanmean(floatn.LATITUDE)))  

    last = 1 #why this?

    #set lat/lon crossover limits
    lat_min = floatn.LATITUDE-lat_tol
    lat_max = floatn.LATITUDE+lat_tol
    lon_min = floatn.LONGITUDE-lon_tol
    lon_max = floatn.LONGITUDE+lon_tol

    #for each float (og float), compare with all other floats (test floats)
    #in lat/lon range of each og float profile, find 
    #1) other og float profiles in this range and 
    #2) test float profiles in this range

# for each of the matched profiles from the main float

# +
#def float_glodap_crossovers():
# -


