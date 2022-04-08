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

# ## Download and process Argo float data

import numpy as np
import glob, os
import pandas as pd
import xarray as xr
import gsw
import ftplib


def get_glodap(save_dir, year=2021):
    url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0237935/GLODAPv2.'+str(year)+'_Merged_Master_File.csv'
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

# +
#def get_argo_bgc_sprof(variables,data_dir,save_dir):
#data_dir = './data/'
#variables = ['NO3','O2','pH']

##get float lists -> eventually want to phase out to get float list directly from DACs, not rely on these saved .csv files
##actually is there a way to do this with Python that never involves needing a list/metadata list? 
##Can we use xarray to make dataset of all floats from ftp then only read in ones with pH/NO3/O2? or use argopy?
#if 'NO3' in variables:
#    filelist = 'Ocean_Ops_NO3_list2.csv'
#    NO3_wmo = pd.read_csv(data_dir+filelist)['REF'].astype('string')
#if 'O2' in variables:
#    filelist = 'Ocean_Ops_O2_list2.csv'
#   O2_wmo = pd.read_csv(data_dir+filelist)['REF'].astype('string')
#
#if 'pH' in variables:
#    filelist = 'Ocean_Ops_pH_list2.csv'
#    pH_wmo = pd.read_csv(data_dir+filelist)['REF'].astype('string')

#wmo_all = np.intersect1d(NO3_wmo,np.intersect1d(O2_wmo,pH_wmo))

#metadata file to get filepaths- is this needed or can we search the ftp file list for matches for each wmo?
#metafn = 'ar_index_global_meta.txt'
#meta = pd.read_csv(data_dir+metafn,header=8)

#file_dir = 'ftp://ftp.ifremer.fr/ifremer/argo/dac/' + partial_dir + wmo + '_Sprof.nc'

# +
#def float_float_crossovers():

# +
#def float_glodap_crossovers():
# -


