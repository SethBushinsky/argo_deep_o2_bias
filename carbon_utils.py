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

# ## Functions fo carbonate system calculations (pH, TALK etc)

import numpy as np
import glob, os
import pandas as pd
import xarray as xr
import gsw
import PyCO2SYS as pyco2
import matlab.engine


# ### LIAR/LIPHR wrapper

# +
def LIPHR_matlab(LIPHR_path,Coordinates,Measurements, MeasIDVec,EstDates,Equations,OAAdjustTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double(Measurements)
    Coordinates = matlab.double(OutputCoordinates.tolist())
    MeasIDVec = matlab.double(MeasIDVec)

    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIAR_path))

    #call MATLAB function
    pHEstimates,UncertaintyEstimates,MinUncertaintyEquation = eng.LIPHR(Coordinates,Measurements,MeasIDVec, 'OAAdjustTF', OAAdjustTF)
    eng.quit()

    #convert matlab double output back to numpy array
    pH_estimates = np.asarray(pHEstimates)

   
    #reshape DIC/TA back to profiles
    DIC = np.reshape(DIC_all,(len(path_list_p),2628))
    DIC_unc = np.reshape(DIC_unc_all,(len(path_list_p),2628))
    TA = np.reshape(TA_all,(len(path_list_p),2628))
    TA_unc = np.reshape(TA_unc_all,(len(path_list_p),2628))

return()


