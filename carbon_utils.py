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

# ## Functions for carbonate system calculations (pH, TALK etc)

import numpy as np
import glob, os
import pandas as pd
import xarray as xr
import gsw
import PyCO2SYS as pyco2
import matlab.engine


# ### LIAR/LIPHR wrapper

def LIPHR_matlab(LIPHR_path,Coordinates,Measurements,MeasIDVec,OAAdjustTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])
    
    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)

    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIPHR_path))

    #call MATLAB function
    results = eng.LIPHR(Coordinates,Measurements,MeasIDVec,'OAAdjustTF', OAAdjustTF)
    eng.quit()

    results = np.asarray(results)   
    return results


def LIAR_matlab(LIAR_path,Coordinates,Measurements, MeasIDVec,VerboseTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)
    
    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIAR_path))

    #call MATLAB function
    results = eng.LIAR(Coordinates,Measurements,MeasIDVec,'VerboseTF',VerboseTF)
    eng.quit()

    #convert matlab double output back to numpy array
    results = np.asarray(results)

    return results
