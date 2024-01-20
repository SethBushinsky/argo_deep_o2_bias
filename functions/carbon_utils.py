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

def LIPHR_matlab(LIPHR_path,Coordinates,Measurements,MeasIDVec,OAAdjustTF=False,VerboseTF=False):
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
    results = eng.LIPHR(Coordinates,Measurements,MeasIDVec,'OAAdjustTF', OAAdjustTF, 'VerboseTF', VerboseTF)
    eng.quit()

    results = np.asarray(results)   
    return results


def ESPER_mixed_matlab(LIPHR_path,DesiredVariables,Coordinates,Measurements,MeasIDVec_ESPER,Equations, Dates, VerboseTF, pHCalcTF):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec_ESPER = matlab.double([MeasIDVec_ESPER])
    Equations = matlab.double([Equations])
    Dates = matlab.double([Dates])
    DesiredVariables = matlab.double([DesiredVariables])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)

    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIPHR_path))

    #call MATLAB function
    results = eng.ESPER_Mixed_wrapper_for_python(DesiredVariables,Coordinates,Measurements,MeasIDVec_ESPER,Equations,
                        Dates, VerboseTF, pHCalcTF)
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

def LINR_matlab(LIAR_path,Coordinates,Measurements, MeasIDVec,VerboseTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)
    
    #need to make sure LIR subfolders added to matlab path
    eng.addpath(eng.genpath(LIAR_path))

    #call MATLAB function
    results = eng.LINR(Coordinates,Measurements,MeasIDVec,'VerboseTF',VerboseTF)
    eng.quit()

    #convert matlab double output back to numpy array
    results = np.asarray(results)

    return results

def sigma0(salinity,temperature,lon,lat,pressure):
    SA = gsw.SA_from_SP(salinity,
                        pressure,
                        lon,
                        lat)

    CT = gsw.CT_from_t(SA,
                       temperature,
                       pressure)

    sigma = gsw.sigma0(SA,CT)
    
    return sigma


def spiciness0(salinity,temperature,lon,lat,pressure):
    SA = gsw.SA_from_SP(salinity,
                        pressure,
                        lon,
                        lat)

    CT = gsw.CT_from_t(SA,
                       temperature,
                       pressure)

    spiciness = gsw.spiciness0(SA,CT)
    
    return spiciness


def co2sys_pH25C(TALK,pH,temperature,salinity,pressure):
    #   *SOCCOM* version modified by Nancy Williams on 10/15/15 according to
    #    Dickson in 9/7/15 e-mail and in Dickson et al. 2007 
    #    changed KF to Perez and Fraga 1987
    #    Last three inputs should be ... 1,10,3)
    results = pyco2.sys(
        par1=TALK, 
        par2=pH,
        par1_type=1,
        par2_type=3,
        temperature=temperature, 
        pressure=pressure, 
        salinity=salinity, 
        temperature_out=25., #fixed 25C temperature
        pressure_out=pressure,
        opt_pH_scale = 1, #total
        opt_k_carbonic=10, #Lueker et al. 2000
        opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
        opt_total_borate=2, # Lee et al. 2010
        opt_k_fluoride=2, # Perez and Fraga 1987
        opt_buffers_mode=1, # used to be "buffers_mode='auto'" but seems to have changed in versions of pyco2?
    )

    pH_25C = results['pH_total_out']
    
    return pH_25C


def co2sys_DIC(TALK,pH,temperature,salinity,pressure):
    #   *SOCCOM* version modified by Nancy Williams on 10/15/15 according to
    #    Dickson in 9/7/15 e-mail and in Dickson et al. 2007 
    #    changed KF to Perez and Fraga 1987
    #    Last three inputs should be ... 1,10,3)
    results = pyco2.sys(
        par1=TALK, 
        par2=pH,
        par1_type=1,
        par2_type=3,
        temperature=temperature, 
        pressure=pressure, 
        salinity=salinity, 
        temperature_out=25., #fixed 25C temperature
        pressure_out=pressure,
        opt_pH_scale = 1, #total
        opt_k_carbonic=10, #Lueker et al. 2000
        opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
        opt_total_borate=2, # Lee et al. 2010
        opt_k_fluoride=2, # Perez and Fraga 1987
        opt_buffers_mode=1, # used to be "buffers_mode='auto'" but seems to have changed in versions of pyco2?
    )

    DIC = results['dic']
    
    return DIC
