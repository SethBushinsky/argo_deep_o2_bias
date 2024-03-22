import xarray as xr
import numpy as np
import functions.carbon_utils as carbon_utils
import PyCO2SYS as pyco2

# float_wmo_list = glodap_offsets_p.index.values
# MeasIDVec = [1, 7, 3, 6]
# MeasIDVec_LIR = [1, 7, 6]
# DesiredVariables = [3]
# MeasIDVec_ESPER = [1, 2, 6] # S, T, O2 - different numbering than v2 LIRs
# Equations = 7 # for ESPER - asking to use equation w/ S, T, and O2 only 


def calc_derived_2_pH_no3_impacts(LIPHR_path, MeasIDVec_LIR, MeasIDVec_ESPER, DesiredVariables, Equations, argo_path_derived, argo_path_interpolated, wmo_n, glodap_offsets_p):
    # wmo_n = float_wmo_list[n]
    if glodap_offsets_p.pH_group[glodap_offsets_p.index==wmo_n].values=='no pH':
        return
    
    # if oxygen offset is nan, then skip
    if np.isnan(glodap_offsets_p.DOXY_ADJUSTED_offset_trimmed[glodap_offsets_p.index==wmo_n].values):
        return

    print(str(wmo_n) + ' started')
    
    argo_derived_n = xr.load_dataset(argo_path_derived+ str(wmo_n) + '_derived.nc')

    nprof = argo_derived_n.LATITUDE.shape[0]
    if nprof==1: # if there is only one profile, probably don't want to trust any offset
        return

    # Check for presence of non-nan pH
    if ~np.any(~np.isnan(argo_derived_n['PH_IN_SITU_TOTAL_ADJUSTED'])):
        return
    
    # load interpolated file
    argo_interp_n = xr.load_dataset(argo_path_interpolated + str(wmo_n) + '_interpolated.nc')

    # initialize coordinate and measurement arrays
    Coordinates_all = np.empty([nprof, 3],dtype=float)
    Coordinates_all[:] = np.nan

    Measurements_all_S_T_O2 = np.empty([nprof, 3],dtype=float)
    Measurements_all_S_T_O2[:] = np.nan

    Measurements_all_S_T = np.empty([nprof, 2],dtype=float)
    Measurements_all_S_T[:] = np.nan
    for p in range(nprof):
        # load interpolated argo_data so that 1500 m data is present 
        argo_profile = argo_interp_n.isel(N_PROF=p)

        # argo_profile = argo_interp_wmo[wmo_n].isel(N_PROF=p)
        argo_profile.load()
        data_1500 = argo_profile.where(argo_profile.PRES_ADJUSTED==1500, drop=True) 
        if data_1500.LATITUDE.size==0:
            continue

        Coordinates_all[p,:] = np.stack((data_1500.LONGITUDE.values, 
                            data_1500.LATITUDE.values, 
                            data_1500.PRES_ADJUSTED.values),
                            axis=1)
        
        Measurements_all_S_T_O2[p,:] = np.stack((data_1500.PSAL_ADJUSTED.values, 
                            data_1500.TEMP_ADJUSTED.values, 
                            data_1500.DOXY_ADJUSTED.values),
                            axis=1)
        
        Measurements_all_S_T[p,:] = np.stack((data_1500.PSAL_ADJUSTED.values, 
                        data_1500.TEMP_ADJUSTED.values),
                        axis=1)
    if np.sum(~np.isnan(Coordinates_all))==0: # if no valid 1500 meter data is present, skip 
        return

    # calculate decimal_year for ESPER
    da = argo_interp_n.juld
    decimal_year = da.dt.year + (da.dt.dayofyear - 1 + (da.dt.hour * 3600 + da.dt.minute * 60 + da.dt.second) / 86400) / (365 + da.dt.is_leap_year)

    # calculate LIPHR pH at 1500 m using original data from all profiles
    orig_pH = carbon_utils.LIPHR_matlab(LIPHR_path,
                                    Coordinates_all.tolist(),
                                    Measurements_all_S_T_O2.tolist(),
                                    MeasIDVec_LIR, 
                                    OAAdjustTF = False)  
    
    # remove data when oxygen is NaN - avoids biasing difference low
    pH_orig_no_nan = np.copy(orig_pH)
    pH_orig_no_nan[np.isnan(Measurements_all_S_T_O2[:,2])] = np.nan

    # save 1500 pH orig
    argo_derived_n['pH_1500_orig_LIPHR'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
    argo_derived_n.pH_1500_orig_LIPHR[:] = pH_orig_no_nan[:,0]

    # calculate ESPER pH at 1500 m using original data from all profiles
    orig_pH_ESPER = carbon_utils.ESPER_mixed_matlab(LIPHR_path,
                                                    DesiredVariables,
                                                    Coordinates_all.tolist(),
                                                    Measurements_all_S_T_O2.tolist(),
                                                    MeasIDVec_ESPER,
                                                    Equations, 
                                                    decimal_year.values.tolist(), 
                                                    0, 0)
    # remove data when oxygen is NaN - avoids biasing difference low
    pH_orig_ESPER_no_nan = np.copy(orig_pH_ESPER)
    pH_orig_ESPER_no_nan[np.isnan(Measurements_all_S_T_O2[:,2])] = np.nan

    # save 1500 pH orig
    argo_derived_n['pH_1500_orig_ESPER'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
    argo_derived_n.pH_1500_orig_ESPER[:] = pH_orig_ESPER_no_nan[:,0]
   
    Equations_ESPER_ST = 8 # for ESPER - asking to use equation w/ S, T 
    MeasIDVec_ESPER_ST = [1, 2] # S, T - different numbering than v2 LIRs

    # calculate ESPER pH at 1500 m without using oxygen
    pH_ESPER_wo_O2 = carbon_utils.ESPER_mixed_matlab(LIPHR_path,
                                                    DesiredVariables,
                                                    Coordinates_all.tolist(),
                                                    Measurements_all_S_T.tolist(),
                                                    MeasIDVec_ESPER_ST,
                                                    Equations_ESPER_ST, 
                                                    decimal_year.values.tolist(), 
                                                    0, 0)
    # remove data when Temp is NaN - avoids biasing difference low
    pH_ESPER_wo_O2_no_nan = np.copy(pH_ESPER_wo_O2)
    pH_ESPER_wo_O2_no_nan[np.isnan(Measurements_all_S_T[:,1])] = np.nan

    # save 1500 pH orig
    argo_derived_n['pH_1500_ESPER_wo_O2'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
    argo_derived_n.pH_1500_ESPER_wo_O2[:] = pH_ESPER_wo_O2_no_nan[:,0]



    # create a second Measurements array that adjusts oxygen according to the mean glodap bias. Subtract the offset to correct it properly
    # create a measurements_offset array that only has S, O2, and T:
    Measurements_S_T_O2_w_o2_offset = np.concatenate(([Measurements_all_S_T_O2[:,0]], 
                                    [Measurements_all_S_T_O2[:,1]],
                                    [Measurements_all_S_T_O2[:,2]] - 
                                    glodap_offsets_p.DOXY_ADJUSTED_offset_trimmed[glodap_offsets_p.index==wmo_n].values))
    Measurements_S_T_O2_w_o2_offset = np.transpose(Measurements_S_T_O2_w_o2_offset)

    # calculate LIPHR pH at 1500 m with oxygen adjustment
    pH_o2_adjust = carbon_utils.LIPHR_matlab(LIPHR_path,
                                    Coordinates_all.tolist(),
                                    Measurements_S_T_O2_w_o2_offset.tolist(),
                                    MeasIDVec_LIR, 
                                    OAAdjustTF = False)  

    # remove data when oxygen is NaN - avoids biasing difference low
    pH_o2_adjust_no_nan = np.copy(pH_o2_adjust)
    pH_o2_adjust_no_nan[np.isnan(Measurements_all_S_T_O2[:,2])] = np.nan

    argo_derived_n['pH_1500_O2_ADJUST_LIPHR'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
    argo_derived_n.pH_1500_O2_ADJUST_LIPHR[:] = pH_o2_adjust_no_nan[:,0]

    # Impact of O2 can be seen in the test_pH minus new_pH average 
    argo_derived_n['PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST'] = \
        argo_derived_n.PH_IN_SITU_TOTAL_ADJUSTED + np.nanmean(pH_o2_adjust_no_nan - pH_orig_no_nan)* \
        ((np.nanmean(Measurements_all_S_T_O2[:,1]+273.15))/(argo_derived_n.TEMP_ADJUSTED + 273.15)) #2023_10_02 adding temperature dependency of k0
    

    # calculate LIPHR pH at 1500 m with oxygen adjustment
    pH_ESPER_o2_adjust = carbon_utils.ESPER_mixed_matlab(LIPHR_path,
                                                    DesiredVariables,
                                                    Coordinates_all.tolist(),
                                                    Measurements_S_T_O2_w_o2_offset.tolist(),
                                                    MeasIDVec_ESPER,
                                                    Equations, 
                                                    decimal_year.values.tolist(), 
                                                    0, 0)

    # remove data when oxygen is NaN - avoids biasing difference low
    pH_ESPER_o2_adjust_no_nan = np.copy(pH_ESPER_o2_adjust)
    pH_ESPER_o2_adjust_no_nan[np.isnan(Measurements_all_S_T_O2[:,2])] = np.nan

    argo_derived_n['pH_1500_O2_ADJUST_ESPER'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
    argo_derived_n.pH_1500_O2_ADJUST_ESPER[:] = pH_ESPER_o2_adjust_no_nan[:,0]

    # Impact of O2 can be seen in the test_pH minus new_pH average 
    argo_derived_n['PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST_ESPER'] = \
        argo_derived_n.PH_IN_SITU_TOTAL_ADJUSTED + np.nanmean(pH_ESPER_o2_adjust_no_nan - pH_orig_ESPER_no_nan)* \
        ((np.nanmean(Measurements_all_S_T_O2[:,1]+273.15))/(argo_derived_n.TEMP_ADJUSTED + 273.15)) #2023_10_02 adding temperature dependency of k0
    
    # Impact of removing O2 can be seen in the test_pH minus new_pH average 
    argo_derived_n['PH_IN_SITU_TOTAL_ADJUSTED_wo_O2_ESPER'] = \
        argo_derived_n.PH_IN_SITU_TOTAL_ADJUSTED + np.nanmean(pH_ESPER_wo_O2_no_nan - pH_orig_ESPER_no_nan)* \
        ((np.nanmean(Measurements_all_S_T[:,1]+273.15))/(argo_derived_n.TEMP_ADJUSTED + 273.15)) #2023_10_02 adding temperature dependency of k0
    

     # do additional nitrate calculations if NITRATE is present
    if 'NITRATE_ADJUSTED' in argo_derived_n.keys():
            
        # calculate LINR NO3 at 1500m using original data from all profiles
        orig_no3 = carbon_utils.LINR_matlab(LIPHR_path,
                                        Coordinates_all.tolist(),
                                        Measurements_all_S_T_O2.tolist(),
                                        MeasIDVec_LIR)  
        
        argo_derived_n['NO3_1500_orig_LINR'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
        argo_derived_n.NO3_1500_orig_LINR[:] = orig_no3[:,0]

        # calculate LINR NO3 at 1500m with corrected oxygen 
        no3_o2_adjust = carbon_utils.LINR_matlab(LIPHR_path,
                                        Coordinates_all.tolist(),
                                        Measurements_S_T_O2_w_o2_offset.tolist(),
                                        MeasIDVec_LIR)  
        # remove data when oxygen is NaN - avoids biasing difference low
        no3_o2_adjust_no_nan = np.copy(no3_o2_adjust)
        no3_o2_adjust_no_nan[np.isnan(Measurements_all_S_T_O2[:,2])] = np.nan

        argo_derived_n['NO3_1500_o2_adjust'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
        argo_derived_n.NO3_1500_o2_adjust[:] = no3_o2_adjust_no_nan[:,0]
        
        argo_derived_n['NITRATE_ADJUSTED_w_O2_ADJUST'] = \
            argo_derived_n.NITRATE_ADJUSTED + np.nanmean(no3_o2_adjust_no_nan - orig_no3)

        MeasIDVec_LIR_no_o2 = [1, 7]
        # calculate LINR NO3 at 1500m without using oxygen
        no3_no_o2 = carbon_utils.LINR_matlab(LIPHR_path,
                                        Coordinates_all.tolist(),
                                        Measurements_all_S_T.tolist(),
                                        MeasIDVec_LIR_no_o2)  
        # remove data when oxygen is NaN - avoids biasing difference low
        no3_no_o2_no_nan = np.copy(no3_no_o2)
        no3_no_o2_no_nan[np.isnan(Measurements_all_S_T[:,1])] = np.nan

        argo_derived_n['NO3_1500_no_o2'] = (['N_PROF'],np.empty(argo_derived_n.PRES_ADJUSTED.shape[0])) #nprof x nlevel
        argo_derived_n.NO3_1500_no_o2[:] = no3_no_o2_no_nan[:,0]
        
        argo_derived_n['NITRATE_ADJUSTED_wo_O2'] = \
            argo_derived_n.NITRATE_ADJUSTED + np.nanmean(no3_no_o2_no_nan - orig_no3)

    argo_derived_n.to_netcdf(argo_path_derived+ str(wmo_n) + '_derived_2.nc')
    print(str(wmo_n)+ ' finished')

def calc_derived_3_pCO2_impacts(argo_path_derived, file):
    print('Starting: ' + file) 
    argo_n = xr.load_dataset(argo_path_derived + file)
    var_list = list(argo_n.data_vars)
    wmo_n = argo_n['wmo']

    if "PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST" not in var_list:
            # print('Skipping')
            return

    if 'TALK_LIAR' not in argo_n.keys():
            return

    if 'NITRATE_ADJUSTED' in argo_n.keys():
            SI = argo_n.NITRATE_ADJUSTED*2.5
            SI.where(~np.isnan(SI), 0)
            PO4 = argo_n.NITRATE_ADJUSTED/16
            PO4.where(~np.isnan(PO4),0)
    else:
            SI = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
            PO4 = np.zeros((argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
   
    # Calculate pCO2 with all original values
    # replace any negative pH values with nans:
    argo_n['PH_IN_SITU_TOTAL_ADJUSTED'] = \
            argo_n.PH_IN_SITU_TOTAL_ADJUSTED.where(argo_n.PH_IN_SITU_TOTAL_ADJUSTED>0)

    results = pyco2.sys(
            par1=argo_n.TALK_LIAR, 
            par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
            par1_type=1,
            par2_type=3,
            temperature=argo_n.TEMP_ADJUSTED, 
            pressure=argo_n.PRES_ADJUSTED, 
            salinity=argo_n.PSAL_ADJUSTED, 
            temperature_out=argo_n.TEMP_ADJUSTED, #fixed 25C temperature
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
    argo_n['pCO2_pH_orig_TALK_orig'] = (['N_PROF','N_LEVELS'],results['pCO2'])  
    argo_n['DIC_pH_orig_TALK_orig'] = (['N_PROF','N_LEVELS'],results['dic'])  

    # Calculate pCO2 including only the impact of adjusting oxygen on pH MLR
    # replace any negative pH values with nans:
    argo_n['PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST'] = \
            argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST.where(argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST>0)

    results = pyco2.sys(
            par1=argo_n.TALK_LIAR, 
            par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST,
            par1_type=1,
            par2_type=3,
            temperature=argo_n.TEMP_ADJUSTED, 
            pressure=argo_n.PRES_ADJUSTED, 
            salinity=argo_n.PSAL_ADJUSTED, 
            temperature_out=argo_n.TEMP_ADJUSTED, #fixed 25C temperature
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
    argo_n['pCO2_pH_O2_TALK_orig'] = (['N_PROF','N_LEVELS'],results['pCO2'])  
    argo_n['DIC_pH_O2_TALK_orig'] = (['N_PROF','N_LEVELS'],results['dic'])  

    # ESPER Version: Calculate pCO2 including only the impact of adjusting oxygen on pH MLR
    # replace any negative pH values with nans:
    argo_n['PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST_ESPER'] = \
            argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST_ESPER.where(argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST_ESPER>0)

    results = pyco2.sys(
            par1=argo_n.TALK_LIAR, 
            par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED_w_O2_ADJUST_ESPER,
            par1_type=1,
            par2_type=3,
            temperature=argo_n.TEMP_ADJUSTED, 
            pressure=argo_n.PRES_ADJUSTED, 
            salinity=argo_n.PSAL_ADJUSTED, 
            temperature_out=argo_n.TEMP_ADJUSTED, #fixed 25C temperature
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
    argo_n['pCO2_pH_O2_ESPER_TALK_orig'] = (['N_PROF','N_LEVELS'],results['pCO2'])  
    argo_n['DIC_pH_O2_ESPER_TALK_orig'] = (['N_PROF','N_LEVELS'],results['dic'])  
 
    # ESPER Version with oxygen removed entirely
    argo_n['PH_IN_SITU_TOTAL_ADJUSTED_wo_O2_ESPER'] = \
            argo_n.PH_IN_SITU_TOTAL_ADJUSTED_wo_O2_ESPER.where(argo_n.PH_IN_SITU_TOTAL_ADJUSTED_wo_O2_ESPER>0)

    results = pyco2.sys(
            par1=argo_n.TALK_LIAR, 
            par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED_wo_O2_ESPER,
            par1_type=1,
            par2_type=3,
            temperature=argo_n.TEMP_ADJUSTED, 
            pressure=argo_n.PRES_ADJUSTED, 
            salinity=argo_n.PSAL_ADJUSTED, 
            temperature_out=argo_n.TEMP_ADJUSTED, #fixed 25C temperature
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
    argo_n['pCO2_pH_wo_O2_ESPER_TALK_orig'] = (['N_PROF','N_LEVELS'],results['pCO2'])  
    argo_n['DIC_pH_wo_O2_ESPER_TALK_orig'] = (['N_PROF','N_LEVELS'],results['dic'])  

    argo_n.to_netcdf(argo_path_derived+ str(wmo_n.values) + '_derived_3.nc')