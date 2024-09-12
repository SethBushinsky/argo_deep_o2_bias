import xarray as xr
import numpy as np
import functions.carbon_utils as carbon_utilities
import PyCO2SYS as pyco2
from scipy import interpolate
import glob, os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

def argo_interp_profiles(argo_path, LIAR_path, argo_path_interpolated, argo_path_derived, argo_file, qc_data_fields, bgc_data_fields, p_interp, \
                          derived_list, interpolation_list, adjustment):
    print('Processing float file '+ argo_file)
    argo_n = xr.load_dataset(argo_path+argo_file)
    argo_n = argo_n.set_coords(('PRES_ADJUSTED','LATITUDE','LONGITUDE','JULD'))

    wmo_n = argo_n.PLATFORM_NUMBER.values.astype(int)[0]
    # wmo_list.append(wmo_n)

    nprof_n = argo_n.dims['N_PROF']

    p_interp_min = p_interp[0]
    p_interp_max = p_interp[-1]

    #   set bad data and possibly bad data to NaN 
    for q in qc_data_fields:      
        if q in argo_n.keys():
            qc_val = argo_n[q+'_QC'].values.astype('float')
            
            # for some reason the .where statement was not filtering out bad values. 
            # This code is now changing QC values of 0 (no qc), 3(probably bad), 4(bad), and 9 (missing value) to nans. 
            # interpolated values are set to nan next for BGC data
            #argo_n[q].where(np.logical_and(qc_val<3.,qc_val>4.))
            argo_n[q].values[np.logical_or(qc_val==4,qc_val==3)]=np.nan
            argo_n[q].values[np.logical_or(qc_val==0,qc_val==9)]=np.nan

            #check for any Inf values not included in QC flag and set to NaN
            argo_n[q].values[np.isinf(argo_n[q]).values] = np.nan
        
    # check for interpolated profile positions (under ice) and set all BGC data to nan
    qc_val = argo_n['POSITION_QC'].values.astype('float')
    for b in bgc_data_fields:
        if b in argo_n.keys() and np.any(qc_val==8):
            naninds = np.argwhere(qc_val==8)[:,0]
            argo_n[b][naninds,:] = np.nan
    
    #Finding and removing all non-delayed mode data
        # sometimes parameters are missing from profiles - 
        # need to loop through all profiles and check which parameters are present
    parameter_array = argo_n.STATION_PARAMETERS.values.astype(str)

    for idx in range(len(parameter_array)):
        prof_parameters = parameter_array[idx]
        # print(prof_parameters)

        # loop through each paramter in the profile 
        for var in prof_parameters:
            var_str = var.strip()
            if len(var_str)==0: # only proceed if the variable exists 
                continue
            var_ind = [idx for idx, s in enumerate(prof_parameters) if s.strip()== var_str]
            # print(var_ind)

            # get parameter data mode values for that profile / variable
            var_data_mode = argo_n.PARAMETER_DATA_MODE[idx,var_ind].values
            # print(var_data_mode)

            decoded_arr = np.array([elem.decode() if isinstance(elem, bytes) else np.nan for elem in var_data_mode.flatten()])
            # print(decoded_arr)
            result = np.where(decoded_arr == 'D', False, True) # true whereever mode is not delayed
            # print(result)
            if result:
                argo_n[var_str +'_ADJUSTED'][idx,:] = np.nan

    # we are currently processing floats that have no valid biogeochemical data. 
    #Should check to see if data in key 
    #original bgc parameters (O2, NO3, pH) is valid and skip the rest if not
    bgc_valid = 0
    for b in bgc_data_fields:
        if b in argo_n.keys() and np.any(~np.isnan(argo_n[b])):
            bgc_valid = bgc_valid+1
    if bgc_valid >=1:
        print(argo_file + ' has valid BGC data')
    else:
        print(argo_file + ' has no valid BGC data')
        return
    
    argo_n['PDENS'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
    argo_n.PDENS[:] = np.nan
    argo_n['spice'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
    argo_n.spice[:] = np.nan

    #initialise interpolated dataset for float
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
    for v in derived_list: # all the variables that will be saved out in the derived and interpolated files
        argo_interp_n[v] = (['N_PROF','N_LEVELS'],np.copy(nan_interp))

    # if reading in adjustment / offset data, load impacts and apply as appropriate
    if adjustment is True:
        impact_n = xr.load_dataset(argo_path + argo_file[0:7] + '_impact.nc')
        argo_n['DOXY_ADJUSTED'] = argo_n['DOXY_ADJUSTED'] - impact_n.mean_O2_offset
        if 'NITRATE_ADJUSTED' in argo_n.keys():
            argo_n['NITRATE_ADJUSTED'] = argo_n['NITRATE_ADJUSTED'] + impact_n.mean_nitrate_impact_change

    #check first if PH_IN_SITU_TOTAL_ADJUSTED exists
    if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
        
        print('Calculating TALK, DIC and pH 25C correction for float '+str(wmo_n))
        
        #initialise pH 25c and DIC variables - could do this only if float has pH
        argo_n['TALK_LIAR'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.TALK_LIAR[:] = np.nan
        argo_n['pH_25C_TOTAL_ADJUSTED'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.pH_25C_TOTAL_ADJUSTED[:] = np.nan
        argo_n['DIC'] = (['N_PROF','N_LEVELS'],np.empty(argo_n.PRES_ADJUSTED.shape)) #nprof x nlevel
        argo_n.DIC[:] = np.nan

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


        results = carbon_utilities.LIAR_matlab(LIAR_path,
                                                Coordinates.tolist(),
                                                Measurements.tolist(),
                                                MeasIDVec,
                                                VerboseTF=False)                                  

        argo_n['TALK_LIAR'] = (['N_PROF','N_LEVELS'],
                                np.reshape(np.asarray(results),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))


        # Keep DIC bc I might want it for crossover comparison
        ##### Calculate float pH at 25C, DIC and apply bias corr
        results = pyco2.sys(
                par1=argo_n.TALK_LIAR, 
                par2=argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
                par1_type=1,
                par2_type=3,
                temperature=argo_n.TEMP_ADJUSTED, 
                pressure=argo_n.PRES_ADJUSTED, 
                salinity=argo_n.PSAL_ADJUSTED, 
                temperature_out=25.,#*np.ones(argo_n.PRES_ADJUSTED.shape), #fixed 25C temperature
                pressure_out=1500., #argo_n.PRES_ADJUSTED, # fixed 1500 db output pressure
                total_silicate=SI,
                total_phosphate=PO4,
                opt_pH_scale = 1, #total
                opt_k_carbonic=10, #Lueker et al. 2000
                opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
                opt_total_borate=2, # Lee et al. 2010
                opt_k_fluoride=2, # Perez and Fraga 1987
                opt_buffers_mode=1,
        )
        argo_n['pH_25C_T_P1500'] = (['N_PROF','N_LEVELS'], results['pH_total_out'])
        argo_n['pH_25C_TOTAL_ADJUSTED'] = (['N_PROF','N_LEVELS'],carbon_utilities.co2sys_pH25C(argo_n.TALK_LIAR,
                                                    argo_n.PH_IN_SITU_TOTAL_ADJUSTED,
                                                    argo_n.TEMP_ADJUSTED,
                                                    argo_n.PSAL_ADJUSTED,
                                                    argo_n.PRES_ADJUSTED))

        # if applying adjustment to pH - apply it to pH 25C, then recalculate pH insitu, then calculate DIC
        if adjustment is True:
            argo_n['pH_25C_TOTAL_ADJUSTED'] = argo_n['pH_25C_TOTAL_ADJUSTED'] + impact_n.mean_pH_impact_change # note that I am not correcting the in situ pH

            results = pyco2.sys(
                par1=argo_n.TALK_LIAR, 
                par2=argo_n.pH_25C_TOTAL_ADJUSTED, # using the impact adjusted pH
                par1_type=1,
                par2_type=3,
                temperature=25, 
                pressure=argo_n.PRES_ADJUSTED, 
                salinity=argo_n.PSAL_ADJUSTED, 
                temperature_out=argo_n.TEMP_ADJUSTED,#*np.ones(argo_n.PRES_ADJUSTED.shape), #fixed 25C temperature
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
            argo_n['DIC'] = (['N_PROF','N_LEVELS'],results['dic'])  
        else: # otherwise, just save DIC with no adjustment 
            argo_n['DIC'] = (['N_PROF','N_LEVELS'],results['dic'])  
            
    ##### now calc potential density, save, and interpolate data for comparison
    for p in range(nprof_n):
        #pressure for profile
        p_prof = argo_n.PRES_ADJUSTED[p,:]
        
        # For interpolated data, shouldn't calculate pdens and spice and then interpolate - 
        # should interpolate psal and temp and then calculate spice and pdens
        # Do both so that you are able to have PDENS and spice in the derived files too (do I need them?)
        argo_n['PDENS'][p,:] = carbon_utilities.sigma0(argo_n.PSAL_ADJUSTED[p,:].values,
                                                    argo_n.TEMP_ADJUSTED[p,:].values,
                                                    argo_n.LONGITUDE[p].values,
                                                    argo_n.LATITUDE[p].values,
                                                    argo_n.PRES_ADJUSTED[p,:].values)
        argo_n['spice'][p,:] = carbon_utilities.spiciness0(argo_n.PSAL_ADJUSTED[p,:].values,
                                                    argo_n.TEMP_ADJUSTED[p,:].values,
                                                    argo_n.LONGITUDE[p].values,
                                                    argo_n.LATITUDE[p].values,
                                                    argo_n.PRES_ADJUSTED[p,:].values)

        #for each profile get pressure values > p_interp_min db
        p100 = p_prof[p_prof>p_interp_min].values
            
        #if only 1 value of pressure or if there is not valid profile data down to p-max, continue loop
        if (len(p100) <= 1) or (np.nanmax(p100)<p_interp_min):
            continue
        
        # # check for the presence of large gaps in the float profile data - can figure out how to deal with them once you know their prevalence 
        # if max(np.diff(p100))>125:
        #     print(np.diff(p100))
        #     data_out = p100.reshape(-1,1)
        #     df = pd.DataFrame(data_out, columns = ['Pressure prior to interpolation'])
        #     df.to_csv(argo_path_interpolated + str(wmo_n) + '_' + str(p) + '.csv', index=False)

        #find which crossover variables exist in main float file
        var_list_n = []
        for vname in interpolation_list:
            if (vname in argo_n.keys()) and (np.any(~np.isnan(argo_n[vname]))):
                var_list_n.append(vname)
                
        argo_interp_n['num_var'][p] = len(var_list_n) 
        
        for var in var_list_n:
            var100 = argo_n[var][p,p_prof>p_interp_min]

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
                (np.nanmax(p100u[~np.isnan(var100u.values)])>p_interp_min) and \
                (np.nanmin(p100u[~np.isnan(var100u.values)])<p_interp_max):
                
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
            
                # check for gaps in the original data greater than 125 m
                gap_index= (np.diff(p100u)>125)

                # if any values of gap_index are true, loop through and set values of interpolated data that are between value of large gaps to nan 
                if any(gap_index):
                    # temp_var = argo_interp_n[var][p,:]
                    # data_out = temp_var.values.reshape(-1,1)
                    # combined_data = np.hstack((p_interp.reshape(-1, 1), data_out))
                    # df = pd.DataFrame(combined_data, columns = ['Pressure', 'Oxygen'])
                    # df.to_csv(argo_path_interpolated + str(wmo_n) + '_' + str(p) + var + '.csv', index=False)

                    # print(argo_interp_n[var][p,:])
                    for idx, gi in enumerate(gap_index):
                        if gi:
                            # print(p100u[idx])
                            # print(p100u[idx+1])
                            argo_interp_n[var][p,np.logical_and(p_interp>p100u[idx],p_interp<p100u[idx+1])] = np.nan
                    # temp_var = argo_interp_n[var][p,:]
                    # data_out = temp_var.values.reshape(-1,1)
                    # combined_data = np.hstack((p_interp.reshape(-1, 1), data_out))
                    # df = pd.DataFrame(combined_data, columns = ['Pressure', 'Oxygen'])
                    # df.to_csv(argo_path_interpolated + str(wmo_n) + '_' + str(p) + var + '_after_removal.csv', index=False)

    #             else: 
                # print('profile data not deep enough to interpolate ' + str(p) + ' ' +  var)
                #                       str(np.nanmax(p100u[~np.isnan(var100u.values)])))
                # print('values greater than min ' + str(var100u[p100u>p_interp_min].values))

    # loop through profiles again to calculate PDENS and spice for interpolated dataset
    for p in range(nprof_n):
        #pressure for profile
        p_prof = argo_interp_n.PRES_ADJUSTED[p,:]

        # For interpolated data, shouldn't calculate pdens and spice and then interpolate - 
        # should interpolate psal and temp and then calculate spice and pdens
        # Do both so that you are able to have PDENS and spice in the derived files too (do I need them?)
        argo_interp_n['PDENS'][p,:] = carbon_utilities.sigma0(argo_interp_n.PSAL_ADJUSTED[p,:].values,
                                                    argo_interp_n.TEMP_ADJUSTED[p,:].values,
                                                    argo_interp_n.LONGITUDE[p].values,
                                                    argo_interp_n.LATITUDE[p].values,
                                                    argo_interp_n.PRES_ADJUSTED[p,:].values)
        argo_interp_n['spice'][p,:] = carbon_utilities.spiciness0(argo_interp_n.PSAL_ADJUSTED[p,:].values,
                                                    argo_interp_n.TEMP_ADJUSTED[p,:].values,
                                                    argo_interp_n.LONGITUDE[p].values,
                                                        argo_interp_n.LATITUDE[p].values,
                                                        argo_interp_n.PRES_ADJUSTED[p,:].values)
        
    #create new dataset with relevant crossover variables only
    argo_n_derived = xr.Dataset()
    argo_n_derived['wmo'] = wmo_n
    argo_n_derived['CYCLE_NUMBER'] = (['N_PROF'],argo_n.CYCLE_NUMBER.values)
    argo_n_derived['LONGITUDE'] = (['N_PROF'],argo_n.LONGITUDE.values)
    argo_n_derived['LATITUDE'] = (['N_PROF'],argo_n.LATITUDE.values)
    argo_n_derived['JULD_LOCATION'] = (['N_PROF'],argo_n.JULD_LOCATION.values)
    for var in derived_list:
        if var in argo_n.keys():
            argo_n_derived[var] = (['N_PROF','N_LEVELS'],argo_n[var].values)
    argo_n_derived.to_netcdf(argo_path_derived+str(wmo_n)+'_derived.nc')

    argo_interp_n.to_netcdf(argo_path_interpolated+str(wmo_n)+'_interpolated.nc')

def glodap_crossover_offsets(argo_path_interpolated, offset_dir, glodap_file_offsets_dir, argo_file, dist, delta_dens, \
                                delta_spice, delta_press, gdap_p, p_interp, plot_profile, var_list_plot, \
                                p_compare_min, p_compare_max):
    
    wmo = argo_file[0:7]
    print('Starting crossover for '+ wmo)

    try:
        argo_interp_n = xr.open_dataset(argo_path_interpolated + wmo + '_interpolated.nc')
    except:
        print('File ' + wmo + '_interpolated.nc' + ' not found:')
        return
    #number of profiles
    nprof = argo_interp_n.LATITUDE.shape[0]

    p_interp_min = p_interp[0]
    p_interp_max = p_interp[-1]
    #initiate offset list
    #number of additional rows for saving metadata items
    num_meta_items = 11
    gdap_offsets =  [[] for _ in range(3*len(var_list_plot)+num_meta_items)]
 
     #check sufficient non-NaN data
    if np.sum(argo_interp_n.num_var>3)==0: #changed to look at all profiles and to only select floats 
        # with more than 3 variables present (i.e. more than T, S, Press)
        print('No non-NAN bgc adjusted data for: '+ wmo)
        return
    
    #  Find lat-lon limits within dist of float location
    lat_tol = dist/ 111.6 
    lon_tol = dist/ (111.6*np.cos(np.deg2rad(np.nanmean(argo_interp_n.LATITUDE))))  
    #set lat/lon crossover limits
    lat_min = argo_interp_n.LATITUDE.values-lat_tol
    lat_max = argo_interp_n.LATITUDE.values+lat_tol
    lon_min = argo_interp_n.LONGITUDE.values-lon_tol
    lon_max = argo_interp_n.LONGITUDE.values+lon_tol

    #find all data in lat-lon limits
    for n in range(nprof): #range(4,5):# range(3, 4): #   #
        #print(argo_interp_n.profile[n].values)

        #index of all gdap profiles within distance range
        if lon_min[n] < 0 and lon_max[n]>0:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n]+360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]))
            
        elif lon_min[n] < 0 and lon_max[n]<0:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_and(gdap_p.LONGITUDE.values>lon_min[n]+360,
                                                  gdap_p.LONGITUDE.values<lon_max[n]+360))  
            
        elif lon_max[n] > 360 and lon_min[n] < 360:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_or(gdap_p.LONGITUDE.values>lon_min[n],
                                                  gdap_p.LONGITUDE.values<lon_max[n]-360))

        elif lon_max[n] > 360 and lon_min[n] > 360:
            match = np.logical_and(np.logical_and(gdap_p.LATITUDE.values>lat_min[n],
                                                  gdap_p.LATITUDE.values<lat_max[n]),
                                   np.logical_and(gdap_p.LONGITUDE.values>lon_min[n]-360,
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
        dens_inds = np.empty(match_inds.shape[0],dtype=int)
        o2_offsets = np.empty(match_inds.shape[0])
        o2_offsets[:] = np.nan
        dens_offsets = np.empty(match_inds.shape[0])
        o2_offsets[:] = np.nan
        

        for m in range(len(match_inds)):
        
            #if no interpolated density (not deep enough) then move on to next profile
            # probably should move the float profile check earlier in this cell, no need to run so much extra code 
            if np.all(np.isnan(argo_interp_n.PDENS.values[n,:])) or np.isnan(gdap_match.PDENS.values[m]):
#                 print('no interpolated density in depth range')
                dens_inds[m]=-1
                continue
                
            #calculate density difference and nearest density index
            dens_diff = np.absolute(argo_interp_n.PDENS.values[n,:]-gdap_match.PDENS.values[m])
            dens_ind = np.nanargmin(dens_diff)

            dens_inds[m] = dens_ind #save indices for plotting 

            
            #don't use matches if the density difference to the closest matching density value exceeds delta_dens
            if dens_diff[dens_ind] > delta_dens:
                continue
                
            #don't use matches at top/bottom of interpolated vector
            if (dens_ind == 0) or (dens_ind == len(p_interp) -1):
                continue
                
            ##optional: add spice criteria
            spice_diff = np.absolute(argo_interp_n.spice.values[n,dens_ind]-gdap_match.spice.values[m])
            
            #don't use matches if spice difference exceeds delta_spice
            if spice_diff > delta_spice:
                continue
                
            ## add depth criteria
            press_diff = np.absolute(argo_interp_n.PRES_ADJUSTED.values[n,dens_ind]-gdap_match.PRES_ADJUSTED.values[m])
            if press_diff> delta_press:
                continue
            #calc offset at the interp profile index for each match for each var to calculate crossover and or plot 
            #check var is present in both float and glodap
            # modifying for parallel - save individual gdap_offsets lists
            for idx, var in enumerate(var_list_plot):
                if var in argo_interp_n.keys():
                    gdap_offset_ind = gdap_match[var].index[m]
                    offset = argo_interp_n[var][n,dens_ind] - gdap_match[var][gdap_offset_ind]
                    
                    #save offsets for optional plotting
                    if var=='DOXY_ADJUSTED':
                        o2_offsets[m] = offset.values
                    elif var=='PDENS':
                        dens_offsets[m] = offset.values
                        
                    #append to full offset list
                    gdap_offsets[idx*3].append(offset.values)
                    
                    #also save absolute glodap value at crossover
                    gdap_offsets[idx*3+1].append(gdap_match[var][gdap_offset_ind])
                    #save absolute float value at crossover
                    gdap_offsets[idx*3+2].append(argo_interp_n[var][n,dens_ind])
#                     print(var)
                #append nan if variable is not there so lists all remain same length?
                else:
                    gdap_offsets[idx*3].append(np.nan) 
                    gdap_offsets[idx*3+1].append(np.nan) 
                    gdap_offsets[idx*3+2].append(np.nan) 
            
            #append metadata to offset list
            gdap_offsets[len(var_list_plot)*3].append(wmo)
            gdap_offsets[len(var_list_plot)*3+1].append(argo_interp_n.profile[n].values)
            gdap_offsets[len(var_list_plot)*3+2].append(argo_interp_n.juld[n].values)
            gdap_offsets[len(var_list_plot)*3+3].append(gdap_match.datetime.values[m])
            gdap_offsets[len(var_list_plot)*3+4].append(argo_interp_n.LONGITUDE[n].values)
            gdap_offsets[len(var_list_plot)*3+5].append(gdap_match.LONGITUDE.values[m])
            gdap_offsets[len(var_list_plot)*3+6].append(argo_interp_n.LATITUDE[n].values)
            gdap_offsets[len(var_list_plot)*3+7].append(gdap_match.LATITUDE.values[m])
            gdap_offsets[len(var_list_plot)*3+8].append(gdap_match.G2cruise.values[m])
            gdap_offsets[len(var_list_plot)*3+9].append(gdap_match.G2station.values[m])
            gdap_offsets[len(var_list_plot)*3+10].append(gdap_match.obs_index.values[m])

        
            #can add additional float metadata variable to list here
     

        ################################
        #optional: plot individual profiles and offsets to check density match up
        if plot_profile == 1:
            # if there are no non-nan values of o2_offsets for this profile, move on
            if np.all(np.isnan(o2_offsets)):
                continue
            
            #if no interpolated density (not deep enough) then move on to next profile
            if np.all(np.isnan(argo_interp_n.PDENS.values[n,:])) or np.isnan(gdap_match.PDENS.values[m]):
                continue
                
            # create a folder for profile plots if one does not exist:
            if not os.path.isdir(offset_dir + 'individual_floats'):
                os.mkdir(offset_dir + 'individual_floats')
                
            if not os.path.isdir(offset_dir + 'individual_floats/' + str(wmo)):
                os.mkdir(offset_dir + 'individual_floats/' + str(wmo))
    
            print('Plotting float '+str(argo_interp_n.wmo[0].values) + ' profile' + str(argo_interp_n.profile[n].values))
            fig = plt.figure(figsize=(10,16))
            
            #DOXY_ADJUSTED
            #dens_inds = dens_inds[dens_inds>0] #only keep valid indices
            
            max_o2 = 0
            min_o2 = 300
            max_dens = 0
            min_dens = 300
            min_pres = 1600
            max_pres = 1700
            
            axn = plt.subplot(3,3,1)
            #plot float profile and matchups
            axn.scatter(argo_interp_n.DOXY_ADJUSTED[n,:],argo_interp_n.PRES_ADJUSTED[n,:],s=4,c='orangered',marker='o')
            # commenting this out, because dens_inds is specific to an individual glodap crossover   
            # if argo_interp_n.DOXY_ADJUSTED[n, dens_inds].size>0:
            #     axn.scatter(argo_interp_n.DOXY_ADJUSTED[n,dens_inds],argo_interp_n.PRES_ADJUSTED[n,dens_inds],s=50,c='k',marker='s')
            #plot glodap profiles and matchups- (split glodap into profiles by cruise and station?)

            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match.DOXY_ADJUSTED.index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    #print(m)
                    axn.scatter(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]],argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                    #plot all glodap data that match these cruise and stations
                    axn.plot(gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'DOXY_ADJUSTED'].values,
                         gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'PRES_ADJUSTED'].values,
                        'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.DOXY_ADJUSTED.values[m],
                            gdap_match.PRES_ADJUSTED.values[m],s=50,c='b',marker='D')
                    
                                    
                    if np.nanmin(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]])<min_o2:
                        min_o2 = np.nanmin(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]])
                    
                    if np.nanmin(gdap_match.DOXY_ADJUSTED.values[m])<min_o2:
                        min_o2 = np.nanmin(gdap_match.DOXY_ADJUSTED.values[m])

                    if np.nanmax(gdap_match.DOXY_ADJUSTED.values[m])>max_o2:
                        max_o2 = np.nanmax(gdap_match.DOXY_ADJUSTED.values[m])
                    
                    if np.nanmax(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]])>max_o2:
                        max_o2 = np.nanmax(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]])
                                        
                    if np.nanmin(argo_interp_n.PDENS[n,dens_inds[m]])<min_dens:
                        min_dens = np.nanmin(argo_interp_n.PDENS[n,dens_inds[m]])
                    
                    if np.nanmax(argo_interp_n.PDENS[n,dens_inds[m]])>max_dens:
                        max_dens = np.nanmax(argo_interp_n.PDENS[n,dens_inds[m]])
                        
                    if np.nanmin(argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]])<min_pres:
                        min_pres = np.nanmin(argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]])
                    
                    if np.nanmax(argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]])>max_pres:
                        max_pres = np.nanmax(argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]])

            plt.gca().invert_yaxis()
            plt.ylim([2000,1200])
            plt.ylabel('pressure')
            plt.xlabel('DOXY_ADJUSTED')
            plt.title(str(argo_interp_n.wmo[0].values))
            plt.xlim([min_o2 - 20, 
                      max_o2 + 20])            
             
             #PDENS
            axn = plt.subplot(3,3,2)
            #plot float profile and matchups
            axn.scatter(argo_interp_n.PDENS[n,:],argo_interp_n.PRES_ADJUSTED[n,:],s=4,c='orangered',marker='o')
            
            #plot glodap profiles and matchups- (split glodap into profiles by cruise and station?)
            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match[0].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
    
                    axn.scatter(argo_interp_n.PDENS[n,dens_inds[m]],argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                                #plot all glodap data that match these cruise and stations
                    axn.plot(gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'PDENS'].values,
                         gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'PRES_ADJUSTED'].values,
                        'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.PDENS.values[m],
                            gdap_match.PRES_ADJUSTED.values[m],s=50,c='b',marker='D')
            plt.ylim([min_pres- 100, max_pres+100])

            plt.gca().invert_yaxis()
#             plt.ylim([2000,1200])
#             plt.xlim([27.5,28.0])
            plt.ylabel('pressure')
            plt.xlabel('PDENS')

            plt.xlim([min_dens-.01, max_dens+.01])

            
           #DOXY vs PDENS
            axn = plt.subplot(3,3,4)
            #plot float profile and matchups
            axn.scatter(argo_interp_n.DOXY_ADJUSTED[n,:],argo_interp_n.PDENS[n,:],s=4,c='orangered',marker='o')

            
            for m in range(len(match_inds)):
                gdap_offset_ind = gdap_match[0].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    axn.scatter(argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]],argo_interp_n.PDENS[n,dens_inds[m]],s=50,c='k',marker='s')

                    cruise = gdap_match.G2cruise[gdap_offset_ind]
                    stat = gdap_match.G2station[gdap_offset_ind]
                    #plot all glodap data that match these cruise and stations
                    axn.plot(gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'DOXY_ADJUSTED'].values,
                             gdap_p.loc[(gdap_p['G2cruise']==cruise) & (gdap_p['G2station']==stat),'PDENS'].values,
                            'm+-',linewidth=0.3)
                    axn.scatter(gdap_match.DOXY_ADJUSTED.values[m],
                                gdap_match.PDENS.values[m],s=50,c='b',marker='D')
    
            plt.gca().invert_yaxis()
            plt.ylabel('PDENS')
            plt.xlabel('DOXY_ADJUSTED')
            plt.xlim([min_o2 - 10, 
                      max_o2 + 10])
            plt.ylim([min_dens-.05, max_dens+.05])

            #print cruise names and dates
            axn = plt.subplot(3,3,3)
            axn.text(0.05,0.95,'G2cruise G2station G2date',fontsize=12)
            yn=0.9
            for m in range(len(match_inds)):
                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover

                    off_ind = gdap_match['G2cruise'].index[m]
                
                    #if m == 0:
                    axn.text(0.05,yn,str(gdap_match.G2cruise[off_ind])+' '
                             +str(gdap_match.G2station[off_ind])+' '
                             +str(gdap_match.datetime[off_ind])+' '
                             +str(round(o2_offsets[m],2)))
                    print()
                    yn=yn-0.05
                    #elif gdap_match.G2cruise[off_ind] != cc:
                    #    axn.text(0.05,yn,str(gdap_match.G2cruise[off_ind])+' '+str(gdap_match.G2station[off_ind])+' '+str(gdap_match.datetime[off_ind]))
                    #    yn=yn-0.05    
                    #cc = gdap_match.G2cruise[off_ind] #current cruise to compare with next 
            
            # map of offsets:
         #   ax = plt.subplot(3,3,5,subplot_kw={'projection': ccrs.PlateCarree()})
            ax = plt.subplot(3,3,5,projection=ccrs.PlateCarree())

            #ax = plt.axes(projection=ccrs.PlateCarree())
            ax.coastlines()
            ax.set_global()
            # find any longitude values over 180, then subtract
            temp_lon = gdap_match.LONGITUDE.values
            temp_lon[temp_lon>180] = temp_lon[temp_lon>180]-360

            max_lon = np.nanmax(temp_lon[~np.isnan(o2_offsets)])
            min_lon = np.nanmin(temp_lon[~np.isnan(o2_offsets)])
            max_lat = np.nanmax(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])
            min_lat = np.nanmin(gdap_match.LATITUDE.values[~np.isnan(o2_offsets)])

            # only plot gdap locations where offsets exist
            #plt.plot(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)],marker='x')
            scat = plt.scatter(temp_lon[~np.isnan(o2_offsets)],gdap_match.LATITUDE.values[~np.isnan(o2_offsets)], c=o2_offsets[~np.isnan(o2_offsets)])
            plt.xlim([min_lon - .05, max_lon+.05])
            plt.ylim([min_lat - .05, max_lat+.05])
            plt.colorbar(scat, ax=ax)
            #plt.show()
            
            # plot concentration vs. offset
            axn = plt.subplot(3,3,6)
#             plt.plot(o2_offsets, argo_interp_n.DOXY_ADJUSTED[n,:], 'bx')
            
            for m in range(len(match_inds)):
#                 gdap_offset_ind = gdap_match[var].index[m]

                if ~np.isnan(o2_offsets[m]): # only plot glodap data if there is a non-nan crossover
                    scat = plt.scatter(o2_offsets[m], argo_interp_n.DOXY_ADJUSTED[n,dens_inds[m]],s=50,
                                c=argo_interp_n.PRES_ADJUSTED[n,dens_inds[m]], marker='s')


    
            plt.colorbar()

            
            #histogram of offsets
            axn = plt.subplot(3,3,7)
            o2_offsets[o2_offsets==np.inf] = np.nan
            if not np.all(np.isnan(o2_offsets)):
                axn.hist(o2_offsets[~np.isnan(o2_offsets)],color='b')
                plt.xlabel('DOXY_OFFSETS')
                
            axn = plt.subplot(3,3,8)
            dens_offsets[dens_offsets==np.inf] = np.nan
            if not np.all(np.isnan(dens_offsets)):
                axn.hist(dens_offsets[~np.isnan(o2_offsets)],color='b')
                plt.xlabel('DENS_OFFSETS')
                plt.xticks(rotation = 45)
            plt.tight_layout()

            plt.savefig(offset_dir+ 'individual_floats/' + str(wmo) + '/' + str(argo_interp_n.wmo[0].values) + ' profile' + str(int(argo_interp_n.profile[n].values)))
            plt.clf()
        
    #convert GLODAP offset lists to xarray and save to netcdf
    glodap_offsets = xr.Dataset(coords={'N_CROSSOVERS':(['N_CROSSOVERS'],np.arange(len(gdap_offsets[0])))})

    glodap_offsets['p_compare_min'] = p_compare_min
    glodap_offsets['p_compare_max'] = p_compare_max
    glodap_offsets['delta_dens'] = delta_dens
    glodap_offsets['delta_spice'] = delta_spice
    glodap_offsets['delta_press'] = delta_press
    glodap_offsets['dist'] = dist

    for idx, var in enumerate(var_list_plot):
        glodap_offsets[var+'_offset'] = (['N_CROSSOVERS'], gdap_offsets[idx*3])
        glodap_offsets[var+'_glodap'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+1])
        glodap_offsets[var+'_float'] = (['N_CROSSOVERS'],gdap_offsets[idx*3+2])
        
    glodap_offsets['main_float_wmo'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3])
    glodap_offsets['main_float_profile'] = (['N_CROSSOVERS'], gdap_offsets[len(var_list_plot)*3+1])
    glodap_offsets['main_float_juld'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+2])
    glodap_offsets['glodap_datetime'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+3])
    glodap_offsets['main_float_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+4])
    glodap_offsets['glodap_longitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+5])
    glodap_offsets['main_float_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+6])
    glodap_offsets['glodap_latitude'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+7])
    glodap_offsets['glodap_cruise'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+8])
    glodap_offsets['glodap_station'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+9])
    glodap_offsets['glodap_obs_index'] = (['N_CROSSOVERS'],gdap_offsets[len(var_list_plot)*3+10])

    glodap_offset_filename = wmo + '_glodap_offset.nc'
    # print(glodap_offsets)
    # delete file if it exists already in case 
    if os.path.exists(glodap_file_offsets_dir+glodap_offset_filename):
        os.remove(glodap_file_offsets_dir+glodap_offset_filename)
    glodap_offsets.to_netcdf(glodap_file_offsets_dir+glodap_offset_filename)
    glodap_offsets.close()

