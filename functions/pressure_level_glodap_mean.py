import xarray as xr
import pandas as pd
import functions.outlier_filter_ESD_test as outlier
import numpy as np
import time

#list of DOXY_ADJUSTED SCIENTIFIC_CALIB_COMMENT substrings to group together for bias analysis

#no calibration bad data
bad_cal_list = ['Sensor issue','out of order','Bad data; not adjustable','Biofouling','unadjustable']

#no calibration, reason unspecified
no_cal_list = ['no adjustment','No QC','none','not applicable']

#blank cal
blank_cal = ['                ']

#air cal
# air_cal_list = ['DOXY_ADJUSTED corrected using co','SVU Foil','Partial pressure','Bittig',
#                 'Adjusted with SAGEO2 using co','Adjusted with SAGEO2 with in-air',
#                 'PPOX converted from DOXY','G determined from float measure']
air_cal_list = ['in air', 'in-air']
#no air cal
noair_cal_surf_list = ['World Ocean Atlas', 'woa', 'WOA', 'no in-air', 'no in air', 
                       'climatology','DOXY_QCs are modified during visual check',
                      'Takeshita']
# ['DOXY_ADJUSTED is computed from','DOXY_ADJUSTED is estimated from',
#                        'Adjustment done on PPOX_DOXY;Tem','Polynomial calibration','Adjustment via comparison of',
#                       'optode multi calibration','RT adjustment by comparison','adjusted with median',
#                       'Adjusted with annual WOA','Adjusted on WOA monthly','Adjusted PPOX at zero for',
#                       'G determined by surface']

noair_cal_subsurf_list = []
# ['1-point multiplicative corr']

noair_cal_funcofdoxy_list = []
# ['Percent saturation corrected as','DOXY_ADJUSTED computed using Ste',
#                         'Oxygen concentration corrected']

noair_cal_unspec_list = []
# ['DOXY_ADJUSTED corrected based','Adjust with WOA monthly','GAIN determined from WOA2013',
#                         'Adjusted with WOA climatology','Adjusted with SAGEO2 based on WO',
#                          'Adjusted with SAGEO2 on WOA','Adjusted with WOA 2018','Takeshita and all, 2013']

noair_cal_withdrift_list = []
# ['Adjustment done on PPOX_DOXY;Tem'] #this is incomplete

noair_cal_combined_list = noair_cal_surf_list+noair_cal_subsurf_list+noair_cal_funcofdoxy_list+noair_cal_unspec_list

def pressure_level_filter(argo_path, output_dir, out_filename, gdap_offsets_n_temp, var_list_plot, year_filt, 
                          pressure_level_min, pressure_level_max, 
                          year_plus_minus, time_delay, float_age_bins):
    time.sleep(time_delay) 
    # set all data not at that pressure level to nan                  
    for var in var_list_plot:
        pressure_index = np.logical_and(gdap_offsets_n_temp['PRES_ADJUSTED_float']>pressure_level_min, 
                                gdap_offsets_n_temp['PRES_ADJUSTED_float']<=pressure_level_max)

        if year_filt==0:
            gdap_offsets_n_temp[var + '_float'] = gdap_offsets_n_temp[var + '_float'].where(pressure_index)
            gdap_offsets_n_temp[var + '_glodap'] = gdap_offsets_n_temp[var + '_glodap'].where(pressure_index)
            gdap_offsets_n_temp[var + '_offset'] = gdap_offsets_n_temp[var + '_offset'].where(pressure_index)
        elif year_filt==1:
            year_index = np.abs(gdap_offsets_n_temp.main_float_juld.dt.year-gdap_offsets_n_temp.glodap_datetime.dt.year)<=year_plus_minus
            gdap_offsets_n_temp[var + '_float'] = \
                gdap_offsets_n_temp[var + '_float'].where(np.logical_and(pressure_index, year_index))
            gdap_offsets_n_temp[var + '_glodap'] = \
                gdap_offsets_n_temp[var + '_glodap'].where(np.logical_and(pressure_index, year_index))
            gdap_offsets_n_temp[var + '_offset'] = \
                gdap_offsets_n_temp[var + '_offset'].where(np.logical_and(pressure_index, year_index))
            
    # then group by wmo and proceed with DOXY_trimmed calculations
    offsets_g = gdap_offsets_n_temp.groupby(gdap_offsets_n_temp.main_float_wmo)

    if pressure_level_min<=400:
        # adding option to filter by time of year as well - for use in surface data
        time_filt = 1
        filt_days = 10
    else:
        time_filt=0


    # DOXY_ADJUSTED_offset_trimmed = []
    # DOXY_ADJUSTED_offset_trimmed = []
    # Create an empty data variable for DOXY offset trimmed with the same dimensions as N_CROSSOVERS
    empty_data = np.empty(len(gdap_offsets_n_temp['N_CROSSOVERS']))
    empty_data[:] = np.nan


    # Create a new xarray DataArray with the empty data and the same coordinates
    new_data_array = xr.DataArray(empty_data, coords={'N_CROSSOVERS': gdap_offsets_n_temp['N_CROSSOVERS']}, dims=['N_CROSSOVERS'])
    gdap_offsets_n_temp['DOXY_ADJUSTED_offset_trimmed'] = new_data_array

    for n,g in offsets_g:

        # run a GESD test using "test_num" number of possible outliers
        test_num = int((len(g.DOXY_ADJUSTED_offset.dropna(dim="N_CROSSOVERS", how="any").values)*.1)) # allowing for ~10 % to be outliers
        ESD_test_out = outlier.ESD_Test(g.DOXY_ADJUSTED_offset.dropna(dim="N_CROSSOVERS", how="any").values, 0.05, test_num, False, True)

        # only trim the data if deep, otherwise apply a day of year test but no other filtering 
        if time_filt==1:
            within_days = np.logical_or(np.abs(g.main_float_juld.dt.dayofyear - g.glodap_datetime.dt.dayofyear)<=filt_days, 
                            np.abs(g.main_float_juld.dt.dayofyear - g.glodap_datetime.dt.dayofyear)>(365-filt_days)) 
            temp_o2_offset = g.DOXY_ADJUSTED_offset.where(within_days)
        else:         # create temp_o2_offest to set all datapoints to nans that the GESD test says are outliers

            temp_o2_offset = g.DOXY_ADJUSTED_offset
            for a in range(0, ESD_test_out[1]):
                temp_o2_offset = temp_o2_offset.where(temp_o2_offset != ESD_test_out[2][a])

        # if there are too few points, set all to nans
        if temp_o2_offset.count()<20:
            temp_o2_offset[:] = np.nan

        # replace nan values with values of temp_o2_offset
        gdap_offsets_n_temp['DOXY_ADJUSTED_offset_trimmed'][gdap_offsets_n_temp['main_float_wmo']==[n]] = temp_o2_offset

        # append each temp_o2_offset to the new DOXY_ADJUSTED_offset_trimmed vector
        # DOXY_ADJUSTED_offset_trimmed.append(temp_o2_offset.values)
        #print(len(DOXY_ADJUSTED_offset_trimmed))
        # break

    # concatenate all vectors within DOXY_ADJUSTED_offset_trimmed (each represents one WMO)
    # result_vector = np.concatenate(DOXY_ADJUSTED_offset_trimmed)
    # # convert to Xarray DataArray
    # result_da = xr.DataArray(result_vector, dims='N_CROSSOVERS', coords={'N_CROSSOVERS': gdap_offsets_n_temp['N_CROSSOVERS']})
    # # add to glodap_offsets
    # gdap_offsets_n_temp['DOXY_ADJUSTED_offset_trimmed']=result_da
    # print(glodap_offsets)

    # calculate mean DOXY_ADJUSTED_offsets for different day ranges
    if len(float_age_bins)>0: # 0 = skip and do not apply
        print('in age section')
        # create new variables that will be the mean offsets for different time ranges
        empty_data = np.empty(len(gdap_offsets_n_temp['N_CROSSOVERS']))
        empty_data[:] = np.nan
        new_data_array = xr.DataArray(empty_data, coords={'N_CROSSOVERS': gdap_offsets_n_temp['N_CROSSOVERS']}, dims=['N_CROSSOVERS'])

        gdap_offsets_n_temp['First_Float_Profile_Date'] = new_data_array.copy()  # for storing the first date of first float profile

        for fa in range(len(float_age_bins)-1):
            gdap_offsets_n_temp['DOXY_ADJUSTED_offset_' + 'age_' + str(float_age_bins[fa])] = new_data_array.copy()
            gdap_offsets_n_temp['DOXY_ADJUSTED_offset_' + 'age_' + str(float_age_bins[fa]) + '_count'] = new_data_array.copy()

        wmo_list = np.unique(gdap_offsets_n_temp.main_float_wmo)
        # loop through all floats
        for wmo_n in wmo_list:
            # load Sprof file, get date of first profile
            argo_n = xr.open_dataset(argo_path + str(wmo_n) + '_Sprof.nc')
            first_profile_date = argo_n.JULD[0].values

            gdap_offsets_n_temp['First_Float_Profile_Date'][gdap_offsets_n_temp.main_float_wmo==wmo_n] = first_profile_date

            # calculate offset time relative to first deployment
            time_since_first_ns = gdap_offsets_n_temp.main_float_juld[gdap_offsets_n_temp.main_float_wmo==wmo_n].values-first_profile_date
            time_since_first_days = time_since_first_ns / np.timedelta64(1, 'D')

            # save out o2 offsets for this float
            temp_o2_offset = gdap_offsets_n_temp.DOXY_ADJUSTED_offset[gdap_offsets_n_temp.main_float_wmo==wmo_n]


            # now loop through float_age_bins to calculate mean offsets for each age range 
            for fa in range(len(float_age_bins)-1):
                # find ages within range
                age_index = np.logical_and(time_since_first_days>=float_age_bins[fa], time_since_first_days<float_age_bins[fa+1])

                # save the mean o2 offset where age_index is true to the correct age variable for that float
                if not np.all(np.isnan(temp_o2_offset.where(age_index).values)):
                    gdap_offsets_n_temp['DOXY_ADJUSTED_offset_' + 'age_' + str(float_age_bins[fa])]\
                        [gdap_offsets_n_temp.main_float_wmo==wmo_n] = np.nanmean(temp_o2_offset.where(age_index).values)
                gdap_offsets_n_temp['DOXY_ADJUSTED_offset_' + 'age_' + str(float_age_bins[fa]) + '_count']\
                                    [gdap_offsets_n_temp.main_float_wmo==wmo_n] = (temp_o2_offset.where(age_index).count())


    # save a copy of gdap_offsets_n_temp for later plotting / analysis:
    gdap_offsets_n_temp.to_netcdf(output_dir+out_filename+ '_all_offsets_depth_grouped_' + \
                                   'year_filt_' + str(year_filt) +'_' + str(year_plus_minus) + '_level_' + str(pressure_level_min) + '.nc')
    
    # then group again by WMO, now with the new variable:
    offsets_g = gdap_offsets_n_temp.groupby(gdap_offsets_n_temp.main_float_wmo)

    glodap_offsets_mean = offsets_g.mean(skipna='True')

    #initialise metadata DataArrays
    varnames = ['o2_calib_comment','o2_calib_equation','o2_calib_coeff','o2_calib_group','o2_calib_air_group',
            'o2_calib_drift','project_name','plat_type','data_centre']
    coord= range(glodap_offsets_mean.main_float_wmo.shape[0])
    for v in varnames:
        glodap_offsets_mean[v] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

    #also create groups by sensor 
    glodap_offsets_mean['pH_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
    glodap_offsets_mean['pH_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
    glodap_offsets_mean['nitrate_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
    glodap_offsets_mean['nitrate_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
    glodap_offsets_mean['DOXY_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)
    glodap_offsets_mean['DOXY_sensor'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

    #group for floats that have some under ice profiles with some air calibrated some not
    glodap_offsets_mean['ice_group'] = xr.DataArray(dims=['main_float_wmo'],coords={'main_float_wmo':coord}).astype(str)

    num_crossovers = glodap_offsets_mean.main_float_wmo.size

    # loop through each wmo in offsets
    for n,w in enumerate(glodap_offsets_mean.main_float_wmo):
        
        #find full float file matching offset
        wmo_n = w.values
        fn = argo_path + str(wmo_n) + '_Sprof.nc'
        float_og = xr.open_dataset(fn)
        
        #also load meta file  for same float
        fn1 = argo_path + str(wmo_n) + '_meta.nc'
        meta_og = xr.open_dataset(fn1)

        #option for multiple calibrations at different times- default to use most recent only for now?
        #O2 calibration
        #if len(float_og.N_CALIB) >1 :
            #print(str(g.values)+' has more than 1 calibration')
        
        no_cal = []
        if 'DOXY_ADJUSTED' in float_og.keys():
            is_doxy = 'DOXY'
            
            #get O2 sensor type
            #check sensor index
            sensors = meta_og.SENSOR.values
            sind = [idx for idx, s in enumerate(sensors) if 'DOXY' in s.decode("utf-8")] 
            #grab sensor model for index
            if len(sind):
                o2_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
            else:
                o2_sensor = 'not listed'

            cal_str = float_og.STATION_PARAMETERS.values.astype(str)[0]
            o2_ind = [idx for idx, s in enumerate(cal_str) if 'DOXY' in s]
            if len(o2_ind):

                
                #get calibration info
                o2_cal_full = float_og.SCIENTIFIC_CALIB_COMMENT.values[:,-1,o2_ind[0]]
                o2_cal_date = float_og.SCIENTIFIC_CALIB_DATE.values[:,-1,o2_ind[0]]
                o2_cal_unique = np.unique(o2_cal_full)
                
                #check if any unique cal comments are bad/no cal
                isbad = []
                for i in range(len(o2_cal_unique)):
                    o2_cal_i = o2_cal_unique[i].decode("utf-8")
                    if any(substring in o2_cal_i for substring in bad_cal_list) or \
                    any(substring in o2_cal_i for substring in no_cal_list) or \
                    any(substring in o2_cal_i for substring in blank_cal):
                        isbad.append(i)
                #remove bad indices
                o2_cal_unique = np.delete(o2_cal_unique,isbad)

                #check if more than one unique calibration comment remaining after excluding bad data
                if len(o2_cal_unique)>1:
                    print('multiple cal comments')  
                    print(wmo_n)
                    #how to integrate multiple comments? Do any contain different groupings?
                    o2_cal = float_og.SCIENTIFIC_CALIB_COMMENT.values[0,-1,o2_ind[0]].decode("utf-8")
                    if np.any(float_og.POSITION_QC.values==8):
                        under_ice = 'some air cal, some under ice'
                        print('under ice')
                    else:
                        under_ice = 'no ice'
                        
                else:
                    o2_cal = float_og.SCIENTIFIC_CALIB_COMMENT.values[0,-1,o2_ind[0]].decode("utf-8")
                o2_eq = float_og.SCIENTIFIC_CALIB_EQUATION.values[0,-1,o2_ind[0]].decode("utf-8")
                o2_co = float_og.SCIENTIFIC_CALIB_COEFFICIENT.values[0,-1,o2_ind[0]].decode("utf-8")
                under_ice = 'no ice'
            else:
                o2_cal = 'empty O2 calib comment'
                o2_eq = 'none'
                o2_co = 'none'

        else:
            is_doxy = 'no DOXY'
            o2_cal = 'no O2'
            o2_eq = 'none'
            o2_co = 'none'
            o2_sensor = 'none'
        

        ######## group calibration comments into categories
        
        #o2 calib comment groups
        if any(substring in o2_cal for substring in bad_cal_list):
            o2_group = 'bad'
        elif any(substring in o2_cal for substring in no_cal_list):
            o2_group = 'no cal'
        elif any(substring in o2_cal for substring in noair_cal_surf_list): # should catch the "no air cal" ones before they get put in air-cal list...
            o2_group = 'no air surface cal'
        elif any(substring in o2_cal for substring in air_cal_list):
            o2_group = 'air cal'
        elif any(substring in o2_cal for substring in noair_cal_subsurf_list):
            o2_group = 'no air subsurface cal'
        elif any(substring in o2_cal for substring in noair_cal_funcofdoxy_list):
            o2_group = 'no air function of DOXY cal'
        elif any(substring in o2_cal for substring in noair_cal_unspec_list):
            o2_group = 'no air unspecified'
        else:
            o2_group = 'ungrouped'
        
        
        #group O2 air cal and no air cal meta groups
        if any(substring in o2_cal for substring in air_cal_list):
            o2_air_group = 'air cal'
        elif any(substring in o2_cal for substring in noair_cal_combined_list):
            o2_air_group = 'no air cal'
        else:
            o2_air_group = 'no cal/bad'
            
        #separate groups for drift/no drift
        if any(substring in o2_cal for substring in noair_cal_withdrift_list):
            o2_group_drift = 'drift corrected'
        else:
            o2_group_drift = 'no drift corrected'
        
        if 'PH_IN_SITU_TOTAL_ADJUSTED' in float_og.keys():
            is_pH = 'pH'
            #get pH sensor type
            #check sensor index
            sensors = meta_og.SENSOR.values
            sind = [idx for idx, s in enumerate(sensors) if 'PH' in s.decode("utf-8")] 
            #grab sensor model for index
            pH_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
        else:
            is_pH = 'no pH'
            pH_sensor = 'none'
        
        if 'NITRATE_ADJUSTED' in float_og.keys():
            is_nitrate = 'nitrate'
            #get nitrate sensor type
            #check sensor index
            sensors = meta_og.SENSOR.values
            sind = [idx for idx, s in enumerate(sensors) if 'NITRATE' in s.decode("utf-8")] 
            #grab sensor model for index
            nitrate_sensor = meta_og.SENSOR_MODEL.values[sind[0]].decode("utf-8")
        else:
            is_nitrate = 'no nitrate'
            nitrate_sensor = 'none'
        

        ###### get sensor information from meta files

        
        
        ##################################################
        #add metadata to glodap_offset Dataset
        glodap_offsets_mean.o2_calib_comment[n] = o2_cal
        glodap_offsets_mean.o2_calib_equation[n] = o2_eq
        glodap_offsets_mean.o2_calib_coeff[n] = o2_co
        glodap_offsets_mean.o2_calib_group[n] = o2_group
        glodap_offsets_mean.o2_calib_air_group[n] = o2_air_group
        glodap_offsets_mean.o2_calib_drift[n] = o2_group_drift
        glodap_offsets_mean.project_name[n] = float_og.PROJECT_NAME.values.astype(str)[0]
        glodap_offsets_mean.plat_type[n] = float_og.PLATFORM_TYPE.values.astype(str)[0]
        glodap_offsets_mean.data_centre[n] = float_og.DATA_CENTRE.values.astype(str)[0]
        glodap_offsets_mean.DOXY_group[n] = is_doxy
        glodap_offsets_mean.DOXY_sensor[n] = o2_sensor
        glodap_offsets_mean.pH_group[n] = is_pH
        glodap_offsets_mean.pH_sensor[n] = pH_sensor
        glodap_offsets_mean.nitrate_group[n] = is_nitrate
        glodap_offsets_mean.nitrate_sensor[n] = nitrate_sensor
        glodap_offsets_mean.ice_group[n] = under_ice
        
        print(str(n) + ' out of ' + str(num_crossovers))
        # break
    glodap_offsets_mean.to_netcdf(output_dir+out_filename+ '_floatmean_withcalibration_depth_grouped_' + \
                                   'year_filt_' + str(year_filt) +'_' + str(year_plus_minus) + '_level_' + str(pressure_level_min) + '.nc')
    # trimmed_means[f'level_{pressure_level_min}'] = glodap_offsets_mean