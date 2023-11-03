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

# ## Plotting script for analysing offsets and applying bias adjustments

import numpy as np
import pandas as pd
import xarray as xr
import glob, os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from datetime import datetime, date, time, timedelta
import time
import gsw
import scipy.stats as stats
import matplotlib.dates as mdates
# import statsmodels.api as sm
from scipy import interpolate


# +
glodap_offsets_filename = 'glodap_offsets_100km_1450_to_2000_100m_0.005dens_0.005spice_4.nc'
# glodap_offsets_filename = 'glodap_offsets_100km_2_to_50_50m_0.1dens_0.1spice_5.nc'
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
        
# Set the paths
output_dir = 'output/'
data_dir = 'data/'

# Check for a glodap_offsets_plots directory, create if it does not exist
offset_dir = output_dir + 'glodap_offset_plots/'
if not os.path.isdir(offset_dir):
    os.mkdir(offset_dir)
individual_plot_dir = offset_dir + 'individual_floats' + glodap_offsets_filename[14:-3] + '/'
if not os.path.isdir(individual_plot_dir):
    os.mkdir(individual_plot_dir)


# +
def test_stat(y, iteration, verboseTF):
    std_dev = np.std(y)
    avg_y = np.mean(y)
    abs_val_minus_avg = abs(y - avg_y)
    max_of_deviations = max(abs_val_minus_avg)
    max_ind = np.argmax(abs_val_minus_avg)
    cal = max_of_deviations/ std_dev
    if verboseTF:
        print('Test {}'.format(iteration))
        print("Test Statistics Value(R{}) : {}".format(iteration,cal))
    return cal, max_ind

def calculate_critical_value(size, alpha, iteration, verboseTF):
    t_dist = stats.t.ppf(1 - alpha / (2 * size), size - 2)
    numerator = (size - 1) * np.sqrt(np.square(t_dist))
    denominator = np.sqrt(size) * np.sqrt(size - 2 + np.square(t_dist))
    critical_value = numerator / denominator
    if verboseTF:
        print("Critical Value(λ{}): {}".format(iteration, critical_value))
    return critical_value

def check_values(R, C, inp, max_index, iteration):
    if R > C:
        print('{} is an outlier. R{} > λ{}: {:.4f} > {:.4f} \n'.format(inp[max_index],iteration, iteration, R, C))
    else:
        print('{} is not an outlier. R{}> λ{}: {:.4f} > {:.4f} \n'.format(inp[max_index],iteration, iteration, R, C))




def ESD_Test(input_series, alpha, max_outliers, verboseTF, NoOutputTF):
    stats = []
    critical_vals = []
    tested_values = []
    max_i = 0
    for iterations in range(1, max_outliers + 1):
        stat, max_index = test_stat(input_series, iterations, verboseTF)
        critical = calculate_critical_value(len(input_series), alpha, iterations, verboseTF)
        if verboseTF:
            check_values(stat, critical, input_series, max_index, iterations)
        tested_values.append(input_series[max_index])
        input_series = np.delete(input_series, max_index)
        critical_vals.append(critical)
        stats.append(stat)
        if stat > critical:
            max_i = iterations
    if ~NoOutputTF:
        print('H0:  there are no outliers in the data')
        print('Ha:  there are up to 10 outliers in the data')
        print('')
        print('Significance level:  α = {}'.format(alpha))
        print('Critical region:  Reject H0 if Ri > critical value')
        print('Ri: Test statistic')
        print('λi: Critical Value')
        print(' ')
    df = pd.DataFrame({'i' :range(1, max_outliers + 1), 'Ri': stats, 'λi': critical_vals , 'Vals': tested_values})
    
    def highlight_max(x):
        if x.i == max_i:
            return ['background-color: yellow']*4
        else:
            return ['background-color: white']*4
    df.index = df.index + 1
    print('Number of outliers {}'.format(max_i))
    print(max_i)
    return  df.style.apply(highlight_max, axis = 1), max_i, tested_values


# -

def regress_confidence_sokal_rohlf(X, Y, alpha):

    q = 1-alpha
    # calculate averages
    X_bar = np.mean(X)
    Y_bar = np.mean(Y)

    #difference of datapoints from mean
    x_small = X - X_bar
    y_small = Y - Y_bar

    # sum of the square differences in x
    sum_x2 = np.nansum(np.square(x_small))

    # multiple differences from mean to get xy and sum
    xy = x_small*y_small
    sum_xy = np.nansum(xy)

    # regression coefficient (slope of the regression)
    b_yx = sum_xy/sum_x2

    # y intercept 
    a = Y_bar - b_yx*X_bar

    # predicted Y based on regression
    Y_hat = b_yx*X + a

    # unexplained sum of squares
    d2_yx = np.sum(np.square(Y-Y_hat))

    # explained variance:
    s2_y_hat = np.sum(np.square(Y_hat - Y_bar))/(len(Y)-1)

    # unexplained variance:
    dof = (len(Y) - 2)
    s2_yx = d2_yx/dof;

    # total variance:
    s2_Y = np.sum(np.square(Y - Y_bar))/(len(Y)-1)

    # coefficient of determination (r2)
    r2 = s2_y_hat/s2_Y

    # standard error of the regression coefficient
    sb = np.sqrt(s2_yx/sum_x2);

    # regression coefficient (slope) different than zero?
    t_s = (b_yx - 0)/sb

    # t stat
    t_stat = stats.t.ppf(1-q/2, dof)

    CI_alpha_slope = [b_yx - t_stat*sb, b_yx + t_stat*sb]

    # standard error of Y_hat for a given value of X (every value)
    sy_hat = (s2_yx*(1/len(Y) + ((X-X_bar)**2)/sum_x2))**0.5

    # upper confidence interval
    y_err = t_stat*sy_hat

    return b_yx, a, r2, CI_alpha_slope, t_stat*sb/2, y_err


# ### 1. Load offsets and look at offset summary plots for each float
# This step helps to 1) identify patterns in float offsets in particular floats so that we can go back and refine offset criteria and 2) make final decisions about whether a float should be adjusted or not.
#
# First load saved glodap offsets

# +
if 'glodap_offsets' in locals():
    glodap_offsets.close()

glodap_offsets = xr.load_dataset(output_dir+glodap_offsets_filename)
# -

# plot histograms of offsets for each main float with all crossovers

glodap_offsets

#group by main float wmo
offsets_g = glodap_offsets.groupby(glodap_offsets.main_float_wmo)

# Determine outliers and remove, creating a "DOXY_ADJUSTED_offset_trimmed" variable

# +
# adding option to filter by time of year as well - for use in surface data
time_filt = 0
filt_days = 10
DOXY_ADJUSTED_offset_trimmed = []
for n,g in offsets_g:
    
    # run a GESD test using "test_num" number of possible outliers
    test_num = int((len(g.DOXY_ADJUSTED_offset.dropna(dim="N_CROSSOVERS", how="any").values)*.1))
    ESD_test_out = ESD_Test(g.DOXY_ADJUSTED_offset.dropna(dim="N_CROSSOVERS", how="any").values, 0.05, test_num, False, True)


    # create temp_o2_offest to set all datapoints to nans that the GESD test says are outliers
    if time_filt==1:
        within_days = np.logical_or(np.abs(g.main_float_juld.dt.dayofyear - g.glodap_datetime.dt.dayofyear)<=filt_days, 
                           np.abs(g.main_float_juld.dt.dayofyear - g.glodap_datetime.dt.dayofyear)>(365-filt_days)) 
        temp_o2_offset = g.DOXY_ADJUSTED_offset.where(within_days)
    else:
        temp_o2_offset = g.DOXY_ADJUSTED_offset
    for a in range(0, ESD_test_out[1]):
        temp_o2_offset = temp_o2_offset.where(temp_o2_offset != ESD_test_out[2][a])
        
    # append each temp_o2_offset to the new DOXY_ADJUSTED_offset_trimmed vector
    DOXY_ADJUSTED_offset_trimmed.append(temp_o2_offset.values)
    #print(len(DOXY_ADJUSTED_offset_trimmed))

# concatenate all vectors within DOXY_ADJUSTED_offset_trimmed (each represents one WMO)
result_vector = np.concatenate(DOXY_ADJUSTED_offset_trimmed)
# convert to Xarray DataArray
result_da = xr.DataArray(result_vector, dims='N_CROSSOVERS', coords={'N_CROSSOVERS': glodap_offsets['N_CROSSOVERS']})
# add to glodap_offsets
glodap_offsets['DOXY_ADJUSTED_offset_trimmed']=result_da
print(glodap_offsets)

# then group again by WMO, now with the new variable:
offsets_g = glodap_offsets.groupby(glodap_offsets.main_float_wmo)
# -



# +
# calculate whether change in ocean oxygen content 
# as determined by glodap could be responsible for float-glodap difference

# create a list to put in True / False if the glodap offsets intersect zero at float mid date
glodap_drift_possible_list = []

for n,g in offsets_g:
    g_ox = g.where(~np.isnan(g.DOXY_ADJUSTED_offset_trimmed), drop=True)
    
    g_ox_sorted = g_ox.sortby("glodap_datetime")
    
    X_series = pd.Series(mdates.date2num(g_ox_sorted.glodap_datetime))
    if len(X_series)<5: # too short to realistically do a regression
        uncert_min = np.nan
        uncert_max = np.nan
    else:
        Y_series = pd.Series(g_ox_sorted.DOXY_ADJUSTED_offset_trimmed.values)
        CI_alpha = 0.95
        b_yx, a, r2, CI_alpha_slope, ttt, y_err = regress_confidence_sokal_rohlf(X_series, Y_series, CI_alpha)

    #     y_est = a+ X_series*b_yx
    #     plt.plot(X_series, y_est, color='blue', label="y = "+ str(b_yx.round(decimals=4)) +"x + " + str(a.round(decimals=2)), linestyle = '-')
    #     plt.fill_between(X_series, y_est-y_err, y_est+y_err, color='blue', alpha=0.5, label='Confidence Interval')    

        # extrapolate regression and CI if needed to intersect float_mid_date
        # should apply the slope of the 
        float_mid_date = mdates.date2num(g.main_float_juld.mean()) # mean float date in number

        if np.max(X_series)<float_mid_date:
            X_extend = X_series.append(pd.Series(float_mid_date))

            y_extend = a+ X_extend*b_yx

            X_series_last_third = X_series.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
            Y_err_last_third = y_err.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
            b_err, a_err, _, _, _, _ = regress_confidence_sokal_rohlf(X_series_last_third, Y_err_last_third, CI_alpha)
            extrap_error = X_extend.iloc[-1]*b_err+a_err
            y_err_extend = y_err.append(pd.Series(extrap_error))

    #         y_err_extend = y_err.append(pd.Series(y_err.iloc[-1]))

#             plt.plot(X_extend, y_extend, color='blue', linestyle = '--')
#             plt.fill_between(X_extend, y_extend-y_err_extend, y_extend+y_err_extend, color='blue', alpha=0.25)    
            uncert_min = y_extend.iloc[-1] - y_err_extend.iloc[-1]
    #         print(uncert_min)
            uncert_max = y_extend.iloc[-1] + y_err_extend.iloc[-1]
    #         print(uncert_max)
        elif np.min(X_series)>float_mid_date:
            X_extend = pd.Series(float_mid_date).append(X_series)
            y_extend = a+ X_extend*b_yx
            X_series_first_half = X_series.iloc[0:np.int64(np.round(len(X_series)*1/2))]
            Y_err_first_half = y_err.iloc[1:np.int64(np.round(len(X_series)*1/2))]
            b_err, a_err, _, _, _, _ = regress_confidence_sokal_rohlf(X_series_last_third, Y_err_last_third, CI_alpha)
            extrap_error = X_extend.iloc[0]*b_err+a_err
            y_err_extend = pd.Series(extrap_error).append(y_err)

#             plt.plot(X_extend, y_extend, color='blue', linestyle = '--')
#             plt.fill_between(X_extend, y_extend-y_err_extend, y_extend+y_err_extend, color='blue', alpha=0.25)    
            uncert_min = y_extend.iloc[0] - y_err_extend.iloc[0]
            uncert_max = y_extend.iloc[0] + y_err_extend.iloc[0]
        else:
            y_float = float_mid_date*b_yx+a
            f = interpolate.interp1d(X_series,y_err)
            y_float_err = f(float_mid_date)

            uncert_min = y_float - y_float_err
            uncert_max = y_float + y_float_err
    #         print(uncert_min)
    #         print(uncert_max)

        
    if np.isnan(uncert_min):
        glodap_drift_possible = np.nan
    elif np.logical_and(uncert_min<0, uncert_max>0):
        glodap_drift_possible = True
    else:
        glodap_drift_possible = False

    glodap_drift_possible_list.append(glodap_drift_possible)
    

# +
# Depth dependency of offsets:
# 1. Bin by depth - 50m bins? 100m bins?
# 2. Save a dataframe with wmo, depth, o2 conc, o2 offset, o2 offset minus mean o2 offset
# That should give you what you need to plot all float offsets vs. depth and concentration
press_bin_width = 10
press_bins = np.linspace(np.int32(glodap_offsets['p_compare_min'].item()), 
                   np.int32(glodap_offsets['p_compare_max'].item()), press_bin_width)

press_response_all_df = pd.DataFrame(columns=['wmo','avg_depth','o2_conc', 'o2_offset', 'o2_offset_minus_mean'])
num_count = 0
for n,g in offsets_g:
    num_count = num_count+1    
    print(str(n)+ ' ' + str(num_count))
    g_ox = g.where(~np.isnan(g.DOXY_ADJUSTED_offset_trimmed), drop=True)

    # loop through press_bins, finding pressures that fit and averaging everything
    press_response_data = []  # List to store data

    for idx in range(0, len(press_bins)-1):
        depth_index = np.logical_and(g_ox['PRES_ADJUSTED_float']>=press_bins[idx], 
                                     g_ox['PRES_ADJUSTED_float']<press_bins[idx+1])
        if sum(depth_index)==0:
            continue

        # make a new dataframe for this pressure bin

        press_response_data.append({
                    'wmo': np.int32(g_ox['main_float_wmo'][0].item()),
                    'avg_depth': np.mean(g_ox['PRES_ADJUSTED_float'][depth_index]).values.tolist(),
                    'o2_conc': np.mean(g_ox['DOXY_ADJUSTED_float'][depth_index]).values.tolist(),
                    'o2_offset': np.mean(g_ox['DOXY_ADJUSTED_offset_trimmed'][depth_index]).values.tolist(),
                    'o2_offset_minus_mean': (np.mean(g_ox['DOXY_ADJUSTED_offset_trimmed'][depth_index]).values - \
                                             np.mean(g_ox['DOXY_ADJUSTED_offset_trimmed']).values).tolist()

            })
    # Convert the list of dictionaries into a DataFrame
    press_response_df = pd.DataFrame(press_response_data)
    # concatenate dataframes
    press_response_all_df = pd.concat([press_response_all_df, press_response_df], ignore_index=True)


# +
plt.plot(press_response_all_df['o2_offset_minus_mean'], press_response_all_df['avg_depth'], 'x')
plt.xlim([-20, 20])
plt.grid()
# plt.yscale(reversed)

print(press_response_all_df['o2_offset_minus_mean'].mean())
for idx in range(0, len(press_bins)-1):
    depth_index = np.logical_and(press_response_all_df['avg_depth']>=press_bins[idx], 
                                     press_response_all_df['avg_depth']<press_bins[idx+1])
    plt.plot(np.mean(press_response_all_df['o2_offset_minus_mean'][depth_index]),
                 (press_bins[idx+1]-press_bins[idx])/2+press_bins[idx], 'ro')
    plt.gca().invert_yaxis()

    plt.xlabel('Pressure binned offsets minus mean for each float')
    plt.ylabel('Pressure')
    plt.savefig(offset_dir+ 'Offsets_v_pressure.png')

# -

press_response_all_df['avg_depth'][press_response_all_df['wmo']==1901212]

# +
sum_true = 0
sum_false = 0
for idx, value in enumerate(glodap_drift_possible_list):
    if value is True:
        sum_true = sum_true+1
    elif value is False:
        sum_false = sum_false+1
        
print('true ' + str(sum_true))
print('false ' + str(sum_false))
# -

glodap_drift_possible_list[:]==0
value

offsets_g

# +
# show_plot = False
#loop through each float
fig = plt.figure(figsize=(16,16))

# create a list to put in True / False if the glodap offsets intersect zero at float mid date
glodap_drift_possible_list = []

for n,g in offsets_g:
    print(n)
    ncross = len(g.DOXY_ADJUSTED_offset)

    #if doxy_adjusted_offset is nan, skip
    if np.all(np.isnan(g.DOXY_ADJUSTED_offset.values)):
        continue

    #add crossover location map
    axn = plt.subplot(4,4,1)
    axn.plot(g.main_float_longitude,g.main_float_latitude,'g+',markersize=10,label='Current float')
    #glodap
    glodap_lon = xr.where(g.glodap_longitude>180,g.glodap_longitude-360.,g.glodap_longitude)
    axn.plot(glodap_lon,g.glodap_latitude,'ro',label = 'Glodap',markersize=10)
    axn.legend()
    axn.set_title('WMO: %d; N: %d; Dist filt %d; Depth min %d max %d\nDepth filt %d; Dens %.4f; Spice %.4f' % 
                  (g.main_float_wmo.values[0],ncross, glodap_offsets["dist"],
                   glodap_offsets["p_compare_min"], 
                   glodap_offsets["p_compare_max"],glodap_offsets["delta_press"],
                   glodap_offsets["delta_dens"],glodap_offsets["delta_spice"]))

    #time
    axn = plt.subplot(4,4,(2,4))
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset,'rx',label='float')
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset_trimmed,'ro',label='float_outliers_removed', markerfacecolor="None")
    axn.plot(g.glodap_datetime,g.DOXY_ADJUSTED_offset,'bx',label='GDAP')
    axn.plot(g.glodap_datetime,g.DOXY_ADJUSTED_offset_trimmed,'bo', markerfacecolor="None")

    # fit regressions
    g_ox = g.where(~np.isnan(g.DOXY_ADJUSTED_offset_trimmed), drop=True)
    m1, b1 = np.polyfit(mdates.date2num(g_ox.main_float_juld),g_ox.DOXY_ADJUSTED_offset_trimmed, 1)
#     t = np.arange(datetime(1980,1,1), datetime(2022,1,1), timedelta(days=1)).astype(datetime)
#     tnum = mdates.date2num(t)
    axn.plot(mdates.date2num(g_ox.main_float_juld), m1*mdates.date2num(g_ox.main_float_juld)+b1, 
             "-", c="r", label="y = "+ str(m1.round(decimals=4)) +"x + " + str(b1.round(decimals=2)))
    axn.plot(g_ox.main_float_juld, m1*mdates.date2num(g_ox.main_float_juld)+b1, c="r")
    
    g_ox_sorted = g_ox.sortby("glodap_datetime")
    
    X_series = pd.Series(mdates.date2num(g_ox_sorted.glodap_datetime))
    Y_series = pd.Series(g_ox_sorted.DOXY_ADJUSTED_offset_trimmed.values)
    CI_alpha = 0.95
    b_yx, a, r2, CI_alpha_slope, ttt, y_err = regress_confidence_sokal_rohlf(X_series, Y_series, CI_alpha)

#     t = np.arange(datetime(1980,1,1), datetime(2022,1,1), timedelta(days=1)).astype(datetime)
#     tnum = mdates.date2num(t)
    
    y_est = a+ X_series*b_yx
    plt.plot(X_series, y_est, color='blue', label="y = "+ str(b_yx.round(decimals=4)) +"x + " + str(a.round(decimals=2)), linestyle = '-')
    plt.fill_between(X_series, y_est-y_err, y_est+y_err, color='blue', alpha=0.5, label='Confidence Interval')    
    
    # extrapolate regression and CI if needed to intersect float_mid_date
    # should apply the slope of the 
    float_mid_date = mdates.date2num(g.main_float_juld.mean()) # mean float date in number
    if len(X_series)<5: # too short to realistically do a regression
        uncert_min = -100
        uncert_max = 100
    elif np.max(X_series)<float_mid_date:
        X_extend = X_series.append(pd.Series(float_mid_date))

        y_extend = a+ X_extend*b_yx
        
        X_series_last_third = X_series.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
        Y_err_last_third = y_err.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
        b_err, a_err, _, _, _, _ = regress_confidence_sokal_rohlf(X_series_last_third, Y_err_last_third, CI_alpha)
        extrap_error = X_extend.iloc[-1]*b_err+a_err
        y_err_extend = y_err.append(pd.Series(extrap_error))
        
#         y_err_extend = y_err.append(pd.Series(y_err.iloc[-1]))
         
        plt.plot(X_extend, y_extend, color='blue', linestyle = '--')
        plt.fill_between(X_extend, y_extend-y_err_extend, y_extend+y_err_extend, color='blue', alpha=0.25)    
        uncert_min = y_extend.iloc[-1] - y_err_extend.iloc[-1]
#         print(uncert_min)
        uncert_max = y_extend.iloc[-1] + y_err_extend.iloc[-1]
#         print(uncert_max)
    elif np.min(X_series)>float_mid_date:
        X_extend = pd.Series(float_mid_date).append(X_series)
        y_extend = a+ X_extend*b_yx
        X_series_first_half = X_series.iloc[0:np.int64(np.round(len(X_series)*1/2))]
        Y_err_first_half = y_err.iloc[1:np.int64(np.round(len(X_series)*1/2))]
        b_err, a_err, _, _, _, _ = regress_confidence_sokal_rohlf(X_series_last_third, Y_err_last_third, CI_alpha)
        extrap_error = X_extend.iloc[0]*b_err+a_err
        y_err_extend = pd.Series(extrap_error).append(y_err)
        
        plt.plot(X_extend, y_extend, color='blue', linestyle = '--')
        plt.fill_between(X_extend, y_extend-y_err_extend, y_extend+y_err_extend, color='blue', alpha=0.25)    
        uncert_min = y_extend.iloc[0] - y_err_extend.iloc[0]
        uncert_max = y_extend.iloc[0] + y_err_extend.iloc[0]
    else:
        y_float = float_mid_date*b_yx+a
        f = interpolate.interp1d(X_series,y_err)
        y_float_err = f(float_mid_date)

        uncert_min = y_float - y_float_err
        uncert_max = y_float + y_float_err
#         print(uncert_min)
#         print(uncert_max)

        
    if np.logical_and(uncert_min<0, uncert_max>0):
        glodap_drift_possible = True
    else:
        glodap_drift_possible = False

    glodap_drift_possible_list.append(glodap_drift_possible)
#     print(glodap_drift_possible)
#     m2, b2 = np.polyfit(mdates.date2num(g_ox.glodap_datetime),g_ox.DOXY_ADJUSTED_offset_trimmed, 1)
    
#     axn.plot(t, m2*tnum+b2, "--", c="b")
#     axn.plot(X_series, m2*mdates.date2num(g_ox.glodap_datetime)+b2, c="b")

    axn.axhline(y=0, color='k', linestyle='--')
    axn.legend()
    joined_dates = mdates.date2num(sorted(np.append(g.glodap_datetime.values, g.main_float_juld.values)))
    axn.set_xlim(mdates.num2date(joined_dates[0] - 365), mdates.num2date(joined_dates[-1] + 365))
    try:
        axn.set_ylim(min(g.DOXY_ADJUSTED_offset.values)  - 5,  max(g.DOXY_ADJUSTED_offset.values) + 5)
    except ValueError:
        pass

    axn.set_title('O2 Offsets vs date, Real ocean O2 change? ' + str(glodap_drift_possible))
#     axn.text(0.01, 0.01, "y = "+ str(m1.round(decimals=4)) +"x + " + str(b1.round(decimals=2)) + "\n"
#              "y = "+ str(m2.round(decimals=4)) +"x + " + str(b2.round(decimals=2)),
#              verticalalignment='bottom', horizontalalignment='left',
#              transform=axn.transAxes)

    # plot offset histograms
    # Temperature
    axn = plt.subplot(4,5,6)
    g_plot = g.TEMP_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('TEMP %.2f +/- %.2f' % (g_mean,g_std))
    plt.grid()

    #Salinity
    axn = plt.subplot(4,5,7)
    g_plot = g.PSAL_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('PSAL %.2f +/- %.2f' % (g_mean,g_std))
    plt.grid()

    # Oxygen histogram
    axn = plt.subplot(4,5,8)
    g_plot = g.DOXY_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values)
    t_stat, p_value = stats.ttest_1samp(a=g_plot, popmean=0, nan_policy='omit') ############
    
    g_std = np.nanstd(g_plot.values).round(decimals=2)

    g_plot_trim = g.DOXY_ADJUSTED_offset_trimmed
    g_mean_trim = np.nanmean(g_plot_trim.values)
    g_std_trim = np.nanstd(g_plot_trim.values).round(decimals=2)

    t_stat_trim, p_value_trim = stats.ttest_1samp(a=g_plot_trim, popmean=0, nan_policy='omit') ############

    CI_99 = stats.norm.interval(alpha=0.99, 
                 loc=np.nanmean(g_plot_trim[~np.isnan(g_plot_trim)].values), 
                 scale = stats.sem(g_plot_trim[~np.isnan(g_plot_trim)].values))
    
    #g_ttest = stats.ttest_1samp(a=g_plot.values[~np.isnan(g_plot)], popmean=g_mean) ############
    if not np.all(np.isnan(g_plot)):
        all_o2_offsets_no_nan = g_plot
        all_o2_offsets_no_nan = all_o2_offsets_no_nan[~np.isnan(all_o2_offsets_no_nan)]
        hist1, bins = np.histogram(all_o2_offsets_no_nan)
        axn.hist(g_plot, bins=bins, color='b',alpha=1, label='all' + ' p: ' + str(p_value.round(5)))
        axn.hist(g_plot_trim, bins=bins,color='r',alpha=0.5, label='trimmed' + ' p: ' + str(p_value_trim.round(5)))

        axn.legend()

        
    axn.set_title('O2 n: %.1f, %.1f +/- %.1f' % (len(g_plot[~np.isnan(g_plot)]),g_mean,g_std))
   # axn.text(0.99, 0.99, "p_value: " + str(p_value.round(5)) + "\n"
   #      + "t-statistic: " + str(t_stat.round(3)),
   #      verticalalignment='top', horizontalalignment='right',
   #      transform=axn.transAxes)
    plt.grid()

    # oxygen histogram - trimmed data only:
    axn = plt.subplot(4,5,9)
    g_plot = g.DOXY_ADJUSTED_offset_trimmed
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot,color='r',alpha=0.5)
        plt.axvline(CI_99[1], color='k', linestyle='--', label='99% CI')
        plt.axvline(CI_99[0], color='k', linestyle='--')
        axn.legend()
    axn.set_title('O2 trim n: %.1f, %.1f +/- %.1f' % (len(g_plot_trim[~np.isnan(g_plot_trim)]), g_mean,g_std))
    plt.grid()
    
    # Nitrate
    axn = plt.subplot(4,5,10)
    g_plot = g.NITRATE_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values).round(decimals=2)
    g_std = np.nanstd(g_plot.values).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('NO3 %.2f +/- %.2f' % (g_mean,g_std))
    plt.grid()

    # pH
    axn = plt.subplot(4,4,9)
    g_plot = g.pH_25C_TOTAL_ADJUSTED_offset
    g_mean = np.nanmean(g_plot.values*1000).round(decimals=2)
    g_std = np.nanstd(g_plot.values*1000).round(decimals=2)
    if not np.all(np.isnan(g_plot)):
        axn.hist(g_plot*1000,color='b',alpha=0.5)
    axn.set_title('PH 25C %.2f +/- %.2f' % (g_mean,g_std))
    plt.grid()


    #O2 vs pressure
    axn = plt.subplot(4,4,10)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_float,'bx')
    axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PRES_ADJUSTED_float,'bo', markerfacecolor="none")

    axn.set_title('Float pres vs DOXY offset')
    plt.gca().invert_yaxis()
    axn.set_ylabel('pres')
    plt.grid()

    axn = plt.subplot(4,4,11)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_glodap,'bx')
    axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PRES_ADJUSTED_glodap,'bo', markerfacecolor="none")

    axn.set_title('GDAP pres vs DOXY offset')
    plt.gca().invert_yaxis()
    plt.grid()

    #vs density
    axn = plt.subplot(4,4,12)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PDENS_float,'bx')
    axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PDENS_float,'bo', markerfacecolor="none")

    axn.set_title('Float Pdens vs DOXY offset')
    axn.set_ylabel('dens')
    plt.grid()

    #axn = plt.subplot(4,4,12)
    #axn.plot(g.DOXY_ADJUSTED_offset,g.PDENS_glodap,'bx')
    #axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PDENS_glodap,'bo', markerfacecolor="none")

    #axn.set_title('GDAP Pdens vs DOXY offset')
    #plt.grid()

    #O2 offset vs o2 concentration
    axn = plt.subplot(4,4,13)
#     axn.plot(g.DOXY_ADJUSTED_offset,g.DOXY_ADJUSTED_float,'bx')
    # trying lots of things:
    float_time_days = mdates.date2num(g.main_float_juld)
    time_since_deploy = float_time_days - float_time_days[0]

#     cruise_nums = np.unique(g.glodap_cruise)
# #     cruise_num_scaled = np.zeros(len(g.glodap_cruise))
#     for idx, cr in enumerate(cruise_nums):
#         cruise_num_scaled[g.glodap_cruise==cruise_nums[idx]] = idx

    obs_inds = np.unique(g.glodap_obs_index)
    obs_inds_scaled = np.zeros(len(g.glodap_obs_index))
    for idx, cr in enumerate(obs_inds):
        obs_inds_scaled[g.glodap_obs_index==obs_inds[idx]] = idx
        
    cn = axn.scatter(g.DOXY_ADJUSTED_offset_trimmed, g.DOXY_ADJUSTED_float,30, obs_inds_scaled, 
                marker='o', cmap='tab20c')
    plt.colorbar(cn,label='?')
    
    axn.set_title('Float O2 vs DOXY offset')
    axn.set_ylabel('O2 concentration')
    plt.grid()

    axn = plt.subplot(4,4,14)
#     axn.plot(g.DOXY_ADJUSTED_offset,g.DOXY_ADJUSTED_float,'bx')
    # trying lots of things:
    float_time_days = mdates.date2num(g.main_float_juld)
    time_since_deploy = float_time_days - float_time_days[0]

    cn = axn.scatter(g.DOXY_ADJUSTED_offset_trimmed, g.DOXY_ADJUSTED_glodap,30, obs_inds_scaled, 
                marker='o', cmap='tab20c')
    plt.colorbar(cn,label='?')
    
    axn.set_title('Glodap O2 vs DOXY offset')
    axn.set_ylabel('O2 concentration')
    plt.grid()

    
    
    # Add density lines FLOAT
    smax = float(max(g.PSAL_ADJUSTED_float)) + 0.1
    smin = float(min(g.PSAL_ADJUSTED_float)) - 0.1
    tmax = float(max(g.TEMP_ADJUSTED_float)) + 0.1
    tmin = float(min(g.TEMP_ADJUSTED_float)) - 0.1
    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(tmin,tmax,num=int(round((tmax-tmin)*10 + 1, 0)))
    si = np.linspace(smin,smax,num=int(round((smax-smin)*10 + 1, 0)))
    dens = np.zeros((len(ti),len(si)))
    # Loop to fill in grid with densities
    for j in range(0,len(ti)):
        for i in range(0, len(si)):
            dens[j,i]=gsw.sigma0(si[i],ti[j])

    axn = plt.subplot(4,4,15)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    axn.clabel(CS, fontsize=12, inline=1)
    axn.scatter(g.PSAL_ADJUSTED_float,g.TEMP_ADJUSTED_float,c=g.DOXY_ADJUSTED_offset)
    axn.set_xlim(min(g.PSAL_ADJUSTED_float) - .1, max(g.PSAL_ADJUSTED_float) + .1)
    axn.set_ylim(min(g.TEMP_ADJUSTED_float) - 0.1, max(g.TEMP_ADJUSTED_float) + 0.1)
    axn.set_title('float T-S')

    # Add density lines GLODAP
    smax = float(max(g.PSAL_ADJUSTED_glodap)) + 0.1
    smin = float(min(g.PSAL_ADJUSTED_glodap)) - 0.1
    tmax = float(max(g.TEMP_ADJUSTED_glodap)) + 0.1
    tmin = float(min(g.TEMP_ADJUSTED_glodap)) - 0.1
    # Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(tmin,tmax,num=int(round((tmax-tmin)*10 + 1, 0)))
    si = np.linspace(smin,smax,num=int(round((smax-smin)*10 + 1, 0)))
    dens = np.zeros((len(ti),len(si)))
    # Loop to fill in grid with densities
    for j in range(0,len(ti)):
        for i in range(0, len(si)):
            dens[j,i]=gsw.sigma0(si[i],ti[j])

    axn = plt.subplot(4,4,16)
    CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
    axn.clabel(CS, fontsize=12, inline=1)
    cn = axn.scatter(g.PSAL_ADJUSTED_glodap,g.TEMP_ADJUSTED_glodap,c=g.DOXY_ADJUSTED_offset)
    axn.set_title('GDAP T-S')
    cx = fig.add_axes([0.72,0.12,0.02,0.17])
    plt.colorbar(cn,cax=cx,label='DOXY OFFSET')
    plt.tight_layout()
    plt.savefig(individual_plot_dir + str(g.main_float_wmo.values[0])+'_v_glodap_v2.png')
    plt.clf()
#     if show_plot is False:
#         plt.close()
    # wait = input("Press Enter to continue.")
    
# -

# Convert the list to an xarray DataArray or Dataset
glodap_drift_possible_dataarray = xr.DataArray(glodap_drift_possible_list, dims=['main_float_wmo'])


# ### 2. Combine metadata (calibration, instrument etc) with offsets so offsets can be sorted and plotted using different criteria
#
# For the following analyis we only need the mean offset for each float, so can take the mean by each float wmo first.

glodap_drift_possible_dataarray

glodap_offsets_mean

# +
glodap_offsets_mean = offsets_g.mean(skipna='True')

# add in boolean array of whether drift/real ocean change can explain glodap/float differences
# glodap_offsets_mean['glodap_drift_possible'] = glodap_drift_possible_dataarray
# -

# We also need to create broad groupings of calibration information provided in SCIENTIFIC_CALIB_COMMENT to enable easy sorting and plotting

# +
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
# -

# Now we can loop through each WMO and combine the calibration + sensor information with the mean offset for each float (note that in some cases individual profiles from a float may have different calibration comments and coeffients, but each float should be able to be grouped by O2 calibration)

# +
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
#save offsets with cal info
glodap_offsets_mean.to_netcdf(output_dir+glodap_offsets_filename[0:-3]+ '_floatmean_withcalibration.nc')


# +
# copying crossovers for floats with both air cal and pH to another directory
glodap_offsets_p = glodap_offsets_mean.to_dataframe()

parameter_a = 'o2_calib_air_group'
parameter_b = 'pH_group'
# offsets_g = glodap_offsets_p.groupby(parameter_a)
offsets_pH = glodap_offsets_p.groupby([parameter_a, parameter_b])
# -

air_cal_ph = offsets_pH.get_group(('air cal', 'pH'))
# print(air_cal_ph.index)
for wmo in air_cal_ph.index:
    print(wmo)
    os.system('cp ' + offset_dir+'individual_floats/' + str(wmo) + '_v_glodap_v2.png ' + \
              offset_dir+'individual_floats/air_cal_w_ph/' + str(wmo) + '_v_glodap_v2.png')  

os.system('cp ' + offset_dir+'individual_floats/' + str(wmo) + '_v_glodap_v2.png ' + \
              offset_dir+'individual_floats/air_cal_w_ph/' + str(wmo) + '_v_glodap_v2.png')  


# +
# print numbers in different calibration group
cal_groups = np.unique(glodap_offsets_mean['o2_calib_group'])
print(cal_groups)
for cg in cal_groups:
    
    print(cg + ' '+ str(np.sum(glodap_offsets_mean['o2_calib_group']==cg).values))
    

glodap_offsets_mean['o2_calib_comment'][glodap_offsets_mean['o2_calib_group']=='ungrouped']
# -

np.sum(glodap_offsets_mean['o2_calib_group']==cg).values

# +
air_subset_1 = glodap_offsets_mean['o2_calib_comment'][glodap_offsets_mean['o2_calib_comment'].str.contains('in air')]

print(air_subset_1)

# air_subset_1[glodap_offsets_mean['o2_calib_air_group'][glodap_offsets_mean['o2_calib_comment'].str.contains('no in air')] \
#     .str.contains('no')]
# -

sum(glodap_offsets_mean['o2_calib_air_group'].str.match("no air cal"))

# ### 3. Example: group offsets by chosen float variables/metadata and plot histograms
#
# In this example, float offsets are grouped by a) floats with/without pH and b) DOXY air calibrated/not air calibrated. Parameters a) and b) can be substituted by any other dataarray in glodap_offsets_mean
#
# First define parameters a) and b) to group offsets

glodap_offsets_mean


parameter_a = 'pH_group'
parameter_b = 'o2_calib_air_group'
#can pandas group by three different things?

# Now group mean offsets by parameters a and b and print statistics

#convert to pandas
glodap_offsets_p = glodap_offsets_mean.to_dataframe()
#do groupby paramater a and b
offsets_g = glodap_offsets_p.groupby([parameter_a,parameter_b])

# Plot pdfs for each group

glodap_offsets_mean

n

# +
#loop through groups and plot histograms and mean/median values of offsets
plt.figure(figsize=(16,10))
for n,group in offsets_g:
    #if n[1] == 'no cal/bad':
    #    continue
    print(n)
    
    #calc mean values
    nmean = np.around(group['DOXY_ADJUSTED_offset'].mean(), decimals=2)
    nmedian = np.around(group['DOXY_ADJUSTED_offset'].median(), decimals=2)
    ncount = group['DOXY_ADJUSTED_offset'].count()
    print(nmean)
    print(nmedian)
    print(ncount)

    #plot histogram
    plt.hist(group['DOXY_ADJUSTED_offset'], bins=np.linspace(-60, 60, 121),
             alpha=0.3,label=str(n)+', mean='+str(nmean) + ', n='+str(ncount))



plt.xlabel('DOXY Offset')
plt.legend()
plt.savefig(offset_dir + 'Glodap_offsets_doxy_all_'+parameter_a+parameter_b+'_v3.png')
# -

# ### 4. Example: map offsets for single float sub-group
#
# In this example plot histogram of float-mean DOXY_ADJUSTED offsets for air-calibrated floats only

parameter = 'o2_calib_air_group'
group_name = 'air cal'

# group offsets by single parameter

offsets_g1 = glodap_offsets_mean.groupby(parameter)

# now map offsets for single group, in this example floats with air calibration only - this should be profile by profile offsets not float by float??

for n, group in offsets_g1:
    if n == group_name:
        fig = plt.figure(figsize=(20,12))
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_global()
        sct = plt.scatter(group.main_float_longitude, 
                group.main_float_latitude, 
                c=group.DOXY_ADJUSTED_offset,cmap='RdBu_r',s=40,vmin=-10,vmax=10)
        cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
        cbar.set_label('O2 offset', labelpad=15, fontsize=14)
        #plt.scatter(glodap_offsets.glodap_longitude,glodap_offsets.glodap_longitude,s=4)
        plt.savefig(offset_dir+ 'map_o2_offsets_air_only_v3.png')
        plt.show()

# ### 5. Example: histogram of all global DOXY offsets, mean of each float

plt.figure()
plt.hist(glodap_offsets_mean['DOXY_ADJUSTED_offset'], bins=np.linspace(-20, 20, 61),label=str(n))
print(np.around(glodap_offsets_mean['DOXY_ADJUSTED_offset'].median().values, decimals=2))
plt.title('All Offsets, first averaged by float')

# ### Other code

plt.figure(figsize=(20,12))
plt.hist(glodap_offsets_mean.DOXY_ADJUSTED_offset, bins=np.linspace(-400, 400, 401))
plt.xlabel('DOXY Offset')
plt.savefig(offset_dir + 'Glodap_offsets_doxy_plus_minus_400.png')


#check location of O2 offsets > 100
offset_100 = glodap_offsets.DOXY_ADJUSTED_offset.where(glodap_offsets.DOXY_ADJUSTED_offset<-100)


fig = plt.figure(figsize=(20,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_global()
sct = plt.scatter(glodap_offsets.main_float_longitude, 
                glodap_offsets.main_float_latitude, 
                c=offset_100,cmap='viridis',s=10,vmin=-300,vmax=-100)
cbar = plt.colorbar(sct, fraction=.08, pad = 0.04, shrink=0.5)
cbar.set_label('O2 offset', labelpad=15, fontsize=14)
plt.scatter(glodap_offsets.glodap_longitude,glodap_offsets.glodap_longitude,s=4)
plt.savefig(offset_dir+ 'map_o2_offsets_100.png')
plt.show()

# +
#nitrate all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.NITRATE_ADJUSTED_offset, bins=np.linspace(-4, 4, 30))
nmean = np.around(glodap_offsets['NITRATE_ADJUSTED_offset'].mean().values, decimals=2)
nmedian = np.around(glodap_offsets['NITRATE_ADJUSTED_offset'].median().values, decimals=2)
plt.grid()
    
plt.title('Mean: ' + str(nmean))
plt.xlabel('NITRATE Offset')
plt.savefig(offset_dir + 'Glodap_offsets_nitrate_plus_minus_400.png')


# +
#DIC all glodap offsets
plt.figure(figsize=(20,12))
plt.hist(glodap_offsets.DIC_offset, bins=np.linspace(-40, 40, 30))
nmean = np.around(glodap_offsets['DIC_offset'].mean().values, decimals=2)
nmedian = np.around(glodap_offsets['DIC_offset'].median().values, decimals=2)
plt.grid()
    
plt.title('Mean: ' + str(nmean))
plt.xlabel('DIC Offset')
plt.savefig(offset_dir + 'Glodap_offsets_DIC.png')
