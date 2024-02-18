import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
import numpy as np
import functions.sokal_rohlf_calculations as SR
import pandas as pd
from scipy import interpolate
import scipy.stats as stats
import gsw
from datetime import datetime


# will want to pass individual float files for plotting, so that they can be distributed to parallel cores
def plot_glodap_crossovers(individual_plot_dir, mean_level_offsets, g, pressure_level_min_max, float_age_bins):
    print('Plotting crossover for: ' + str(g.main_float_wmo.values[0]))
    # show_plot = False
    #loop through each float
    label_size = 15
    fig = plt.figure(figsize=(26,26))
    plt.rcParams.update({'font.size': 15})  

    # create a list to put in True / False if the glodap offsets intersect zero at float mid date
    glodap_drift_possible_list = []

    ncross = len(g.DOXY_ADJUSTED_offset)

    #if doxy_adjusted_offset is nan, skip
    if np.all(np.isnan(g.DOXY_ADJUSTED_offset.values)):
        print('No non-nan DOXY offsets')
        return
    # print('At step 1: ' + str(g.main_float_wmo.values[0]))

    plt_rows = 4
    plot_cols = 6

    #add crossover location map
    axn = plt.subplot(plt_rows,4,1)
    axn.plot(g.main_float_longitude,g.main_float_latitude,'g+',markersize=10,label='Current float')
    #glodap
    glodap_lon = xr.where(g.glodap_longitude>180,g.glodap_longitude-360.,g.glodap_longitude)
    axn.plot(glodap_lon,g.glodap_latitude,'ro',label = 'Glodap',markersize=10)
    axn.legend()
    axn.set_title('WMO: %s; N: %d; Dist filt %d; Depth min %d max %d\nDepth filt %d; Dens %.4f; Spice %.4f; Plot depth range: %d to %d' % 
                (g.main_float_wmo.values[0],ncross, g["dist"][0],
                g["p_compare_min"][0], 
                g["p_compare_max"][0],g["delta_press"][0],
                g["delta_dens"][0],g["delta_spice"][0], pressure_level_min_max[0], pressure_level_min_max[1]), fontsize=label_size)
    plt.ylabel('Latitude', fontsize=label_size)
    plt.xlabel('Longitude', fontsize=label_size)

    #time
    axn = plt.subplot(plt_rows,10,(4,10))
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset,'rx',label='float')
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset_trimmed,'ro',label='float_outliers_removed', markerfacecolor="None")
    axn.plot(g.glodap_datetime,g.DOXY_ADJUSTED_offset,'bx',label='GDAP')
    axn.plot(g.glodap_datetime,g.DOXY_ADJUSTED_offset_trimmed,'bo', markerfacecolor="None")

    # fit regressions
    g_ox = g.where(~np.isnan(g.DOXY_ADJUSTED_offset_trimmed), drop=True)
    if g_ox.dims['N_CROSSOVERS']!=0:
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
        b_yx, a, r2, CI_alpha_slope, ttt, y_err = SR.regress_confidence_sokal_rohlf(X_series, Y_series, CI_alpha)

    #     t = np.arange(datetime(1980,1,1), datetime(2022,1,1), timedelta(days=1)).astype(datetime)
    #     tnum = mdates.date2num(t)
        print(a)
        print(b_yx)
        if not np.all(np.isnan(a)):
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
                X_extend = pd.concat([X_series, pd.Series(float_mid_date)])
                y_extend = a+ X_extend*b_yx

                X_series_last_third = X_series.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
                Y_err_last_third = y_err.iloc[np.int64(np.round(len(X_series)*1/2)):-1]
                b_err, a_err, _, _, _, _ = SR.regress_confidence_sokal_rohlf(X_series_last_third, Y_err_last_third, CI_alpha)
                extrap_error = X_extend.iloc[-1]*b_err+a_err
                y_err_extend = pd.concat([y_err, pd.Series(extrap_error)])

        #         y_err_extend = y_err.append(pd.Series(y_err.iloc[-1]))

                plt.plot(X_extend, y_extend, color='blue', linestyle = '--')
                plt.fill_between(X_extend, y_extend-y_err_extend, y_extend+y_err_extend, color='blue', alpha=0.25)    
                uncert_min = y_extend.iloc[-1] - y_err_extend.iloc[-1]
        #         print(uncert_min)
                uncert_max = y_extend.iloc[-1] + y_err_extend.iloc[-1]
        #         print(uncert_max)
            elif np.min(X_series)>float_mid_date:
                X_extend = pd.concat([pd.Series(float_mid_date), X_series])
                y_extend = a+ X_extend*b_yx
                X_series_first_half = X_series.iloc[0:np.int64(np.round(len(X_series)*1/2))]
                Y_err_first_half = y_err.iloc[1:np.int64(np.round(len(X_series)*1/2))]
                b_err, a_err, _, _, _, _ = SR.regress_confidence_sokal_rohlf(X_series_first_half, Y_err_first_half, CI_alpha)
                extrap_error = X_extend.iloc[0]*b_err+a_err
                y_err_extend = pd.concat([pd.Series(extrap_error), y_err])

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
        else:
            uncert_min = -100
            uncert_max = 100

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

    axn.set_title('O2 Offsets vs GLODAP date', fontsize=label_size) #, Real ocean O2 change? ' + str(glodap_drift_possible))
    #     axn.text(0.01, 0.01, "y = "+ str(m1.round(decimals=4)) +"x + " + str(b1.round(decimals=2)) + "\n"
    #              "y = "+ str(m2.round(decimals=4)) +"x + " + str(b2.round(decimals=2)),
    #              verticalalignment='bottom', horizontalalignment='left',
    #              transform=axn.transAxes)
    plt.ylabel(r'$\Delta$O$_{2,off}$  ($\mu$mol kg$^{-1})$', fontsize=label_size)

    # Print meta data
    axn = plt.subplot(plt_rows,4,5)
    plt.text(0,3, 'Project: ' + str(mean_level_offsets.sel(main_float_wmo=g.main_float_wmo.values[0]).project_name.values))
    plt.text(0,2, 'O2 Sensor: ' + str(mean_level_offsets.sel(main_float_wmo=g.main_float_wmo.values[0]).DOXY_sensor.values))
    plt.text(0,1, 'Data Center: ' + str(mean_level_offsets.sel(main_float_wmo=g.main_float_wmo.values[0]).data_centre.values))
    plt.text(0,0, 'Cal. Group: ' + str(mean_level_offsets.sel(main_float_wmo=g.main_float_wmo.values[0]).o2_calib_air_group.values))
    axn.axis('off')
    axn.set_ylim([-1, 4])


    # Float time series
    axn = plt.subplot(plt_rows,10,(14,20))
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset,'rx',label='float')
    axn.plot(g.main_float_juld,g.DOXY_ADJUSTED_offset_trimmed,'ro',label='float_outliers_removed', markerfacecolor="None")

    numeric_timestamp = g.First_Float_Profile_Date[0].values

    time_binned_offset = np.empty((len(float_age_bins)-1, 2), dtype=object)
    for fa in range(len(float_age_bins)-1):
        timestamp_seconds = numeric_timestamp / 1e9 + ((float_age_bins[fa]+float_age_bins[fa+1])/2)*24*60*60
        datetime_object = datetime.utcfromtimestamp(timestamp_seconds)
        time_binned_offset[fa] = (datetime_object, g['DOXY_ADJUSTED_offset_' + 'age_' + str(float_age_bins[fa])].values[0])

    axn.axhline(y=0, color='k', linestyle='--')
    axn.plot(time_binned_offset[:, 0], time_binned_offset[:, 1],'k-s',label='float', linewidth=2)
    axn.set_title('Offsets vs. float time only with time-binned means overlaid (black squares)', fontsize=label_size)
    plt.ylabel(r'$\Delta$O$_{2,off}$  ($\mu$mol kg$^{-1})$', fontsize=label_size)
    plt.tight_layout()

    hist_col_index = 13


    # plot offset histograms
    # Temperature
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index)
    g_plot = g.TEMP_ADJUSTED_offset
    if not np.all(np.isnan(g_plot.values)):
        g_mean = np.nanmean(g_plot.values).round(decimals=2)
        g_std = np.nanstd(g_plot.values).round(decimals=2)
        axn.hist(g_plot,color='b',alpha=0.5)
        axn.set_title('TEMP %.2f +/- %.2f' % (g_mean,g_std), fontsize=label_size)
    plt.grid()
    plt.xlabel(r'$\Delta$$\theta$$_{off}$ ($^{\circ}$ C)', fontsize=label_size)
    plt.ylabel('Count')


    #Salinity
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index+1)
    g_plot = g.PSAL_ADJUSTED_offset
    if not np.all(np.isnan(g_plot.values)):
        g_mean = np.nanmean(g_plot.values).round(decimals=2)
        g_std = np.nanstd(g_plot.values).round(decimals=2)
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('PSAL %.2f +/- %.2f' % (g_mean,g_std), fontsize=label_size)
    plt.xlabel(r'$\Delta$S$_{off}$ ($^{\circ}$ C)', fontsize=label_size)

    plt.grid()

    # Oxygen histogram
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index+2)
    g_plot = g.DOXY_ADJUSTED_offset
    if not np.all(np.isnan(g_plot.values)):
        g_mean = np.nanmean(g_plot.values)
        t_stat, p_value = stats.ttest_1samp(a=g_plot, popmean=0, nan_policy='omit') ############

        g_std = np.nanstd(g_plot.values).round(decimals=2)

        g_plot_trim = g.DOXY_ADJUSTED_offset_trimmed
        if not np.all(np.isnan(g_plot_trim.values)):
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
            if not np.all(np.isnan(g_plot_trim.values)):
                axn.hist(g_plot_trim, bins=bins,color='r',alpha=0.5, label='trimmed' + ' p: ' + str(p_value_trim.round(5)))

            axn.legend()

            
        axn.set_title('O2 n: %.1f, %.1f +/- %.1f' % (len(g_plot[~np.isnan(g_plot)]),g_mean,g_std), fontsize=label_size)
    plt.xlabel(r'$\Delta$O$_{2,off}$ ($\mu$mol kg$^{-1}$)', fontsize=label_size)

    # axn.text(0.99, 0.99, "p_value: " + str(p_value.round(5)) + "\n"
    #      + "t-statistic: " + str(t_stat.round(3)),
    #      verticalalignment='top', horizontalalignment='right',
    #      transform=axn.transAxes)
    plt.grid()


    # oxygen histogram - trimmed data only:
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index+3)
    if not np.all(np.isnan(g_plot.values)):
        if not np.all(np.isnan(g_plot_trim)):
            axn.hist(g_plot_trim,color='r',alpha=0.5)
            plt.axvline(CI_99[1], color='k', linestyle='--', label='99% CI')
            plt.axvline(CI_99[0], color='k', linestyle='--')
            axn.legend()
            axn.set_title('O2 trim n: %.1f, %.1f +/- %.1f' % (len(g_plot_trim[~np.isnan(g_plot_trim)]), g_mean_trim, g_std_trim), fontsize=label_size)
    plt.grid()
    plt.xlabel(r'$\Delta$O$_{2,off}$ ($\mu$mol kg$^{-1}$)', fontsize=label_size)


    # Nitrate
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index+4)
    g_plot = g.NITRATE_ADJUSTED_offset
    if not np.all(np.isnan(g_plot)):
        g_mean = np.nanmean(g_plot.values).round(decimals=2)
        g_std = np.nanstd(g_plot.values).round(decimals=2)
        axn.hist(g_plot,color='b',alpha=0.5)
    axn.set_title('NO3 %.2f +/- %.2f' % (g_mean,g_std), fontsize=label_size)
    plt.grid()
    plt.xlabel(r'$\Delta$NO$_3$$^-$$_{off}$ ($\mu$mol kg$^{-1}$)', fontsize=label_size)


    # pH
    axn = plt.subplot(plt_rows,plot_cols,hist_col_index+5)
    axn.set_title('PH 25C %.2f +/- %.2f' % (g_mean,g_std))
    g_plot = g.pH_25C_TOTAL_ADJUSTED_offset
    if not np.all(np.isnan(g_plot)):
        g_mean = np.nanmean(g_plot.values*1000).round(decimals=2)
        g_std = np.nanstd(g_plot.values*1000).round(decimals=2)
        axn.hist(g_plot*1000,color='b',alpha=0.5)
    plt.grid()
    plt.xlabel(r'$\Delta$pH$_{off}$ (mpH)', fontsize=label_size)


    #O2 vs pressure
    axn = plt.subplot(plt_rows,4,13)
    axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_float,'bx')
    axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PRES_ADJUSTED_float,'bo', markerfacecolor="none")

    axn.set_title('Float pres vs DOXY offset', fontsize=label_size)
    plt.gca().invert_yaxis()
    axn.set_ylabel('Press. (db)')
    plt.xlabel(r'$\Delta$O$_{2,off}$ ($\mu$mol kg$^{-1}$)', fontsize=label_size)

    plt.grid()

    # axn = plt.subplot(plt_rows,4,15)
    # axn.plot(g.DOXY_ADJUSTED_offset,g.PRES_ADJUSTED_glodap,'bx')
    # axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PRES_ADJUSTED_glodap,'bo', markerfacecolor="none")

    # axn.set_title('GDAP pres vs DOXY offset')
    # plt.gca().invert_yaxis()
    # plt.grid()

    # #vs density
    # axn = plt.subplot(plt_rows,4,16)
    # axn.plot(g.DOXY_ADJUSTED_offset,g.PDENS_float,'bx')
    # axn.plot(g.DOXY_ADJUSTED_offset_trimmed, g.PDENS_float,'bo', markerfacecolor="none")

    # axn.set_title('Float Pdens vs DOXY offset')
    # axn.set_ylabel('dens')
    # plt.grid()


    #O2 offset vs o2 concentration
    # axn = plt.subplot(plt_rows,4,(plt_rows-1)*4+1)
    # #     axn.plot(g.DOXY_ADJUSTED_offset,g.DOXY_ADJUSTED_float,'bx')
    # # trying lots of things:
    # float_time_days = mdates.date2num(g.main_float_juld)
    # time_since_deploy = float_time_days - float_time_days[0]

    # #     cruise_nums = np.unique(g.glodap_cruise)
    # # #     cruise_num_scaled = np.zeros(len(g.glodap_cruise))
    # #     for idx, cr in enumerate(cruise_nums):
    # #         cruise_num_scaled[g.glodap_cruise==cruise_nums[idx]] = idx

    obs_inds = np.unique(g.glodap_obs_index)
    obs_inds_scaled = np.zeros(len(g.glodap_obs_index))
    for idx, cr in enumerate(obs_inds):
        obs_inds_scaled[g.glodap_obs_index==obs_inds[idx]] = idx
        
    # cn = axn.scatter(g.DOXY_ADJUSTED_offset_trimmed, g.DOXY_ADJUSTED_float,30, obs_inds_scaled, 
    #             marker='o', cmap='tab20c')
    # plt.colorbar(cn,label='?')

    # axn.set_title('Float O2 vs DOXY offset')
    # axn.set_ylabel('O2 concentration')
    # plt.grid()

    axn = plt.subplot(plt_rows,4,(plt_rows-1)*4+2)
    #     axn.plot(g.DOXY_ADJUSTED_offset,g.DOXY_ADJUSTED_float,'bx')
    # trying lots of things:
    float_time_days = mdates.date2num(g.main_float_juld)
    time_since_deploy = float_time_days - float_time_days[0]

    cn = axn.scatter(g.DOXY_ADJUSTED_offset_trimmed, g.DOXY_ADJUSTED_glodap,30, obs_inds_scaled, 
                marker='o', cmap='tab20c')
    plt.colorbar(cn,label='GLODAP Samples')

    axn.set_title('Glodap O2 vs DOXY offset', fontsize=label_size)
    axn.set_ylabel('O2 concentration')
    plt.xlabel(r'$\Delta$O$_{2,off}$ ($\mu$mol kg$^{-1}$)', fontsize=label_size)

    plt.grid()


    if not np.all(np.isnan(g.TEMP_ADJUSTED_float)):

        # Add density lines FLOAT
        ASAL_ADJUSTED_float = gsw.SA_from_SP(g.PSAL_ADJUSTED_float, g.PRES_ADJUSTED_float, g.main_float_longitude, g.main_float_latitude)
        CT_ADJUSTED_float = gsw.CT_from_t(ASAL_ADJUSTED_float, g.TEMP_ADJUSTED_float, g.PRES_ADJUSTED_float)
        smax = float(np.nanmax(ASAL_ADJUSTED_float)) + 0.1
        smin = float(np.nanmin(ASAL_ADJUSTED_float)) - 0.1
        tmax = float(np.nanmax(CT_ADJUSTED_float)) + 0.1
        tmin = float(np.nanmin(CT_ADJUSTED_float)) - 0.1
        # Create temp and salt vectors of appropiate dimensions
        ti = np.linspace(tmin,tmax,num=int(round((tmax-tmin)*10 + 1, 0)))
        si = np.linspace(smin,smax,num=int(round((smax-smin)*10 + 1, 0)))
        dens = np.zeros((len(ti),len(si)))
        # Loop to fill in grid with densities
        for j in range(0,len(ti)):
            for i in range(0, len(si)):
                dens[j,i]=gsw.sigma0(si[i],ti[j])

        axn = plt.subplot(plt_rows,4,(plt_rows-1)*4+3)
        CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
        axn.clabel(CS, fontsize=12, inline=1)
        axn.scatter(ASAL_ADJUSTED_float,CT_ADJUSTED_float,c=g.DOXY_ADJUSTED_offset)
        axn.set_xlim(np.nanmin(ASAL_ADJUSTED_float) - .1, np.nanmax(ASAL_ADJUSTED_float) + .1)
        axn.set_ylim(np.nanmin(CT_ADJUSTED_float) - 0.1, np.nanmax(CT_ADJUSTED_float) + 0.1)
        axn.set_title('Float T-S', fontsize=label_size)
        plt.ylabel(r'$\theta$', fontsize=label_size)
        plt.xlabel('S', fontsize=label_size)

        # Add density lines GLODAP
        ASAL_ADJUSTED_glodap = gsw.SA_from_SP(g.PSAL_ADJUSTED_glodap, g.PRES_ADJUSTED_glodap, g.glodap_longitude, g.glodap_latitude)
        CT_ADJUSTED_glodap = gsw.CT_from_t(ASAL_ADJUSTED_glodap, g.TEMP_ADJUSTED_glodap, g.PRES_ADJUSTED_glodap)

        smax = float(np.nanmax(ASAL_ADJUSTED_glodap)) + 0.1
        smin = float(np.nanmin(ASAL_ADJUSTED_glodap)) - 0.1
        tmax = float(np.nanmax(CT_ADJUSTED_glodap)) + 0.1
        tmin = float(np.nanmin(CT_ADJUSTED_glodap)) - 0.1
        # Create temp and salt vectors of appropiate dimensions
        ti = np.linspace(tmin,tmax,num=int(round((tmax-tmin)*10 + 1, 0)))
        si = np.linspace(smin,smax,num=int(round((smax-smin)*10 + 1, 0)))
        dens = np.zeros((len(ti),len(si)))
        # Loop to fill in grid with densities
        for j in range(0,len(ti)):
            for i in range(0, len(si)):
                dens[j,i]=gsw.sigma0(si[i],ti[j])

        axn = plt.subplot(plt_rows,4,(plt_rows-1)*4+4)
        CS = plt.contour(si,ti,dens, linestyles='dashed', colors='k')
        axn.clabel(CS, fontsize=12, inline=1)
        cn = axn.scatter(ASAL_ADJUSTED_glodap,CT_ADJUSTED_glodap,c=g.DOXY_ADJUSTED_offset)
        axn.set_title('GDAP T-S', fontsize=label_size)
        # cx = fig.add_axes([0.72,0.12,0.02,0.17])
        plt.colorbar(cn,label=r'$\Delta$O$_{2,off}$ ($\mu$mol kg$^{-1}$)')
        # plt.tight_layout()
        # print('Ready to save: ' + str(g.main_float_wmo.values[0]))
        plt.ylabel(r'$\theta$', fontsize=label_size)
        plt.xlabel('S', fontsize=label_size)
    plt.tight_layout()

    plt.savefig(individual_plot_dir + str(g.main_float_wmo.values[0])+'_v_glodap_v3.png')
    plt.close()
    #     break
    #     if show_plot is False:
    #         plt.close()
        # wait = input("Press Enter to continue.")