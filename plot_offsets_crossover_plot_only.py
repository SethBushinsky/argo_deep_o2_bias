import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.dates as mdates
import numpy as np




def plot_glodap_crossovers(glodap_offsets, offsets_g):
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
    #     break
    #     if show_plot is False:
    #         plt.close()
        # wait = input("Press Enter to continue.")