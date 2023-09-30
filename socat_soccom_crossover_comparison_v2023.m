% SOCAT SOCCOM crossover comparison v 2023


plot_dir = [home_dir 'Work/Presentations/2023_09 BGC Argo Science Meeting/plots/'];

SOCCOM_dir = [data_dir 'ARGO_O2_Floats/Global/SOCCOM/2022_05_19_Snapshot_LoRes_LIAR/'];
load([SOCCOM_dir 'SO_calc_28-Jun-2022_w_calc_param_pco2_Global_SOCCOM_only.mat'])



%%
% still need to update socat to 2023
socat = load([data_dir 'Data_Products/SOCAT/v2022/SOCATv2022_Global_ABCD_WOCE_2.mat']);
%% socat calculations

socat.calc_pCO2 = pCO2_from_fCO2(socat.fCO2rec, socat.SST);

socat.GMT_Matlab = datenum(socat.yr, socat.mon, socat.day, socat.hh, 0, 0);
socat.pdens = sw_pden(socat.sal, socat.SST, 0, 0);


%%
%% pre-calculate SOCAT distance
Argo_socat_dist = [];
day_range = 1;
distance_range = 25; % km
explore_plots=0;

for f =  1:length(SO_SNs)
    
    if ~isfield(Argo.(SO_SNs{f}), 'pCO2_LIAR')
        continue
    end
    if sum(sum(~isnan(Argo.(SO_SNs{f}).pCO2_LIAR(:,1:5))))==0
        %         disp('blue')
        continue
    end
    %     Argo_socat_dist.(SO_SNs{f}).socat_dist_all = NaN(length(Argo.(SO_SNs{f}).GMT_Matlab), length(socat.latitude));
    Argo_socat_dist.(SO_SNs{f}).profile_time_match = [];
    %     Argo_socat_dist.(SO_SNs{f}).socat_dist_all = [];
    
    Argo_socat_dist.(SO_SNs{f}).mean_dist = [];
    Argo_socat_dist.(SO_SNs{f}).mean_time = [];
    Argo_socat_dist.(SO_SNs{f}).mean_pCO2 = [];
    Argo_socat_dist.(SO_SNs{f}).mean_sst = [];
    Argo_socat_dist.(SO_SNs{f}).mean_sal = [];
    Argo_socat_dist.(SO_SNs{f}).mean_lat = [];
    Argo_socat_dist.(SO_SNs{f}).mean_lon = [];
    Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr = [];
    Argo_socat_dist.(SO_SNs{f}).expocode = {};
    
    for p=1:length(Argo.(SO_SNs{f}).GMT_Matlab)
        
        time_match = (socat.GMT_Matlab>= Argo.(SO_SNs{f}).GMT_Matlab(p)-day_range)  & (socat.GMT_Matlab<=  Argo.(SO_SNs{f}).GMT_Matlab(p)+day_range);
        if sum(time_match)>0
            lon_temp_argo = Argo.(SO_SNs{f}).Lon(p);
            lon_temp_argo(lon_temp_argo>180) = lon_temp_argo(lon_temp_argo>180)-360;
            
            lon_temp_socat = socat.longitude(time_match);
            lon_temp_socat(lon_temp_socat>180) = lon_temp_socat(lon_temp_socat>180)-360;
            
            socat_dist = distance('gc',Argo.(SO_SNs{f}).Lat(p), lon_temp_argo, socat.latitude(time_match), lon_temp_socat).*111.19; % distance in km
            
            if sum(socat_dist<=distance_range)
%                 disp(p)
                
                
                Argo_socat_dist.(SO_SNs{f}).mean_dist(end+1,1) = nanmean(socat_dist(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_dist(end,2) = nanstd(socat_dist(socat_dist<distance_range));
                
                temp_time = socat.GMT_Matlab(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_time(end+1,1) = nanmean(temp_time(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_time(end,2) = nanstd(temp_time(socat_dist<distance_range));
                
                temp_pCO2 = socat.calc_pCO2(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end+1,1) = nanmean(temp_pCO2(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end,2) = nanstd(temp_pCO2(socat_dist<distance_range));
                
                temp_sst = socat.SST(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_sst(end+1,1) = nanmean(temp_sst(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_sst(end,2) = nanstd(temp_sst(socat_dist<distance_range));
                
                temp_sal = socat.sal(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_sal(end+1,1) = nanmean(temp_sal(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_sal(end,2) = nanstd(temp_sal(socat_dist<distance_range));
                
                temp_lat = socat.latitude(time_match);
                Argo_socat_dist.(SO_SNs{f}).mean_lat(end+1,1) = nanmean(temp_lat(socat_dist<distance_range));
                Argo_socat_dist.(SO_SNs{f}).mean_lat(end,2) = nanstd(temp_lat(socat_dist<distance_range));
                
                temp_lon = socat.longitude(time_match);
                filt_temp_lon = temp_lon(socat_dist<distance_range);
                
                if (max(filt_temp_lon) - min(filt_temp_lon)) > 300
                    filt_temp_lon(filt_temp_lon<100) = filt_temp_lon(filt_temp_lon<100) + 360;
                end
                filt_temp_lon_mean = nanmean(filt_temp_lon);
                if filt_temp_lon_mean>360
                    filt_temp_lon_mean = filt_temp_lon_mean-360;
                end
                Argo_socat_dist.(SO_SNs{f}).mean_lon(end+1,1) = nanmean(filt_temp_lon_mean);
                Argo_socat_dist.(SO_SNs{f}).mean_lon(end,2) = nanstd(filt_temp_lon);
                
                
                % recalculate the SOCAT pCO2 to account for the temperature
                % difference with the float
                ArgALK_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                ArgSST_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                Press_CO2sys = NaN(length(temp_pCO2(socat_dist<distance_range)),1);
                Press_CO2sys(:) = 5;
                
                ArgALK_CO2sys(:) = Argo.(SO_SNs{f}).TALK_LIAR(p,1);
                ArgSST_CO2sys(:) = Argo.(SO_SNs{f}).Temp_C(p,1);
                
                [DATA,~,~]= CO2SYSSOCCOM(ArgALK_CO2sys, temp_pCO2(socat_dist<distance_range) , ...
                    1,4,temp_sal(socat_dist<distance_range), temp_sst(socat_dist<distance_range), ...
                    ArgSST_CO2sys,...
                    Press_CO2sys,Press_CO2sys,20,1.8,1,10,3);
                
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr(end+1,1) = nanmean(DATA(:,19));
                Argo_socat_dist.(SO_SNs{f}).mean_pCO2_sstcorr(end,2) = nanstd(DATA(:,19));
                
                Argo_socat_dist.(SO_SNs{f}).profile_time_match(end+1) = p;
                
                temp_expocode = socat.Expocode(time_match);
                temp_expocode = temp_expocode(socat_dist<distance_range);
                
                Argo_socat_dist.(SO_SNs{f}).expocode{end+1,1} =  unique(temp_expocode);
                
                %                 Argo_socat_dist.(SO_SNs{f}).socat_dist_all(end+1,:) = NaN(1, length(socat.latitude));
                %                 Argo_socat_dist.(SO_SNs{f}).socat_dist_all(end,time_match) = socat_dist;
                if explore_plots==1
                    clf
                    subplot(4,2,1); plot(socat_dist(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_dist(end,1), 'ro');
                    
                    subplot(4,2,2); plot(temp_time(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_time(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).GMT_Matlab(p), 'm^')
                    
                    subplot(4,2,3); plot(temp_pCO2(socat_dist<distance_range+50), 'gs'); hold on
                    plot(DATA(:,19), 'x'); hold on;
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_pCO2(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).pCO2_LIAR(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,4); plot(temp_sst(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_sst(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).Temp_C(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,5); plot(temp_sal(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_sal(end,1), 'ro');
                    plot(Argo.(SO_SNs{f}).Sal(p, Argo.(match_SNs{f}).Press_db(p,:)<=10), 'm^')
                    
                    subplot(4,2,6); plot(temp_lat(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_lat(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).Lat(p), 'm^')
                    
                    subplot(4,2,7); plot(temp_lon(socat_dist<distance_range), 'x'); hold on
                    plot(Argo_socat_dist.(SO_SNs{f}).mean_lon(end,1), 'ro');
                    
                    plot(Argo.(SO_SNs{f}).Lon(p), 'm^')
                    
                    pause
                end
                
            end
        end
    end
    

    
    disp(f)
    %         save([data_dir 'Argo_SOCAT_dist_' num2str(f) '.mat'], 'Argo_socat_dist')
    
end

clear explore_plots temp_expocode DATA ArgALK_CO2sys Press_CO2sys temp_pCO2 temp_sal temp_sst
clear temp_lon temp_lat filt_temp_lon filt_temp_lon_mean lon_temp_argo lon_temp_socat
clear clear socat_dist f time_match ArgSST_CO2sys 
%% creating the difference table
match_SNs = fieldnames(Argo_socat_dist);
diff_table = NaN(150,17);
diff_index = 0;

pCO2_type = 'pCO2_LIAR';

for f=1:length(match_SNs)
    for q = 1:length(Argo_socat_dist.(match_SNs{f}).profile_time_match)
        diff_index = diff_index+1;
        matched_profiles = Argo_socat_dist.(match_SNs{f}).profile_time_match;
        %         disp(f)
        %         disp(q)
        argo_match_depths = Argo.(match_SNs{f}).Press_db(matched_profiles(q),:)<=10;
        %                 argo_match_depths = 1:4;
        
        %         plot(Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,1), nanmean(Argo.(match_SNs{f}).pCO2_LIAR(matched_profiles(q),argo_match_depths)), 's')
        %         diff_table(diff_index,1) = Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,1);
        diff_table(diff_index,1) = Argo_socat_dist.(match_SNs{f}).mean_pCO2_sstcorr(q,1);
        
        diff_table(diff_index,2) =  nanmean(Argo.(match_SNs{f}).(pCO2_type)(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,3) = Argo_socat_dist.(match_SNs{f}).mean_sst(q,1);
        diff_table(diff_index,4) = nanmean(Argo.(match_SNs{f}).Temp_C(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,5) = Argo_socat_dist.(match_SNs{f}).mean_sal(q,1);
        diff_table(diff_index,6) = nanmean(Argo.(match_SNs{f}).Sal(matched_profiles(q),argo_match_depths));
        
        diff_table(diff_index,7) = Argo_socat_dist.(match_SNs{f}).mean_dist(q,1);
        
        diff_table(diff_index,8) = Argo_socat_dist.(match_SNs{f}).mean_time(q,1);
        diff_table(diff_index,9) = (Argo.(match_SNs{f}).GMT_Matlab(matched_profiles(q),1));
        
        diff_table(diff_index,10) = Argo_socat_dist.(match_SNs{f}).mean_lat(q,1);
        diff_table(diff_index,11) = (Argo.(match_SNs{f}).Lat(matched_profiles(q),1));
        
        diff_table(diff_index,12) = Argo_socat_dist.(match_SNs{f}).mean_lon(q,1);
        diff_table(diff_index,13) = (Argo.(match_SNs{f}).Lon(matched_profiles(q),1));
        
        if abs(diff_table(diff_index,12) - diff_table(diff_index,13))>10  && abs(diff_table(diff_index,12) - diff_table(diff_index,13))<355
            disp(f); disp(diff_index); disp(' ');
        end
        diff_table(diff_index,14) = Argo_socat_dist.(match_SNs{f}).mean_time(q,1) - Argo.(match_SNs{f}).GMT_Matlab(matched_profiles(q),1);
        
        diff_table(diff_index,15) = Argo_socat_dist.(match_SNs{f}).mean_pCO2(q,2);
        diff_table(diff_index,16) = nanstd(Argo.(match_SNs{f}).(pCO2_type)(matched_profiles(q),argo_match_depths));
        
        % float deployment
        diff_table(diff_index,16) = Argo.(match_SNs{f}).GMT_Matlab(1);
        diff_table(diff_index,17) = str2double(match_SNs{f}(2:end));
        
    end
end

clear mathed_profiles argo_match_depths f q 
%% Figure S3 
density_threshold = 0.05;


plot_filename = ['Fig_S3_SOCAT_Float_1_to_1_SOCATv2022_AB ' num2str(day_range) ' days and '  num2str(distance_range) ' km SST Corr ' pCO2_type ' density_comp ' num2str(density_threshold)];
plasma_map = plasma;

clf
set(gcf, 'units', 'inches')
paper_w = 12; paper_h =5;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);

cmap = brewermap(8, 'Paired');
set(gcf, 'colormap', plasma_map)
d1 = subplot(1,2,1);
hold on
grid on

socat_minus_argo = diff_table(:,1) - diff_table(:,2);


plot([250 440], [250 440], '-k')
plot([250 440], [240 430], '--k')
plot([250 440], [260 450], '--k')

xlabel('SOCAT pCO_2 (\muatm)')
ylabel('Float pCO_2 (\muatm)')
set(gca, 'ylim', [250 440], 'xlim', [250 440])
d2 = subplot(1,2,2);

hold on
plot([0 0], get(gca, 'ylim'), '--k')


socat_dens = sw_dens(diff_table(:,5), diff_table(:,3), 1);
float_dens = sw_dens(diff_table(:,6), diff_table(:,4), 1);
dens_diff = abs(socat_dens - float_dens);
%
plot(d1, diff_table(dens_diff<=density_threshold,1), diff_table(dens_diff<=density_threshold,2), ...
    'marker', 's','markeredgecolor', 'k', 'markerfacecolor', 'b', 'linestyle', 'none')
hold on
title(d1, 'SOCCOM and SOCAT Matchup comparison')
% plot(d1, diff_table(dens_diff>density_threshold,1), diff_table(dens_diff>density_threshold,2), 'ro', 'markersize', 15)
diff_filtered_dens_difference = diff_table(:,1) - diff_table(:,2);
diff_filtered_dens_difference(dens_diff>density_threshold)=nan;
%
cla(gca)

hold on
h2 = histogram(diff_filtered_dens_difference);


h2.BinWidth = 5;
h2.FaceColor = cmap(4,:);
% set(gca, 'ylim', [0 20])
plot([0 0], get(gca, 'ylim'), '--k')
orig_x_lim = get(d2, 'xlim');
grid on
set(d2, 'xlim', [-max(abs(orig_x_lim)) max(abs(orig_x_lim))])
xlabel('SOCAT minus Argo (\muatm)')
title('Differences')
disp(['Green: remove > ' num2str(density_threshold) 'kg/m3 diff. Mean filtered: ' ...
    num2str(nanmean(diff_filtered_dens_difference),3) ', std: ' num2str(nanstd(diff_filtered_dens_difference),3)])
disp( ['original n: ' num2str(sum(~isnan(socat_minus_argo))) ' ' ...
    'filtered n: ' num2str(sum(~isnan(diff_filtered_dens_difference)))])

print(gcf,'-dpng', '-r300', [plot_dir plot_filename '.png'])


clear d2 h2 orig_x_lim d1

% map of offsets:

plot_filename = ['Map_SOCAT_Float_1_to_1_SOCATv2022_AB ' num2str(day_range) ' days and '  num2str(distance_range) ' km SST Corr ' pCO2_type ' density_comp ' num2str(density_threshold)];

c_map = brewermap(20, 'RdBu');
set(gcf, 'colormap', flipud(c_map))
clf
set(gcf, 'units', 'inches')
paper_w = 12; paper_h =10;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h

d1 =        subplot(1,1,1);

m1= axesm ('pcarree', 'Frame', 'off', 'Grid', 'on','Origin',[0,0], 'MapLatLimit',[-80 80], 'MapLonLimit', [-180 180]);
plabel on; mlabel on; hold on
%         setm(d1)
geoshow(landareas,'FaceColor',[.6 .6 .6],'EdgeColor',[.2 .2 .2]);
setm(m1, 'MlabelParallel', 'south');

scatterm(d1, diff_table(:,11), diff_table(:,13), 60, ...
   diff_filtered_dens_difference, 'filled');  %clim(limits(s,:))
c1 = colorbar;
clim([-30 30])

%         geoshow([-75 0], [lon_line lon_line], 'displaytype', 'line', 'linestyle', '--', 'color', [.3 .3 .3])
ylabel(c1, 'pCO2 diff (\muatm)')
title(['SOCAT minus Argo'])
m1.Box = 'off';

print(gcf,'-dpng', '-r300', [plot_dir plot_filename 'density_comp ' num2str(density_threshold) '.png'])
