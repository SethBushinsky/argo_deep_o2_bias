% v3 - using the global BGC argo dataset
% v2 - switching to a density matching method, 2020_11_06 (Biden gains the
% lead!)

% Float crossover assessment: looking for comparisons to shipboard data as
% well as float-float comparisons

% SOCCOM_float_directory = [home_dir 'Data/ARGO_O2_Floats/Global/SOCCOM/2020_08_30_Snapshot_LoRes_LIAR/'];
% load([SOCCOM_float_directory 'SO_calc_14-Oct-2020_w_calc_param_o2.mat']);


% SOCCOM_float_directory = [home_dir 'Data/ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2020_10_23/'];
% SOCCOM_float_directory = [home_dir 'Data/ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2020_11_16/'];
% SOCCOM_float_directory = [data_dir 'ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2022_08_25/'];
SOCCOM_float_directory = [data_dir 'Data_Products/BGC_ARGO_GLOBAL/2022_08_25/'];

% load([SOCCOM_float_directory 'Argo_BGC_2020-10-26.mat']);
% load([SOCCOM_float_directory 'Argo_BGC_2020-11-17.mat']);
% load([SOCCOM_float_directory 'Argo_BGC_2022-09-06.mat']);
load([SOCCOM_float_directory 'Argo_BGC_2022-11-16.mat']);

%%
% gdap = load([home_dir 'Data/Data_Products/GLODAP/GLODAPv2.2020_Merged_Master_File.mat']);
gdap = load([data_dir 'Data_Products/GLODAP/GLODAPv2.2021_Merged_Master_File.mat']);

clear gdap_SO
gdap.PDENS = sw_pden(gdap.G2salinity,gdap.G2temperature,gdap.G2pressure, 0);
gdap.GMT_Matlab = datenum(gdap.G2year, gdap.G2month, gdap.G2day, gdap.G2hour, gdap.G2minute, zeros(size(gdap.G2minute)));


gdap_fields_no_flags = {'G2longitude', 'G2latitude', 'GMT_Matlab', 'G2year', 'G2month', 'G2cruise', 'G2station', 'G2pressure', 'G2temperature', 'PDENS'};
gdap_fields_w_flags = {'G2salinity','G2oxygen', 'G2nitrate', 'G2tco2', 'G2talk', 'G2phts25p0'};

gdap.G2longitude(gdap.G2longitude<0) = gdap.G2longitude(gdap.G2longitude<0)+360;

% npac_gdap_index = gdap.G2latitude>20 & gdap.G2latitude<60 & gdap.G2longitude>100 & gdap.G2longitude<240 & ~isnan(gdap.GMT_Matlab);
SO_gdap_index = gdap.G2latitude<80 & ~isnan(gdap.GMT_Matlab);

for g = 1:length(gdap_fields_no_flags)
    gdap_SO.(gdap_fields_no_flags{g}) = gdap.(gdap_fields_no_flags{g})(SO_gdap_index);
end
clear g

for g = 1:length(gdap_fields_w_flags)
    temp_data = gdap.(gdap_fields_w_flags{g});
    temp_data(gdap.([gdap_fields_w_flags{g} 'f'])~=2) = nan;
    
    gdap_SO.(gdap_fields_w_flags{g}) =temp_data(SO_gdap_index);
    
end
clear g SO_gdap_index temp_data
gdap_SO.G2MLD = NaN(size(gdap_SO.GMT_Matlab));

% Gdap mixed layer calcs
cruise_numbers = unique(gdap_SO.G2cruise);
stations_per_cruise = NaN(size(cruise_numbers));
for c = 1:length(cruise_numbers)
    cruise_index = find(gdap_SO.G2cruise==cruise_numbers(c));
    
    station_numbers =  unique(gdap_SO.G2station(cruise_index));
    
    stations_per_cruise(c) = numel(station_numbers);
    for s=1:length(station_numbers)
        station_index = find(gdap_SO.G2cruise==cruise_numbers(c) & gdap_SO.G2station==station_numbers(s));
        
        temp_T_0 = gdap_SO.G2temperature(station_index);
        temp_S_0 = gdap_SO.G2salinity(station_index);
        temp_P_0 = gdap_SO.G2pressure(station_index);
        
        temp_T = temp_T_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        temp_S = temp_S_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        temp_P = temp_P_0(~isnan(temp_T_0) & ~isnan(temp_S_0) & ~isnan(temp_P_0));
        
        try
        gdap_SO.G2MLD(station_index) = mld_dbm(temp_T, temp_S, ...
            temp_P, 0);
        catch
        end
        %         clear temp_P temp_P_0 temp_S temp_S_0 temp_T temp_T_0 station_index
    end
    
    %     clear cruise_index station_numbers stations_per_cruise s
end

clear c cruise_numbers stations_per_cruise cruise_index station_numbers stations_per_cruise s temp_P temp_P_0 temp_S temp_S_0 temp_T temp_T_0 station_index
gdap_SO.G2o2sat = GGo2_units(gdap_SO.G2temperature , gdap_SO.G2salinity , 'umol');
gdap_SO.G2PTMP = sw_ptmp(gdap_SO.G2salinity, gdap_SO.G2temperature, gdap_SO.G2pressure, 0);
gdap_SO.G2deltaO2 = gdap_SO.G2oxygen - gdap_SO.G2o2sat;

clear gdap
%% not running for now -  calculate LIPHR pH at Glodap points below 1480 m and above 2020m

Coordinates = [gdap_SO.G2longitude gdap_SO.G2latitude, gdap_SO.G2pressure];

Measurements = [gdap_SO.G2salinity, gdap_SO.G2temperature, gdap_SO.G2nitrate, gdap_SO.G2oxygen];
    
MeasIDVec = [1 7 3 6];

[pHEstimates,UncertaintyEstimates,MinUncertaintyEquation]= ...
    LIPHR(Coordinates,Measurements,MeasIDVec, 'OAAdjustTF', false)  ;                                  

gdap_SO.pH_in_situ_total = pHEstimates;

gdap_SO.pH_in_situ_total(isnan(gdap_SO.G2phts25p0))=nan;

%% change Glodap names to Argo names
name_convert = {'G2longitude' 'LONGITUDE'; 'G2latitude', 'LATITUDE'; 'G2pressure', 'PRES_ADJUSTED'; 'G2temperature', 'TEMP_ADJUSTED'; 'G2salinity', 'PSAL_ADJUSTED';...
    'G2oxygen' 'DOXY_ADJUSTED'; 'G2nitrate' 'NITRATE_ADJUSTED';'G2tco2' 'DIC' ; 'G2talk' 'TALK_LIAR' ; 'G2MLD' 'MLD'; 'G2o2sat' 'o2sat' ; 'G2PTMP' 'PTMP';'pH_in_situ_total' 'PH_IN_SITU_TOTAL_ADJUSTED'};

for g = 1:length(name_convert)
    gdap_SO.(name_convert{g,2}) = gdap_SO.(name_convert{g,1});
    gdap_SO.(name_convert{g,1}) = [];
end
%% - skip for now

% gdap pH 25C
[DATA,~,~]=CO2SYSSOCCOM_smb(2300.*ones(length(gdap_SO.TEMP_ADJUSTED),1), gdap_SO.PH_IN_SITU_TOTAL_ADJUSTED, ...
    1,3, gdap_SO.PSAL_ADJUSTED, gdap_SO.TEMP_ADJUSTED, 25.*ones(length(gdap_SO.TEMP_ADJUSTED),1),...
    gdap_SO.PRES_ADJUSTED, gdap_SO.PRES_ADJUSTED, zeros(length(gdap_SO.TEMP_ADJUSTED),1), zeros(length(gdap_SO.TEMP_ADJUSTED),1),1,10,3);
gdap_SO.pH_25C_TOTAL = DATA(:,37);

% set pH to nan where there was no original pH data from GLODAP
gdap_SO.pH_25C_TOTAL(isnan(gdap_SO.G2phts25p0))=nan;


%%
qc_data_fields = {'TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED',  'PRES_ADJUSTED'};


for f = 1:length(SNs)
    
    % set bad data and possibly bad data to NaN (also infs)
    for qc = 1:length(qc_data_fields)
        if isfield(Argo.(SNs{f}), qc_data_fields{qc})
            bad_index =  Argo.(SNs{f}).([qc_data_fields{qc} '_QC'])=='3' | Argo.(SNs{f}).([qc_data_fields{qc} '_QC'])=='4';
            Argo.(SNs{f}).(qc_data_fields{qc})(bad_index) = nan;
        
            inf_index = isinf(Argo.(SNs{f}).(qc_data_fields{qc}));

            if sum(reshape(inf_index,[],1))>0
                disp(f)
                Argo.(SNs{f}).(qc_data_fields{qc})(inf_index)=nan;
            end
        end
    end
end

clear bad_index
%% Calculate pH adjustments and DIC for floats
% find pH_in
last = 1;
%%
tic

for f = last:length(SNs)
    last = f;
%     % set bad data and possibly bad data to NaN
%     for qc = 1:length(qc_data_fields)
%         if isfield(Argo.(SNs{f}), qc_data_fields{qc})
%             bad_index =  Argo.(SNs{f}).([qc_data_fields{qc} '_QC'])=='3' | Argo.(SNs{f}).([qc_data_fields{qc} '_QC'])=='4';
%             Argo.(SNs{f}).(qc_data_fields{qc})(bad_index) = nan;
%         end
%     end
    
   if isfield(Argo.(SNs{f}), 'PH_IN_SITU_TOTAL_ADJUSTED') && sum(~isnan(reshape(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED,[],1)))~=0
  
        
%        % calculate potential temperature   
%        Argo.(SNs{f}).PTMP = sw_ptmp(Argo.(SNs{f}).PSAL_ADJUSTED, Argo.(SNs{f}).TEMP_ADJUSTED, Argo.(SNs{f}).PRES_ADJUSTED, 1);
% 

       % Apply bias correction
       Argo.(SNs{f}).pH_insitu_corr = NaN(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       Argo.(SNs{f}).TALK_LIAR = NaN(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       Argo.(SNs{f}).pH_25C_TOTAL = NaN(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       Argo.(SNs{f}).bias_corr = NaN(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       Argo.(SNs{f}).DIC = NaN(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       
       if isfield(Argo.(SNs{f}), 'NITRATE_ADJUSTED')
           SI = Argo.(SNs{f}).NITRATE_ADJUSTED.*2.5;
           SI(isnan(SI)) = 0;
           PO4 = Argo.(SNs{f}).NITRATE_ADJUSTED./16;
           PO4(isnan(PO4)) = 0;
       else
           SI = zeros(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
           PO4 = zeros(size(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED));
       end
       
        LAT_MATR = repmat(Argo.(SNs{f}).LATITUDE, 1, size(Argo.(SNs{f}).PRES_ADJUSTED,2));
       LON_MATR = repmat(Argo.(SNs{f}).LONGITUDE, 1, size(Argo.(SNs{f}).PRES_ADJUSTED,2));
       
       Temp_PRES = reshape(Argo.(SNs{f}).PRES_ADJUSTED,[],1);
       Temp_LAT = reshape(LAT_MATR,[],1);
       Temp_LON = reshape(LON_MATR,[],1);
       
       Coordinates = [Temp_LON, Temp_LAT, Temp_PRES];
       
       Temp_PSAL = reshape(Argo.(SNs{f}).PSAL_ADJUSTED,[],1);
       Temp_TEMP = reshape(Argo.(SNs{f}).TEMP_ADJUSTED,[],1);
       Temp_DOXY = reshape(Argo.(SNs{f}).DOXY_ADJUSTED,[],1);
       
       if isfield(Argo.(SNs{f}), 'NITRATE_ADJUSTED')
           
           Temp_NITR = reshape(Argo.(SNs{f}).NITRATE_ADJUSTED,[],1);
           Measurements = [Temp_PSAL, Temp_TEMP, Temp_NITR, Temp_DOXY];
           MeasIDVec =[1 7 3 6];
       else
           Measurements = [Temp_PSAL, Temp_TEMP, Temp_DOXY];
           MeasIDVec =[1 7 6];
       end
       [AlkalinityEstimates,~,~]= ...
           LIAR(Coordinates,Measurements,MeasIDVec, 'VerboseTF', false);
       Argo.(SNs{f}).TALK_LIAR =   reshape(AlkalinityEstimates,size(Argo.(SNs{f}).PRES_ADJUSTED,1),size(Argo.(SNs{f}).PRES_ADJUSTED,2));
       
       
       for p=1:length(Argo.(SNs{f}).GMT_Matlab)
           
           % skip a profile if pH is above 10.  There seem to be pH's above
           % 10 that causing CO2SYS to hang up and probably not even worth
           % considering otherwise
           if sum(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED(p,:)>10)>0 || sum(~isnan(Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED(p,:)))==0
               continue
           end
           disp(['Profile ' num2str( p) ' of ' num2str(length(Argo.(SNs{f}).GMT_Matlab))])
           
           % Calculate Alkalinity
           
%            Coordinates = [Argo.(SNs{f}).LONGITUDE(p).*ones(size(Argo.(SNs{f}).PRES_ADJUSTED(p,:)))', ...
%                Argo.(SNs{f}).LATITUDE(p).*ones(size(Argo.(SNs{f}).PRES_ADJUSTED(p,:)))', Argo.(SNs{f}).PRES_ADJUSTED(p,:)'];
%            
%            if isfield(Argo.(SNs{f}), 'NITRATE_ADJUSTED')
%                
%                Measurements = [Argo.(SNs{f}).PSAL_ADJUSTED(p,:)', Argo.(SNs{f}).TEMP_ADJUSTED(p,:)', Argo.(SNs{f}).NITRATE_ADJUSTED(p,:)', Argo.(SNs{f}).DOXY_ADJUSTED(p,:)'];
%                MeasIDVec =[1 7 3 6];
%            else
%                Measurements = [Argo.(SNs{f}).PSAL_ADJUSTED(p,:)', Argo.(SNs{f}).TEMP_ADJUSTED(p,:)', Argo.(SNs{f}).DOXY_ADJUSTED(p,:)'];
%                MeasIDVec =[1 7 6];
%            end
%                
%            [AlkalinityEstimates,~,~]= ...
%                LIAR(Coordinates,Measurements,MeasIDVec, 'VerboseTF', false);
%            
%            Argo.(SNs{f}).TALK_LIAR(p,:) = AlkalinityEstimates;

           
           % Calculate pH 25C without a bias correction
           
           [DATA,~,~]=CO2SYSSOCCOM(Argo.(SNs{f}).TALK_LIAR(p,:), Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED(p,:) , ...
               1,3,Argo.(SNs{f}).PSAL_ADJUSTED(p,:),Argo.(SNs{f}).TEMP_ADJUSTED(p,:),25.*ones(length(Argo.(SNs{f}).TEMP_ADJUSTED(p,:)),1),...
               Argo.(SNs{f}).PRES_ADJUSTED(p,:),Argo.(SNs{f}).PRES_ADJUSTED(p,:),SI(p,:),PO4(p,:),1,10,3);
       
           Argo.(SNs{f}).pH_25C_TOTAL(p,:) = DATA(:,37);
           Argo.(SNs{f}).DIC(p,:) = DATA(:,2);


           if sum(Argo.(SNs{f}).PRES_ADJUSTED(p,:)>1480 & Argo.(SNs{f}).PRES_ADJUSTED(p,:)<1520)>0
               correction = -0.034529.*Argo.(SNs{f}).pH_25C_TOTAL(p,Argo.(SNs{f}).PRES_ADJUSTED(p,:)>1480 & Argo.(SNs{f}).PRES_ADJUSTED(p,:)<1520) + 0.26709;
           else
               correction = -0.034529.*Argo.(SNs{f}).pH_25C_TOTAL(p,Argo.(SNs{f}).PRES_ADJUSTED(p,:)>970 & Argo.(SNs{f}).PRES_ADJUSTED(p,:)<1520) + 0.26709;
           end
           if ~isempty(correction)
               Argo.(SNs{f}).bias_corr(p) = nanmean(correction);
               Argo.(SNs{f}).pH_insitu_corr(p,:) = Argo.(SNs{f}).PH_IN_SITU_TOTAL_ADJUSTED(p,:) + Argo.(SNs{f}).bias_corr(p);
           end
           
           
       end
       disp([SNs{f} ' ' num2str(f/length(SNs)*100) ' % done'])

   end

end

disp('finished')
toc
%%  Looking for biases for each float
% Plot_dir = [home_dir 'Work/Proposals/2020_11 NOAA GOMO Float BGC Change Synthesis/plots/all_float_crossovers_float_v_float/'];
% Plot_dir = [home_dir 'Work/Proposals/2020_11 NOAA GOMO Float BGC Change Synthesis/plots/all_float_crossovers_float_v_float/'];
project_dir = [home_dir 'Work/Projects/2021_07_Float_BGC_QC_NOAA/code/'];
Plot_dir = [project_dir '../plots_matlab/Matlab_Offsets/'];

dist = 50; % km
clear offsets

meta_data = {'GMT_Matlab', 'LATITUDE', 'LONGITUDE'};
% comp_data = {'TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'pH_25C_TOTAL', 'PDENS', 'PRES_ADJUSTED', 'DIC'};
comp_data = {'TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PDENS', 'PRES_ADJUSTED'};

plot_on=0;
last = 1;
interp_press_range = 1200:1:2001;

%% while testing comparisons b/t python and matlab, only load SNs that are in the python glodap_offsets nc file:


gdap_offsets.info = ncinfo([project_dir 'output/glodap_offsets.nc']);

for v = 1:length(gdap_offsets.info.Variables)
    gdap_offsets.(gdap_offsets.info.Variables(v).Name) = ncread([project_dir 'output/glodap_offsets.nc'], ...
        gdap_offsets.info.Variables(v).Name);
end

list_python_wmos = unique(gdap_offsets.main_float_wmo);

% make an alternate list of float SNs:

orig_SNs = SNs;

SNs = cell(length(list_python_wmos),1);

for f = 1:length(list_python_wmos)
    SNs{f} = ['f' num2str(list_python_wmos(f))];
end
%%

for q= last:length(SNs) 
    %
    last=q;
    disp([num2str(q) ' ' SNs{q}])
    
    plot_filename = ['Float v Float offsets ' SNs{q} ' dens match'];
    
    comp_data_to_run = [];
    for cd = 1:length(comp_data)
        if isfield(Argo.(SNs{q}), comp_data{cd})
            if sum(~isnan(reshape(Argo.(SNs{q}).(comp_data{cd}),[],1)))>0
                comp_data_to_run = [comp_data_to_run cd];
            end
        end
    end
    
    if length(comp_data_to_run)<4
        disp(['No non-NAN bgc adjusted data for: ' SNs{q}])
        continue
    end
    
    % search other floats, looking for data within 200 km
    
    % some code adapted from Float_vs_GLODAPv2_Z_Date.m from J. Plant
    
    % SDN LAT LON
    %     t1 = track(:,3) == -1e10; % SET MISSING VALUES IN LAT to NaN
    %     track(t1,3:4) = NaN;
    %     clear t1
    %
    % GLODAP LON -180/+180 => CONVERT LON TO -180 to +180 IF NECESSARY
    %     t1 = track(:,4) > 180; % WOA2013 -180 + 180
    %     track(:,4) = track(:,4) - (t1*360);
    
    % CONVERT TOL TO DEGREES LAT AND LON
    lat_tol = dist/ 111.6; % aprox degrees lat
    lon_tol = dist/ (111.6*cosd(nanmean(Argo.(SNs{q}).LATITUDE))) ; % aprox deg lon
    
    lat_range = [Argo.(SNs{q}).LATITUDE-lat_tol Argo.(SNs{q}).LATITUDE+lat_tol];
    lon_range = [Argo.(SNs{q}).LONGITUDE-lon_tol Argo.(SNs{q}).LONGITUDE+lon_tol];
    
    %
    %     % GET LAT BOUNDS - Do 2nd subset
    %     lat_bnds = [min(track(:,3)) - lat_tol, max(track(:,3) + lat_tol)];
    %     if lat_bnds(1) < -90
    %         lat_bnds(1) = -90;
    %     elseif lat_bnds(2) > 90;
    %         lat_bnds(2) = 90;
    %     end
    %
    %  % GET LON BOUNDS - 2nd subset
    %     lon_bnds = [min(track(:,4)) - lon_tol, max(track(:,4) + lon_tol)];
    %     if lon_bnds(1) < -180
    %         disp('longitude bounds - tol crosses +180 / -180 merdian')
    %         cross180 = 1;
    %         lon_bnds(1) = lon_bnds(1)+360;
    %         lon_bnds = sort(lon_bnds); % this will reverse order
    %     elseif lon_bnds(2) > 180;
    %         disp('longitude bounds + tol crosses +180 / -180 merdian')
    %         cross180 = 1;
    %         lon_bnds(2) = lon_bnds(2)-360;
    %         lon_bnds = sort(lon_bnds); % this will reverse order
    %     end
    
    % find any other float profile that crosses this float and save the float SN in the
    % match list, along with the profiles that fall within the lat and lon test
    
    match_list.float_p = {};
    match_list.test_SN = {};
    match_list.test_p = [];
    
    %
    for f = 1:length(SNs)
        
        if f==q || ~isfield(Argo.(SNs{f}), 'PDENS')
            continue
        end
        for p = 1:length(Argo.(SNs{f}).GMT_Matlab)
            
            % is latitude in range
            lat_test = Argo.(SNs{f}).LATITUDE(p)>=lat_range(:,1) &  Argo.(SNs{f}).LATITUDE(p)<=lat_range(:,2);
            
            lon_test = Argo.(SNs{f}).LONGITUDE(p)>=lon_range(:,1) &  Argo.(SNs{f}).LONGITUDE(p)<=lon_range(:,2);
            
            if sum(lat_test & lon_test)>0
                % profiles in the original float that match test float's
                % profile location
                match_list.float_p{end+1,1} = find(lat_test==1 & lon_test==1);
                match_list.test_SN{end+1,1} = SNs{f};
                match_list.test_p(end+1,1) = p;
            end
            
        end
    end
    %
    for md = 1:length(meta_data)
        match_list.(meta_data{md}) = cell(size(match_list.float_p,1),1);
    end
    
    for cd = 1:length(comp_data)
        match_list.(comp_data{cd}) = cell(size(match_list.float_p,1),1);
    end
    %
    % if isempty(match_list)
    %     continue
    % end
    %
    if plot_on==1
        clf
        set(gcf, 'units', 'inches')
        paper_w = 12; paper_h =8;
        set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h
        
        
        d(1) = subplot(2,3,1);
        hold on; title(SNs{q})
        xlabel('Lon'); ylabel('Lat');
        d(2) = subplot(2,3,2);
        hold on; title('Temp'); grid on
        d(3) = subplot(2,3,3);
        hold on; title('O2'); grid on
        d(4) = subplot(2,3,4);
        hold on; title('pH 25C'); grid on
        d(5) = subplot(2,3,5);
        hold on; title('Nitrate'); grid on
        d(6) = subplot(2,3,6);
        hold on; title('Pot. Dens'); grid on
    end
    
    
    
    
    % loop through each profile from main float
    for m = 1:size(match_list.float_p,1)
        % add a cell to the match_list array - this is for saving crossover
        
        if plot_on==1
            % plot all test float locations and overlapping profiles
            p1 = plot(d(1), Argo.(match_list.test_SN{m}).LONGITUDE, Argo.(match_list.test_SN{m}).LATITUDE, '.k');
            p2 = plot(d(1), Argo.(match_list.test_SN{m}).LONGITUDE(match_list.test_p(m)), Argo.(match_list.test_SN{m}).LATITUDE(match_list.test_p(m)), 'b^', 'linewidth', 2);
        end
        
        
        % loop through each of the matched profiles from the main float
        for k = 1:length(match_list.float_p{m})
            % create an interpolated set of data to match to:
            temp_press = Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k), Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k),:)>100);
            
            if length(temp_press)==1 %|| length(unique(temp_press))~=length(temp_press)
                %                 if length(unique(temp_press))~=length(temp_press)
                %                     disp(['non unique pressures ' SNs{q} ' profile ' num2str(match_list.float_p{m}(k))])
                %                 end
                continue
            end
            
            %
            % interp_val has the comparison data for the main float profile
            % interpolated to a 1 db vector
            for cd = comp_data_to_run
                temp_val = Argo.(SNs{q}).(comp_data{cd})(match_list.float_p{m}(k), Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k),:)>100);
                if sum(~isnan(temp_val))>1 % need more than 1 value in order to interpolate
                    % checks for repeating values in the pressure array,
                    % assumes this is a glitch in data transmission and
                    % corrects
                    new_press = temp_press;
                    if length(unique(temp_press))~=length(temp_press)
                        temp_array = [new_press' temp_val'];
                        unique_array = unique(temp_array, 'rows');
                        new_press=unique_array(:,1)';
                        temp_val = unique_array(:,2)';
                    end
                    try
                        interp_val.(comp_data{cd}) = interp1(new_press(~isnan(new_press) & ~isnan(temp_val)), temp_val(~isnan(new_press)& ~isnan(temp_val)), interp_press_range);
                    catch
                        disp(['Catch during interpolation of ' SNs{q} ' profile ' num2str(p)])
                        interp_val.(comp_data{cd}) = NaN(1,length(interp_press_range));

                    end
                else
                    interp_val.(comp_data{cd}) = NaN(1,length(interp_press_range));
                end
            end
            interp_dens = sw_pden(interp_val.PSAL_ADJUSTED, interp_val.TEMP_ADJUSTED, interp_press_range, 0);
            
            % instead of cycling through a index from the main float which is now on an interpolated vector, cycle
            % through an index from the test float
            test_index = find(Argo.(match_list.test_SN{m}).PRES_ADJUSTED(match_list.test_p(m),:)>1480);
            %             press_index = find(Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k),:)>1480);
            
            %     loop through pressures between 1500 and 2000 db and find closest test
            %     float samples
            for y = 1:length(test_index)
                
                %                 % finding float matches based on pressure
                %                 matched_press_index = find(min(abs(Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k),press_index(y)) - ...
                %                     Argo.(match_list.test_SN{m}).PRES_ADJUSTED(match_list.test_p(m),:)))== ...
                %                     abs(Argo.(SNs{q}).PRES_ADJUSTED(match_list.float_p{m}(k),press_index(y)) - ...
                %                     Argo.(match_list.test_SN{m}).PRES_ADJUSTED(match_list.test_p(m),:)));
                
                % finding the closest density from the interpolated fields
                interp_index = find(min(abs(interp_dens - Argo.(match_list.test_SN{m}).PDENS(match_list.test_p(m),test_index(y))))==abs(interp_dens - Argo.(match_list.test_SN{m}).PDENS(match_list.test_p(m),test_index(y))));
                
                for md = 1:length(meta_data)
                    temp_val = Argo.(SNs{q}).(meta_data{md})(match_list.float_p{m}(k));
                    var_vector = NaN(1,length(interp_index));
                    var_vector(:) = temp_val;
                    match_list.(meta_data{md}){m} = [ match_list.(meta_data{md}){m} var_vector];
                    
                end
                
                
                for cd = comp_data_to_run
                    if isfield(Argo.(match_list.test_SN{m}), comp_data{cd})
                        %
                        %                         match_list.(comp_data{cd}){m} = [match_list.(comp_data{cd}){m} Argo.(SNs{q}).(comp_data{cd})(match_list.float_p{m}(k), press_index(y)) -  ...
                        %                             Argo.(match_list.test_SN{m}).(comp_data{cd})(match_list.test_p(m), matched_press_index)];
                        %
                        match_list.(comp_data{cd}){m} = [match_list.(comp_data{cd}){m} interp_val.(comp_data{cd})(interp_index) -  ...
                            Argo.(match_list.test_SN{m}).(comp_data{cd})(match_list.test_p(m), test_index(y))];
                    else
                        match_list.(comp_data{cd}){m} = [match_list.(comp_data{cd}){m} NaN(1,length(interp_index))];
                    end
                end
                
                
            end
        end
        
    end
    
    clear temp_val var_vector
    
    if plot_on==1
        p3 = plot(d(1), Argo.(SNs{q}).LONGITUDE, Argo.(SNs{q}).LATITUDE, 'm.', 'linewidth', 2);
        
        legend(d(1), [p1 p2 p3], 'Comparison floats', 'Matched profiles', 'Current float', 'location', 'southoutside')
        % test_vars = {'TEMP_ADJUSTED' 4; 'o2_umol_kg' 5; 'pH_25C' 6; 'Nitrate' 7; 'PDENS' 8};
    end
    
    for t = 1:length(meta_data)
        offsets.(SNs{q}).(meta_data{t}) = [];
    end
    for t = 1:length(comp_data)
        offsets.(SNs{q}).(comp_data{t}) = [];
    end
    
    for m = 1:size(match_list.float_p,1)
        
        for t = 1:length(meta_data)
            offsets.(SNs{q}).(meta_data{t}) = [offsets.(SNs{q}).(meta_data{t}) match_list.(meta_data{t}){m}];
            if t==1
                var_length = length(offsets.(SNs{q}).(meta_data{t}));
            else
                if length(offsets.(SNs{q}).(meta_data{t}))~= var_length
                    %                 disp(m)
                    %                 break
                end
            end
        end
        
    end
    
    for m = 1:size(match_list.float_p,1)
        
        for t = 1:length(comp_data)
            offsets.(SNs{q}).(comp_data{t}) = [offsets.(SNs{q}).(comp_data{t}) match_list.(comp_data{t}){m}];
            if t==1
                var_length = length(offsets.(SNs{q}).(comp_data{t}));
            else
                if length(offsets.(SNs{q}).(comp_data{t}))~= var_length
                    %                 disp(m)
                    %                 break
                end
            end
        end
        
    end
    
    
    if plot_on==1
        for t = 1:length(comp_data)
            if ~isempty(offsets.(SNs{q}).(comp_data{t}))
                histogram(d(t+1), offsets.(SNs{q}).(comp_data{t}))
                
                histogram(d(t+1), offsets.(SNs{q}).(comp_data{t})(abs(offsets.(SNs{q}).PDENS)<0.03))
                
                xlabel(d(t+1), ['Mean: ' num2str(nanmean(offsets.(SNs{q}).(comp_data{t})(abs(offsets.(SNs{q}).PDENS)<0.03)),3)])
            end
            plot(d(t+1), [0 0], get(d(t+1), 'ylim'), 'k')
            
        end
        % pause
        clear p1 p2 p3
        print(gcf, '-dpng', '-r400', [Plot_dir  plot_filename '.png' ])
        pause
    end
end

%%  Glodap crossovers
c_map = brewermap(6, 'Set1');

gdap_fields = fieldnames(gdap_SO);
plot_on=0;
plot_final=1;

var_to_plot = 'DOXY_ADJUSTED';
c_p = 3;

[X_S, Y_T] = meshgrid(33:.1:36, -2:.5:16);
dens_grid = sw_pden(X_S, Y_T, 1500, 1)-1000;


for q= 1:length(SNs)
    %     last=q;
    
    disp([num2str(q) ' ' SNs{q}])
    
    if ~isfield(Argo.(SNs{q}), 'PDENS')
        continue
    end
    
    
    comp_data_to_run = [];
    for cd = 1:length(comp_data)
        if isfield(Argo.(SNs{q}), comp_data{cd})
            if sum(~isnan(reshape(Argo.(SNs{q}).(comp_data{cd}),[],1)))>0
                comp_data_to_run = [comp_data_to_run cd];
            end
        end
    end
    
    
    if length(comp_data_to_run)<5
        disp(['No non-NAN bgc adjusted data for: ' SNs{q}])
        continue
    end
    
    
    
    % CONVERT TOL TO DEGREES LAT AND LON
    lat_tol = dist/ 111.6; % aprox degrees lat
    lon_tol = dist/ (111.6*cosd(nanmean(Argo.(SNs{q}).LATITUDE))) ; % aprox deg lon
    
    lat_range = [min(Argo.(SNs{q}).LATITUDE-lat_tol) max(Argo.(SNs{q}).LATITUDE+lat_tol)];
    lon_range = [min(Argo.(SNs{q}).LONGITUDE-lon_tol) max(Argo.(SNs{q}).LONGITUDE+lon_tol)];
    
    
    lat_test = gdap_SO.LATITUDE>=lat_range(:,1) &  gdap_SO.LATITUDE<=lat_range(:,2);
    
    % check if the float crossed the prime meridian
    if lon_range(1)<0
        lon_test = gdap_SO.LONGITUDE>=360+lon_range(:,1) | gdap_SO.LONGITUDE<=lon_range(:,2);
    elseif lon_range(2)>360
        lon_test = gdap_SO.LONGITUDE>=lon_range(:,1) | gdap_SO.LONGITUDE<=lon_range(:,2)-360;
    else
        lon_test = gdap_SO.LONGITUDE>=lon_range(:,1) &  gdap_SO.LONGITUDE<=lon_range(:,2);
    end
    
    % check if any glodap data is close to this float
    if sum(lat_test & lon_test)==0
        continue
        
    end
    % create a subset of the glodap database with just the data that are
    % close
    clear temp_gdap
    for g = 1:length(gdap_fields)
        if ~isempty(gdap_SO.(gdap_fields{g}))
            temp_gdap.(gdap_fields{g}) = gdap_SO.(gdap_fields{g})(lat_test & lon_test);
        end
    end
    %
    clear lat_test lon_test
    % go through each float profile and find any glodap data that is close
    
    % prep offsets to enter glodap data
    for t = 1:length(meta_data)
        offsets.(SNs{q}).gdap.([meta_data{t} '_float']) = [];
        offsets.(SNs{q}).gdap.([meta_data{t} '_gdap']) = [];

    end
    for t = 1:length(comp_data)
        offsets.(SNs{q}).gdap.([comp_data{t} '_offset']) = [];
        offsets.(SNs{q}).gdap.([comp_data{t} '_float']) = [];
        offsets.(SNs{q}).gdap.([comp_data{t} '_gdap']) = [];

    end
    
    
    
    for p = 1:length(Argo.(SNs{q}).GMT_Matlab)
        clear interp_val
        lat_range = [Argo.(SNs{q}).LATITUDE(p)-lat_tol Argo.(SNs{q}).LATITUDE(p)+lat_tol];
        lon_range = [Argo.(SNs{q}).LONGITUDE(p)-lon_tol Argo.(SNs{q}).LONGITUDE(p)+lon_tol];
        
        lat_test = temp_gdap.LATITUDE>=lat_range(:,1) &  temp_gdap.LATITUDE<=lat_range(:,2);
        if lon_range(1)<0
            lon_test = temp_gdap.LONGITUDE>=360+lon_range(:,1) | temp_gdap.LONGITUDE<=lon_range(:,2);
        elseif lon_range(2)>360
            lon_test = temp_gdap.LONGITUDE>=lon_range(:,1) | temp_gdap.LONGITUDE<=lon_range(:,2)-360;
        else
            lon_test = temp_gdap.LONGITUDE>=lon_range(:,1) &  temp_gdap.LONGITUDE<=lon_range(:,2);
        end
        
        if sum(lon_test & lat_test)==0
            continue
        end
        if plot_on==1
            clf
        end
        % get pressures of the float profile between 1500 and 2000 m
        
        %         press_index = find(Argo.(SNs{q}).PRES_ADJUSTED(p,:)>1480);
        
        temp_press = Argo.(SNs{q}).PRES_ADJUSTED(p,Argo.(SNs{q}).PRES_ADJUSTED(p,:)>1050);
        if length(temp_press)==1
            continue
        end
        %         temp_T = Argo.(SNs{q}).TEMP_ADJUSTED(p,:);
        %         temp_S = Argo.(SNs{q}).PSAL_ADJUSTED(p,:);
        %
        %         interp_temp = interp1(temp_press(~isnan(temp_press)), temp_T(~isnan(temp_press)), interp_press_range);
        %         interp_sal = interp1(temp_press(~isnan(temp_press)), temp_S(~isnan(temp_press)), interp_press_range);
        %
        
        % interp_val has the comparison data for the main float profile
        % interpolated to a 1 db vector
        for cd = comp_data_to_run
            temp_val = Argo.(SNs{q}).(comp_data{cd})(p,Argo.(SNs{q}).PRES_ADJUSTED(p,:)>1050);
            temp_press2 = temp_press;
            if sum(~isnan(temp_val))>1 % need at least two values to run interp1
                % checks for repeating values in the pressure array,
                % assumes this is a glitch in data transmission and
                % corrects
                if length(unique(temp_press2))~=length(temp_press2)
                    
                    temp_array = [temp_press2' temp_val'];
                    unique_array = unique(temp_array, 'rows');
                    temp_press2=unique_array(:,1)';
                    temp_val = unique_array(:,2)';
                end
                try
                    interp_val.(comp_data{cd}) = interp1(temp_press2(~isnan(temp_press2) & ~isnan(temp_val)), temp_val(~isnan(temp_press2)& ~isnan(temp_val)), interp_press_range);
                catch
                    disp(['Catch during interpolation of ' SNs{q} ' profile ' num2str(p)])
                end
            else
                interp_val.(comp_data{cd}) = NaN(1,length(interp_press_range));
            end
            

        end
        
        % check whether there are any biogeochemical comp_data with numbers:
        data_true = 0;
        for cd = comp_data_to_run
            if strncmp(comp_data{cd}, 'DOXY', 4) || ...
                    strncmp(comp_data{cd}, 'NITRATE', 4) || strncmp(comp_data{cd}, 'pH_25', 4)   || strncmp(comp_data{cd}, 'DIC', 3)
                if sum(~isnan(interp_val.(comp_data{cd})))>0
                    data_true= 1;
                end
            end
        end
        if data_true==0
            continue
        end
        interp_dens = sw_pden(interp_val.PSAL_ADJUSTED, interp_val.TEMP_ADJUSTED, interp_press_range, 0);
        if plot_on==1

            % plot var vs. pressure
            d1 = subplot(3,2,1);
            plot(d1, Argo.(SNs{q}).(comp_data{c_p})(p,:), Argo.(SNs{q}).PRES_ADJUSTED(p,:), 'o')
            hold on
            plot(d1, interp_val.(comp_data{c_p}), interp_press_range, '-x')
%             disp(p)
            
            d2 = subplot(3,2,2);
            plot(d2, Argo.(SNs{q}).PDENS(p,:), Argo.(SNs{q}).PRES_ADJUSTED(p,:), 'o')
            hold on
            plot(d2, interp_dens, interp_press_range, '-x')
            %            disp(p)
            
            d3 = subplot(3,2,3);
            plot(d3, Argo.(SNs{q}).(comp_data{c_p})(p,:), Argo.(SNs{q}).PDENS(p,:), 'o')
            hold on
            plot(d3, interp_val.(comp_data{c_p}), interp_dens, '-x')

           
        end
        if sum(~isnan(interp_dens))==0
            continue
        end
        
        % find all glodap points between 1480 and 2010 and within the
        % latitude/longitude boundaries
        gdap_press_index = find(temp_gdap.PRES_ADJUSTED>1480 & temp_gdap.PRES_ADJUSTED<2010 & lon_test & lat_test);
        gdap_cruise_list = [];
        if isempty(gdap_press_index)
            continue
        end
        for y = 1:length(gdap_press_index)
            %             gdap_index = temp_gdap.PRES_ADJUSTED>Argo.(SNs{q}).PRES_ADJUSTED(p,press_index(y))-10 & temp_gdap.PRES_ADJUSTED<Argo.(SNs{q}).PRES_ADJUSTED(p,press_index(y))+10 & lon_test & lat_test;
            if isnan(temp_gdap.PDENS(gdap_press_index(y)))
                continue
            end
            
              
            % find the closest density between interpolated float data and
            % glodap potential density
            interp_index = find(min(abs(interp_dens - temp_gdap.PDENS(gdap_press_index(y))))==abs(interp_dens - temp_gdap.PDENS(gdap_press_index(y))));
            
            % don't use matches if they are at one end or the other of the
            % interpolated density vector, only if they are contained
            % within - % note that only using top or bottom doesn't work,
            % because density doesn't necessarily line up. needs to be
            % max/min or look for greater than a density instep:
          
            if interp_index==1 || interp_index==find(~isnan(interp_dens),1,'last') || abs(temp_gdap.PDENS(gdap_press_index(y)) - interp_dens(interp_index))>0.0005
                continue
            end
            
            
            %plot glodap profiles used for comparison
            if plot_on==1
                
                % find points in that profile
                gdap_profile_index = temp_gdap.G2cruise==temp_gdap.G2cruise(gdap_press_index(y)) & temp_gdap.G2station==temp_gdap.G2station(gdap_press_index(y));
                plot(d1, temp_gdap.DOXY_ADJUSTED(gdap_profile_index), temp_gdap.PRES_ADJUSTED(gdap_profile_index), '+-m')
                plot(d2, temp_gdap.PDENS(gdap_profile_index), temp_gdap.PRES_ADJUSTED(gdap_profile_index), '+-m')
                
                plot(d3, temp_gdap.DOXY_ADJUSTED(gdap_profile_index), temp_gdap.PDENS(gdap_profile_index), '+-m')
                
            end
            
            for md = 1:length(meta_data)
                temp_val =  Argo.(SNs{q}).(meta_data{md})(p);
                val_vector = NaN(length(interp_index),1);
                val_vector(:) = temp_val;
                offsets.(SNs{q}).gdap.([meta_data{md} '_float'])  = ...
                    [offsets.(SNs{q}).gdap.([meta_data{md} '_float'])  ; val_vector];
                
                temp_val =  temp_gdap.(meta_data{md})(gdap_press_index(y));
                val_vector = NaN(length(interp_index),1);
                val_vector(:) = temp_val;
                offsets.(SNs{q}).gdap.([meta_data{md} '_gdap'])  = ...
                    [offsets.(SNs{q}).gdap.([meta_data{md} '_gdap'])  ; val_vector];
                
            end
            
            for cd = comp_data_to_run
                
                offsets.(SNs{q}).gdap.([comp_data{cd} '_offset']) = [offsets.(SNs{q}).gdap.([comp_data{cd} '_offset']) ; ...
                    interp_val.(comp_data{cd})(interp_index) - temp_gdap.(comp_data{cd})(gdap_press_index(y))];
                
                offsets.(SNs{q}).gdap.([comp_data{cd} '_float']) = [offsets.(SNs{q}).gdap.([comp_data{cd} '_float']) ; ...
                    interp_val.(comp_data{cd})(interp_index)];
                
                offsets.(SNs{q}).gdap.([comp_data{cd} '_gdap']) = [offsets.(SNs{q}).gdap.([comp_data{cd} '_gdap']); ...
                    temp_gdap.(comp_data{cd})(gdap_press_index(y))];
           
                
            end
            if plot_on==1
                plot(d1, interp_val.(comp_data{c_p})(interp_index), interp_val.PRES_ADJUSTED(interp_index), 'sk', 'markersize', 15)
                
                %plot glodap match
                plot(d1, temp_gdap.(comp_data{c_p})(gdap_press_index(y)), temp_gdap.PRES_ADJUSTED(gdap_press_index(y)), 'db', 'markersize', 15)
                
                
                % density plot
                plot(d3, interp_val.(comp_data{c_p})(interp_index), interp_val.PDENS(interp_index), 'sk', 'markersize', 15)
                plot(d3, temp_gdap.(comp_data{c_p})(gdap_press_index(y)), temp_gdap.PDENS(gdap_press_index(y)), 'db', 'markersize', 15)
            end
            
            if isempty(gdap_cruise_list) || (sum(gdap_cruise_list(:,1)==temp_gdap.G2cruise(gdap_press_index(y)) & gdap_cruise_list(:,2)==temp_gdap.G2station(gdap_press_index(y)))==0)
                gdap_cruise_list = [gdap_cruise_list ; temp_gdap.G2cruise(gdap_press_index(y)) temp_gdap.G2station(gdap_press_index(y)) temp_gdap.GMT_Matlab(gdap_press_index(y))];
            end
%             pause
        end
        
        if plot_on==1
            set([d1 d2], 'ydir', 'reverse', 'ylim', [1200 2000])
            set(d3, 'ylim', [min(interp_dens)-.02 max(interp_dens)+.02])
            ylabel(d1, 'Pressure'); xlabel(d1, comp_data{c_p}, 'interpreter', 'none');
            ylabel(d2, 'Pressure'); xlabel(d2, 'PDENS');
            ylabel(d3, 'PDENS'); xlabel(d3, comp_data{c_p}, 'interpreter', 'none');
            title(d1, [SNs{q} ' Prof: ' num2str(p)])
            
            d4 = subplot(3,2,4);
            hold on
            text(0.1, .9, 'G2cruise G2station Date')
            
            for gg = 1:size(gdap_cruise_list,1)
                text(0.1, 0.9-gg/10, [num2str(gdap_cruise_list(gg,1)) ' ' num2str(gdap_cruise_list(gg,2)) ' ' datestr(gdap_cruise_list(gg,3), 'YYYY-mm-dd')])
            end
            
            % histograms
            prof_index = offsets.(SNs{q}).gdap.GMT_Matlab_float==Argo.(SNs{q}).GMT_Matlab(p);
            d5 = subplot(3,2,5);
            hold on
            histogram(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])(prof_index));
            title([comp_data{c_p} '_offset'], 'interpreter', 'none')
            d6 = subplot(3,2,6);
            histogram(offsets.(SNs{q}).gdap.PDENS_offset(prof_index));
            title('PDENS offset')
            
            if ~isfolder([Plot_dir SNs{q}(2:end)])
                mkdir([Plot_dir SNs{q}(2:end)])
            end
            
            
            plot_filename = [SNs{q}(2:end) '_Prof_' num2str(p) '_' comp_data{c_p}];
            print(gcf, '-dpng', '-r200', [Plot_dir '/' SNs{q}(2:end) '/' plot_filename '.png' ])
            
        end
        
    end
    
    if isempty(offsets.(SNs{q}).gdap.GMT_Matlab_float)
        continue
    end
   % make a plot for entire float's offsets:
   set(gcf, 'colormap', turbo)
   if plot_final==1
       clf
       set(gcf, 'units', 'inches')
       paper_w = 12; paper_h =12;
       set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);
       
       % map
       subplot(5,3,1); hold on
       plot(offsets.(SNs{q}).gdap.LONGITUDE_float, offsets.(SNs{q}).gdap.LATITUDE_float, '+', 'color', c_map(3,:))
       plot(offsets.(SNs{q}).gdap.LONGITUDE_gdap, offsets.(SNs{q}).gdap.LATITUDE_gdap, 'o', 'color', c_map(5,:))
       
       title(SNs{q}(2:end))
       subplot(5,3,2); hold on
       histogram(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))
       title([comp_data{c_p}(1:4) ' mean ' num2str(nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])),2) ...
           ' \pm ' num2str(nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])),2)])
       
       subplot(5,3,3); hold on
       scatter(offsets.(SNs{q}).gdap.LONGITUDE_float, offsets.(SNs{q}).gdap.LATITUDE_float, 30, ...
           offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), 'o', 'filled')
       caxis([nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))-2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])) ...
           nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))+ 2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))])
       colorbar
       title('Map of offsets')
       
       subplot(5,3,4:6); hold on
       plot([min(offsets.(SNs{q}).gdap.GMT_Matlab_float) max(offsets.(SNs{q}).gdap.GMT_Matlab_float)], [0 0], '--k')
       plot(offsets.(SNs{q}).gdap.GMT_Matlab_float, offsets.(SNs{q}).gdap.DOXY_ADJUSTED_offset, 'x');
       
       date_vec = datevec(offsets.(SNs{q}).gdap.GMT_Matlab_float);
       set(gca, 'xtick', datenum(min(date_vec(:,1)):1:max(date_vec(:,1)),1,1))
       datetick('x', 'YY','keepticks', 'keeplimits')
       title('Offsets vs. float date')
       
       subplot(5,3,7:9)
       plot([min(offsets.(SNs{q}).gdap.GMT_Matlab_gdap) max(offsets.(SNs{q}).gdap.GMT_Matlab_gdap)], [0 0], '--k')
       hold on
       plot(offsets.(SNs{q}).gdap.GMT_Matlab_gdap, offsets.(SNs{q}).gdap.DOXY_ADJUSTED_offset, 'x');
        date_vec = datevec(offsets.(SNs{q}).gdap.GMT_Matlab_gdap);
       set(gca, 'xtick', datenum(min(date_vec(:,1)):1:max(date_vec(:,1)),1,1))
       datetick('x', 'YY','keepticks', 'keeplimits')
      title('Offsets vs. gdap date')
       
       subplot(5,3,10)
       plot(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), offsets.(SNs{q}).gdap.PRES_ADJUSTED_float, 'x');
       ylabel('Pres')
       title('Float Pres vs offset')
       set(gca, 'ydir', 'reverse')
       
       subplot(5,3,11)
       plot(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), offsets.(SNs{q}).gdap.PRES_ADJUSTED_gdap, 'x');
       title('GDAP Pres vs offset')
       set(gca, 'ydir', 'reverse')

       subplot(5,3,13)
       plot(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), offsets.(SNs{q}).gdap.PDENS_float, 'x');
       ylabel('PDENS')
       title('Float PDENS vs offset')
       
       subplot(5,3,14)
       plot(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), offsets.(SNs{q}).gdap.PDENS_gdap, 'x');
       title('GDAP PDENS vs offset')
       

       T_lims = [min([offsets.(SNs{q}).gdap.TEMP_ADJUSTED_float; offsets.(SNs{q}).gdap.TEMP_ADJUSTED_gdap]) max([offsets.(SNs{q}).gdap.TEMP_ADJUSTED_float; offsets.(SNs{q}).gdap.TEMP_ADJUSTED_gdap])];
       S_lims = [min([offsets.(SNs{q}).gdap.PSAL_ADJUSTED_float; offsets.(SNs{q}).gdap.PSAL_ADJUSTED_gdap]) max([offsets.(SNs{q}).gdap.PSAL_ADJUSTED_float; offsets.(SNs{q}).gdap.PSAL_ADJUSTED_gdap])];

       subplot(5,3,12)
       hold on
                     contour(X_S, Y_T, dens_grid, 'ShowText', 'off', 'levellist', 24.1:.01:29, 'color', [.6 .6 .6]);

%        contour(X_S, Y_T, dens_grid, 'k', 'ShowText', 'on');
scatter(offsets.(SNs{q}).gdap.PSAL_ADJUSTED_float, offsets.(SNs{q}).gdap.TEMP_ADJUSTED_float, 20, offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), 'filled')
colorbar
set(gca, 'xlim', S_lims, 'ylim', T_lims)
title('Offsets plotted w/ float T and S')
caxis([nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))-2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])) ...
    nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))+ 2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))])
subplot(5,3,15)
hold on
contour(X_S, Y_T, dens_grid, 'ShowText', 'off', 'levellist', 24.1:.01:29, 'color', [.6 .6 .6]);
%        contour(X_S, Y_T, dens_grid, 'k', 'ShowText', 'on');
scatter(offsets.(SNs{q}).gdap.PSAL_ADJUSTED_gdap, offsets.(SNs{q}).gdap.TEMP_ADJUSTED_gdap, 20, offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']), 'filled')
colorbar
set(gca, 'xlim', S_lims, 'ylim', T_lims)

title('Offsets plotted w/ gdap T and S')
caxis([nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))-2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset'])) ...
    nanmean(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))+ 2*nanstd(offsets.(SNs{q}).gdap.([comp_data{c_p} '_offset']))])

       plot_filename = [SNs{q}(2:end) '_' comp_data{c_p}];
       print(gcf, '-dpng', '-r200', [Plot_dir '/' plot_filename '_v2.png' ])
   end
   
   disp([SNs{q} ' ' num2str(q/length(SNs)*100) ' % done'])

end
%%
clear val_vector gdap_index temp_val cd md y p lat_test lon_test t g comp_data_to_run interp_dens interp_index interp_sal interp_temp interp_val

%% 
save([home_dir 'Work/Projects/2021_07_Float_BGC_QC_NOAA/code/output/temp_Matlab_glodap_offsets.mat'], 'offsets', 'comp_data', 'meta_data');
%% counting crossovers
gdap_count.O2 = 0;
gdap_count.NO3 = 0;
gdap_count.pH = 0;

for f = 1:length(SNs)
    if ~isfield(offsets, SNs{f})
        continue
    end
    if ~isfield(offsets.(SNs{f}), 'gdap')
        continue
    end
    gdap_count.O2 = gdap_count.O2 + length(offsets.(SNs{f}).gdap.DOXY_ADJUSTED);
    gdap_count.NO3 = gdap_count.NO3 + length(offsets.(SNs{f}).gdap.NITRATE_ADJUSTED);
    gdap_count.pH = gdap_count.pH + length(offsets.(SNs{f}).gdap.pH_25C_TOTAL);

end
%% Histograms with both gdap and float data
hist_colors = brewermap(2, 'Set1');
last=1;
%%

for q = last: length(SNs)
    last=q;
    if ~isfield(offsets, SNs{q})
        continue
    end
    no_gdap = 0; no_argo = 0;
    
    
    if ~isfield(offsets.(SNs{q}), 'gdap')
        no_gdap = 1;
    elseif isempty(offsets.(SNs{q}).gdap.GMT_Matlab)
        no_gdap = 1;
    end
    if ~isfield(offsets.(SNs{q}), 'GMT_Matlab')
        no_argo=1;
    elseif isempty(offsets.(SNs{q}).GMT_Matlab)
        no_argo=1;
    end
    
    if no_argo==1 && no_gdap==1
        disp([num2str(q) ' ' SNs{q} ' no offset data'])
        continue
    end
    
    %     if ~isfield(offsets.(SNs{q}), 'gdap') &&  ~isfield(offsets.(SNs{q}), 'GMT_Matlab')
    %         disp([num2str(q) ' ' SNs{q} ' no offset data'])
    %         continue
    %     elseif isempty(offsets.(SNs{q}).GMT_Matlab)
    %         if isempty(offsets.(SNs{q}).gdap.GMT_Matlab)
    %             disp([num2str(q) ' ' SNs{q} ' no offset data'])
    %             continue
    %         end
    %     end
    
    
    plot_filename = ['Float v Float and GDAP offsets ' SNs{q} ' dens match'];
    
    clf
    set(gcf, 'units', 'inches')
    paper_w = 12; paper_h =8;
    set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h
    
    for cd = 1:length(comp_data)
        d(cd) = subplot(3,3,cd);
        hold on
    end
    %
    
    for t = 1:length(comp_data)
        argo_present=0;
        gdap_present=0;
        %         if ~isempty(offsets.(SNs{q}).(comp_data{t}))
        
        if isfield(offsets.(SNs{q}), comp_data{t}) && ~isempty(offsets.(SNs{q}).(comp_data{t}))
            if sum(~isnan(offsets.(SNs{q}).(comp_data{t})))>0
                if strcmp(comp_data{cd}, 'PDENS')
                    h1 = histogram(d(t), offsets.(SNs{q}).(comp_data{t})(:));
                    
                else
                    h1 = histogram(d(t), offsets.(SNs{q}).(comp_data{t})(abs(offsets.(SNs{q}).PDENS)<0.03));
                end
                set(h1, 'FaceColor', hist_colors(1,:));
                
            end
            argo_present=1;
        end
        
        if isfield(offsets.(SNs{q}), 'gdap')
            if isfield(offsets.(SNs{q}).gdap, comp_data{t}) && sum(~isnan(offsets.(SNs{q}).gdap.(comp_data{t})))>0
                if strcmp(comp_data{cd}, 'PDENS')
                    h2 = histogram(d(t), offsets.(SNs{q}).gdap.(comp_data{t})(:));
                else
                    h2 = histogram(d(t), offsets.(SNs{q}).gdap.(comp_data{t})(abs(offsets.(SNs{q}).gdap.PDENS)<0.03));
                end
                set(h2, 'Facecolor', hist_colors(2,:));
                
                gdap_present=1;
                
            end
        end
        if argo_present==1 && gdap_present==1
            xlabel(d(t), ['A Mean: ' num2str(nanmean(offsets.(SNs{q}).(comp_data{t})(abs(offsets.(SNs{q}).PDENS)<0.03)),3) ...
                '; G Mean: ' num2str(nanmean(offsets.(SNs{q}).gdap.(comp_data{t})(abs(offsets.(SNs{q}).gdap.PDENS)<0.03)),3) ])
        elseif argo_present==1 && gdap_present==0
            xlabel(d(t), ['A Mean: ' num2str(nanmean(offsets.(SNs{q}).(comp_data{t})(abs(offsets.(SNs{q}).PDENS)<0.03)),3) ])
        elseif argo_present==0 && gdap_present==1
            xlabel(d(t), ['G Mean: ' num2str(nanmean(offsets.(SNs{q}).gdap.(comp_data{t})(abs(offsets.(SNs{q}).gdap.PDENS)<0.03)),3) ])
        end
        
        plot(d(t), [0 0], get(d(t), 'ylim'), 'k')
        
        title(d(t), comp_data{t}, 'interpreter', 'none')
        
        
        
        if t==1
            try
                legend([h1 h2], 'Float', 'Gdap')
            catch
            end
        end
    end
    
    if no_gdap==0 && no_argo==0
        min_lat_lim = min([offsets.(SNs{q}).LATITUDE'; offsets.(SNs{q}).gdap.LATITUDE]);
        max_lat_lim = max([offsets.(SNs{q}).LATITUDE'; offsets.(SNs{q}).gdap.LATITUDE]);
        min_lon_lim = min([offsets.(SNs{q}).LONGITUDE'; offsets.(SNs{q}).gdap.LONGITUDE]);
        max_lon_lim = max([offsets.(SNs{q}).LONGITUDE'; offsets.(SNs{q}).gdap.LONGITUDE]);
        
    elseif no_argo==1
        min_lat_lim = min([offsets.(SNs{q}).gdap.LATITUDE; offsets.(SNs{q}).gdap.LATITUDE]);
        max_lat_lim = max([offsets.(SNs{q}).gdap.LATITUDE; offsets.(SNs{q}).gdap.LATITUDE]);
        min_lon_lim = min([offsets.(SNs{q}).gdap.LONGITUDE; offsets.(SNs{q}).gdap.LONGITUDE]);
        max_lon_lim = max([offsets.(SNs{q}).gdap.LONGITUDE; offsets.(SNs{q}).gdap.LONGITUDE]);
        
        
    else
        min_lat_lim = min([offsets.(SNs{q}).LATITUDE'; offsets.(SNs{q}).LATITUDE']);
        max_lat_lim = max([offsets.(SNs{q}).LATITUDE'; offsets.(SNs{q}).LATITUDE']);
        min_lon_lim = min([offsets.(SNs{q}).LONGITUDE'; offsets.(SNs{q}).LONGITUDE']);
        max_lon_lim = max([offsets.(SNs{q}).LONGITUDE'; offsets.(SNs{q}).LONGITUDE']);
        
        
        
    end
    subplot(3,3,8:9)
    
    if (max_lon_lim - min_lon_lim) > 240 % probably wraps around the prime meridian
        temp_min_lon_lim = min_lon_lim;
        min_lon_lim = max_lon_lim;
        max_lon_lim = temp_min_lon_lim + 360;
    end
    
    hold on
    coast = load('coast');
    axesm('lambertstd','MapLatLimit',[min_lat_lim-2 max_lat_lim+2], 'MapLonLimit', [min_lon_lim-2 max_lon_lim+2], 'MapParallels', [min_lat_lim-2 max_lat_lim+2])
    axis off; framem on; gridm on; mlabel on; plabel on;
    setm(gca, 'fontsize', 16)
    setm(gca,'MLabelParallel',min_lat_lim-2, 'fontsize', 16)
    
    coast_gray = [.4 .4 .4];
    geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)
    marker_size = 10;
    
    % scatterm(offsets.(SNs{q}).LATITUDE, offsets.(SNs{q}).LONGITUDE, 5, offsets.(SNs{q})
    if no_argo==0
        g1 = geoshow(offsets.(SNs{q}).LATITUDE, offsets.(SNs{q}).LONGITUDE, 'displaytype', 'multipoint', 'marker', '+', 'markeredgecolor', 'k','markerfacecolor', 'k',  'markersize', marker_size);
    end
    if no_gdap==0
        g2 = geoshow(offsets.(SNs{q}).gdap.LATITUDE, offsets.(SNs{q}).gdap.LONGITUDE, 'displaytype', 'multipoint', 'marker', '^', 'markeredgecolor', 'k','markerfacecolor', 'none',  'markersize', marker_size);
    end
    
    try
        legend([g1 g2], 'Float', 'GLODAP', 'location', 'eastoutside')
    catch
    end
    
    title(SNs{q})
    
    print(gcf, '-dpng', '-r400', [Plot_dir 'GLOBAL_dens_match/' plot_filename '.png' ])
    %     pause
end


%% assembling data for map plots
clear plotall
off_SNs = fieldnames(offsets);

argo_on = 1;
gdap_on = 1;

for cd = 1:length(comp_data)
    disp(comp_data{cd})
    
    plotall.argo.(comp_data{cd}).LATITUDE = [];
    plotall.argo.(comp_data{cd}).LONGITUDE = [];
    if argo_on==1
        
        plotall.offsets.argo.(comp_data{cd}).LATITUDE = [];
        plotall.offsets.argo.(comp_data{cd}).LONGITUDE = [];
        plotall.offsets.argo.(comp_data{cd}).val = [];
    end
    
    if gdap_on==1
        plotall.offsets.gdap.(comp_data{cd}).LATITUDE = [];
        plotall.offsets.gdap.(comp_data{cd}).LONGITUDE = [];
        plotall.offsets.gdap.(comp_data{cd}).val = [];
    end
    
    for f = 1:length(off_SNs)
        plotall.argo.(comp_data{cd}).LATITUDE = [plotall.argo.(comp_data{cd}).LATITUDE ; Argo.(off_SNs{f}).LATITUDE];
        plotall.argo.(comp_data{cd}).LONGITUDE = [plotall.argo.(comp_data{cd}).LONGITUDE ; Argo.(off_SNs{f}).LONGITUDE];
        
        if argo_on==1
            if ~isempty(offsets.(off_SNs{f}).(comp_data{cd}))
                
                
                temp_index = ~isnan(offsets.(off_SNs{f}).(comp_data{cd}));
                plotall.offsets.argo.(comp_data{cd}).LATITUDE = [plotall.offsets.argo.(comp_data{cd}).LATITUDE ; offsets.(off_SNs{f}).LATITUDE(temp_index)'];
                plotall.offsets.argo.(comp_data{cd}).LONGITUDE = [plotall.offsets.argo.(comp_data{cd}).LONGITUDE; offsets.(off_SNs{f}).LONGITUDE(temp_index)'];
                
                plotall.offsets.argo.(comp_data{cd}).val  = [plotall.offsets.argo.(comp_data{cd}).val ; (offsets.(off_SNs{f}).(comp_data{cd})(temp_index))'];
            end
        end
        if gdap_on==1
            if isfield(offsets.(off_SNs{f}),  'gdap')
                temp_index = ~isnan(offsets.(off_SNs{f}).gdap.(comp_data{cd}));
                plotall.offsets.gdap.(comp_data{cd}).LATITUDE = [plotall.offsets.gdap.(comp_data{cd}).LATITUDE; offsets.(off_SNs{f}).gdap.LATITUDE(temp_index)];
                plotall.offsets.gdap.(comp_data{cd}).LONGITUDE = [plotall.offsets.gdap.(comp_data{cd}).LONGITUDE ; offsets.(off_SNs{f}).gdap.LONGITUDE(temp_index)];
                
                plotall.offsets.gdap.(comp_data{cd}).val = [plotall.offsets.gdap.(comp_data{cd}).val; (offsets.(off_SNs{f}).gdap.(comp_data{cd})(temp_index))];
            end
        end
    end
end
%% plotting all crossover differences
argo_on = 0;
gdap_on = 1;
for cd = [3 4 5 8]
    
    if strcmp(comp_data{cd}, 'DOXY_ADJUSTED')
        c_lims = [0 10];
        units = '\mumol kg^-^1';
        var_name = 'Oxygen';
        color_map = brewermap(15, 'oranges');
    elseif strcmp(comp_data{cd}, 'NITRATE_ADJUSTED')
        c_lims = [0 1.5];
        units = '\mumol kg^-^1';
       var_name = 'Nitrate';
       color_map = brewermap(15, 'BuPu');
    elseif strcmp(comp_data{cd}, 'pH_25C_TOTAL')
        c_lims = [0 .03];
        units = 'pH';
        var_name = 'pH_2_5_C_-_t_o_t_a_l';
        color_map = brewermap(15, 'Greens');
    elseif strcmp(comp_data{cd}, 'DIC')
        c_lims = [0 10];
        units = '\mumol kg^-^1';
        var_name = 'DIC';
        color_map = brewermap(13, 'Blues');
    end
    color_map = color_map(4:end,:);

    plot_filename = ['GLOBAL Float Crossover Diffs ' comp_data{cd} ' argo=' num2str(argo_on) ' gdap=' num2str(gdap_on) ' dens match'];
    clf
    set(gcf, 'units', 'inches')
    paper_w = 12; paper_h =8;
    set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h
    
    
    subplot(1,1,1)
    hold on
    coast = load('coast');
    a1 = axesm('Mercator','MapLatLimit',[-75 75]);
    axis off; framem on; gridm off; mlabel on; plabel on;
    setm(gca,'MLabelParallel',75 , 'fontsize', 17)
    setm(gca, 'MLabelLocation', 60, 'PLabelLocation', 30)
%     a1.Color = [.8 .8 .8];
    
    coast_gray = [.4 .4 .4];
    geoshow(coast.lat,coast.long,'DisplayType','polygon', 'facecolor', coast_gray)
    marker_size = 4;
    
    % plotting background tracks

    % for f = 1:length(SNs)
    %     if ~isempty(offsets.(off_SNs{f}).(comp_data{cd}))
    
    %         geoshow(Argo.(SNs{f}).LATITUDE, Argo.(SNs{f}).LONGITUDE, 'displaytype', 'multipoint', 'marker', '.', 'markeredgecolor', 'k','markerfacecolor', 'k',  'markersize', marker_size)
    %     end
    % end
    if argo_on==1
        geoshow(plotall.offsets.argo.(comp_data{cd}).LATITUDE, plotall.offsets.argo.(comp_data{cd}).LONGITUDE, 'displaytype', 'multipoint', 'marker', '.', 'markeredgecolor', 'k','markerfacecolor', 'k',  'markersize', marker_size)
    end
    if gdap_on==1
        gdap_index = gdap_SO.PRES_ADJUSTED>1480 & gdap_SO.PRES_ADJUSTED<2010 & ~isnan(gdap_SO.(comp_data{cd}));
        
        geoshow(gdap_SO.LATITUDE(gdap_index), gdap_SO.LONGITUDE(gdap_index), 'displaytype', 'multipoint', 'marker', 's', 'markeredgecolor', 'k','markerfacecolor', 'k',  'markersize', 1.2)
    end
    %
    %
    % plotting crossover differences
    % %
    % for f = 1:length(off_SNs)
    %     if ~isempty(offsets.(off_SNs{f}).(comp_data{cd}))
    %         temp_index = ~isnan(offsets.(off_SNs{f}).(comp_data{cd}));
    %         scatterm(offsets.(off_SNs{f}).LATITUDE(temp_index), offsets.(off_SNs{f}).LONGITUDE(temp_index), 25, abs(offsets.(off_SNs{f}).(comp_data{cd})(temp_index)), 'filled', 'markeredgecolor', 'k');
    %
    %         if isfield(offsets.(off_SNs{f}),  'gdap')
    %             temp_index = ~isnan(offsets.(off_SNs{f}).gdap.(comp_data{cd}));
    %             scatterm(offsets.(off_SNs{f}).gdap.LATITUDE(temp_index), offsets.(off_SNs{f}).gdap.LONGITUDE(temp_index), 25, abs(offsets.(off_SNs{f}).gdap.(comp_data{cd})(temp_index)), 'filled', 'markeredgecolor', 'k','marker', 's');
    %         end
    %     end
    % end
    if argo_on==1
        scatterm(plotall.offsets.argo.(comp_data{cd}).LATITUDE, plotall.offsets.argo.(comp_data{cd}).LONGITUDE, 25, abs(plotall.offsets.argo.(comp_data{cd}).val), 'filled');
    end
    if gdap_on==1
        scatterm(plotall.offsets.gdap.(comp_data{cd}).LATITUDE, plotall.offsets.gdap.(comp_data{cd}).LONGITUDE, 25, abs(plotall.offsets.gdap.(comp_data{cd}).val), 'filled', 'marker', 's');
    end
    

    
    c1 = colorbar;
    ylabel(c1, units)
    set(c1, 'fontsize', 20)
    caxis(c_lims)
    set(gcf, 'colormap', color_map)
    title(var_name, 'fontsize', 22)
    
    print(gcf, '-dpdf', '-r400', [Plot_dir '../' plot_filename '.pdf' ])
end
%% Histograms of differences
plot_filename = 'Global mean GLODAP differences';

clf
set(gcf, 'units', 'inches')
paper_w = 9; paper_h =4;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h


cd = 3; cutoff = 50;
d1 = subplot(1,3,1);

histogram(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff),50)
grid on; hold on

val_mean = nanmean(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));
xlabel(['Mean: ' num2str(val_mean,2)])
title(comp_data{cd}, 'interpreter', 'none')
val_std = nanstd(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));
disp([comp_data{cd} ' Std: ' num2str(val_std,2)])

x = -cutoff:.1:cutoff;
y_norm = normpdf(x,0,val_std);
num_pts = numel(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));

plot(d1, x, y_norm*num_pts*3, 'r', 'linewidth', 2)
%
cd = 4; cutoff = 6;
d2 = subplot(1,3,2);

histogram(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff), 50)
grid on; hold on

val_mean = nanmean(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));
xlabel(['Mean: ' num2str(val_mean,2)])
title(comp_data{cd}, 'interpreter', 'none')
val_std = nanstd(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));
disp([comp_data{cd} ' Std: ' num2str(val_std,2)])

x = -cutoff:.1:cutoff;
y_norm = normpdf(x,0,val_std);
num_pts = numel(plotall.offsets.gdap.(comp_data{cd}).val(abs(plotall.offsets.gdap.(comp_data{cd}).val)<cutoff));

plot(d2, x, y_norm*num_pts*.3, 'r', 'linewidth', 2)

print(gcf, '-dpng', '-r400', [Plot_dir '../' plot_filename '.png' ])


%% Looking for regions with lots of data for long-term changes

lat_range = -70:2:70;
lon_range = 0:2:360;

[XX, YY] = meshgrid(lon_range, lat_range);

grid_count.gdap = zeros(length(lon_range), length(lat_range));
grid_count.argo = zeros(length(lon_range), length(lat_range));

for lo = 1:length(lon_range)-1
    disp(lon_range(lo))
    for la = 1:length(lat_range)-1
        
        for f = 1:length(SNs)
            lat_index = Argo.(SNs{f}).LATITUDE>=lat_range(la) & Argo.(SNs{f}).LATITUDE<lat_range(la+1);
            lon_index = Argo.(SNs{f}).LONGITUDE>=lon_range(lo) & Argo.(SNs{f}).LONGITUDE<lon_range(lo+1);
            
            argo_index = find(lat_index & lon_index);
            if isempty(argo_index)
                continue
            else
                grid_count.argo(lo, la) = grid_count.argo(lo, la) + length(argo_index);
            end
        end
        
        gdap_index = gdap_SO.LATITUDE>=lat_range(la) & gdap_SO.LATITUDE<lat_range(la+1) & gdap_SO.LONGITUDE>=lon_range(lo) & gdap_SO.LONGITUDE<lon_range(lo+1) & gdap_SO.PRES_ADJUSTED<30;
        
        grid_count.gdap(lo, la) = grid_count.gdap(lo, la) + sum(gdap_index);
    end
end

%%

grid_count.both =  zeros(length(lon_range), length(lat_range));

grid_count.both(grid_count.argo>10) = grid_count.both(grid_count.argo>10) + 1;

grid_count.both(grid_count.gdap>10) = grid_count.both(grid_count.gdap>10) + 1;

%%
coast = load('coast');

coast.long_360 = coast.long;
coast.long_360(coast.long_360<0) = coast.long_360(coast.long_360<0)+360;
clear regions

regions.SEPac.LONGITUDE =  [254 260];
regions.SEPac.LATITUDE = [-40 -28];
regions.GOA.LONGITUDE = [213 226] ;
regions.GOA.LATITUDE = [45 56] ;
regions.TAS.LONGITUDE = [138 149];
regions.TAS.LATITUDE = [-51 -44]; ...
regions.KC.LONGITUDE = [131 138];
regions.KC.LATITUDE = [20 30];
regions.SAtl.LONGITUDE = [318 326];
regions.SAtl.LATITUDE = [-55 -45];

RRs = fieldnames(regions);

clf
d1 = subplot(3,1,1);
pcolor(XX,YY, grid_count.argo'); shading flat; colorbar
hold on
plot(coast.long_360, coast.lat, 'k.')
caxis([0 100])

d2 = subplot(3,1,2);
pcolor(XX,YY, grid_count.gdap'); shading flat; colorbar
hold on
plot(coast.long_360, coast.lat, 'k.')
caxis([0 100])
set(gcf, 'colormap', brewermap(30, 'Reds'))

d3 = subplot(3,1,3);
pcolor(XX,YY, grid_count.both'); shading flat; colorbar
hold on
plot(coast.long_360, coast.lat, 'k.')
% caxis([0 200])
set(gcf, 'colormap', brewermap(30, 'Reds'))

for r= 1:size(RRs,1)
    plot(d1, [regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(1)], ...
        [regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(1)], 'k', 'linewidth', 2)
    
    plot(d2, [regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(1)], ...
        [regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(1)], 'k', 'linewidth', 2)
    
    plot(d3, [regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(2) regions.(RRs{r}).LONGITUDE(1) regions.(RRs{r}).LONGITUDE(1)], ...
        [regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(1) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(2) regions.(RRs{r}).LATITUDE(1)], 'k', 'linewidth', 2)
    
end
%% find data within geographic ranges, collect all info
clear r_prop
%%

for r = 3%1:length(RRs)
    
    for m = 1:length(meta_data)
        r_prop.(RRs{r}).argo.(meta_data{m}) = [];
    end
    
    for c = 1:length(comp_data)
        r_prop.(RRs{r}).argo.(comp_data{c}) = [];
    end
    
    
    
    for q = 1:length(SNs)
        if q/15==round(q/15)
            disp([num2str(q/length(SNs)*100,3) '% done'])
        end
        prof_index = find(Argo.(SNs{q}).LATITUDE>= regions.(RRs{r}).LATITUDE(1) & Argo.(SNs{q}).LATITUDE<= regions.(RRs{r}).LATITUDE(2) & ...
            Argo.(SNs{q}).LONGITUDE>= regions.(RRs{r}).LONGITUDE(1) & Argo.(SNs{q}).LONGITUDE<= regions.(RRs{r}).LONGITUDE(2));
        
        if isempty(prof_index)
            continue
        end
        for p = 1:length(prof_index)
            
            for m = 1:length(meta_data) % add metadata vectors that are the length of the profile data
                r_prop.(RRs{r}).argo.(meta_data{m}) = [r_prop.(RRs{r}).argo.(meta_data{m}) ; Argo.(SNs{q}).(meta_data{m})(prof_index(p)).*ones(size(Argo.(SNs{q}).PRES_ADJUSTED,2),1)];
            end
            
            for m = 1:length(comp_data)
                if isfield(Argo.(SNs{q}), comp_data{m})
                    r_prop.(RRs{r}).argo.(comp_data{m}) = [r_prop.(RRs{r}).argo.(comp_data{m}) ; Argo.(SNs{q}).(comp_data{m})(prof_index(p),:)'];
                else
                    r_prop.(RRs{r}).argo.(comp_data{m}) = [r_prop.(RRs{r}).argo.(comp_data{m}) ; NaN(size(Argo.(SNs{q}).PRES_ADJUSTED,2),1)];
                end
            end
        end
    end
    
    gdap_index = gdap_SO.LATITUDE>= regions.(RRs{r}).LATITUDE(1) & gdap_SO.LATITUDE<= regions.(RRs{r}).LATITUDE(2) & ...
        gdap_SO.LONGITUDE>= regions.(RRs{r}).LONGITUDE(1) & gdap_SO.LONGITUDE<= regions.(RRs{r}).LONGITUDE(2);
    
    for m = 1:length(meta_data)
        r_prop.(RRs{r}).gdap.(meta_data{m}) = gdap_SO.(meta_data{m})(gdap_index);
    end
    
    for c = 1:length(comp_data)
        r_prop.(RRs{r}).gdap.(comp_data{c}) =  gdap_SO.(comp_data{c})(gdap_index);
    end
    
    
end
%% monthly means
DS = {'argo';'gdap'};

r = 3;
dens_bin = [27.34 27.35];
t = 2;

date_vec = datevec(r_prop.(RRs{r}).(DS{t}).GMT_Matlab);

start_year = min(date_vec(:,1));
end_year = max(date_vec(:,1));

monthly_date = datenum(start_year,1:12.*(end_year - start_year+1),15);

r_prop.(RRs{r}).(DS{t}).avg.GMT_Matlab = monthly_date';

for m = 2:length(meta_data)
    r_prop.(RRs{r}).(DS{t}).avg.(meta_data{m}) = NaN(length(monthly_date),2);
end

for c = 1:length(comp_data)
    r_prop.(RRs{r}).(DS{t}).avg.(comp_data{c}) = NaN(length(monthly_date),2);
end

for d = 1:length( monthly_date)
    mon_vec = datevec(monthly_date(d));
    index = r_prop.(RRs{r}).(DS{t}).PDENS-1000>= dens_bin(1) & r_prop.(RRs{r}).(DS{t}).PDENS-1000<= dens_bin(2) &  ...
        r_prop.(RRs{r}).(DS{t}).GMT_Matlab>= datenum(mon_vec(1), mon_vec(2), 1) &  r_prop.(RRs{r}).(DS{t}).GMT_Matlab<= datenum(mon_vec(1), mon_vec(2)+1, 1);
    
    for m = 2:length(meta_data)
        
        r_prop.(RRs{r}).(DS{t}).avg.(meta_data{m})(d,1) =  nanmean(r_prop.(RRs{r}).(DS{t}).(meta_data{m})(index));
        r_prop.(RRs{r}).(DS{t}).avg.(meta_data{m})(d,2) =  nanstd(r_prop.(RRs{r}).(DS{t}).(meta_data{m})(index));
        
    end
    
    for m = 1:length(comp_data)
        
        r_prop.(RRs{r}).(DS{t}).avg.(comp_data{m})(d,1) =  nanmean(r_prop.(RRs{r}).(DS{t}).(comp_data{m})(index));
        r_prop.(RRs{r}).(DS{t}).avg.(comp_data{m})(d,2) =  nanstd(r_prop.(RRs{r}).(DS{t}).(comp_data{m})(index));
        
    end
end

%% Calculate SPICE

for r = 1:length(RRs)
    for t = 1:2
       r_prop.(RRs{r}).(DS{t}).SPICE = spice(r_prop.(RRs{r}).(DS{t}).TEMP_ADJUSTED, r_prop.(RRs{r}).(DS{t}).PSAL_ADJUSTED);
        
    
    end
end
    %%
clf
plot(r_prop.SEPac.argo.LONGITUDE, r_prop.SEPac.argo.LATITUDE, 'x')
hold on
plot(r_prop.SEPac.gdap.LONGITUDE, r_prop.SEPac.gdap.LATITUDE, 'or')

%%
% pres_bin = [300 1000];
r =2;
dens_bin = 27.3:.01:27.55; % 10 m is ~.01 kg/m3
DS = {'argo';'gdap'};
symbols = {'s', '^'};
marker_edge_color = {[1 1 1];[0 0 0]};
for p = 1:length(dens_bin)-1
    clf
    d1 = subplot(4,2,2); hold on; grid on; title('Temp')
    d2 = subplot(4,2,4);hold on; grid on; title('Sal')
    d3 = subplot(4,2,6);hold on; grid on; title('O2')
    d4 = subplot(4,2,[1 3]); hold on; grid on
    
    d5 = subplot(4,2,7); hold on; grid on; title('NO3')
    d6 = subplot(4,2,8); hold on; grid on; title('pH 25C Total')
    d7 = subplot(4,2,5); hold on; grid on; title('DIC')

    for t = 1:length(DS)
        %     index = r_prop.(RRs{r}).(DS{t}).PRES_ADJUSTED>= pres_bin(1) & r_prop.(RRs{r}).(DS{t}).PRES_ADJUSTED<= pres_bin(2);
        index = r_prop.(RRs{r}).(DS{t}).PDENS-1000>= dens_bin(p) & r_prop.(RRs{r}).(DS{t}).PDENS-1000<= dens_bin(p+1);
        plot(d1, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).TEMP_ADJUSTED(index), 'marker', symbols{t}, 'linestyle', 'none')
        plot(d2, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).PSAL_ADJUSTED(index), 'marker', symbols{t}, 'linestyle', 'none')
        
        scatter(d3, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).DOXY_ADJUSTED(index), 100,r_prop.(RRs{r}).(DS{t}).SPICE(index),'filled', 'marker', symbols{t}); colorbar(d3)
        plot(d5, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).NITRATE_ADJUSTED(index), 'marker', symbols{t}, 'linestyle', 'none')
        plot(d6, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).pH_25C_TOTAL(index), 'marker', symbols{t}, 'linestyle', 'none')
        plot(d7, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).DIC(index), 'marker', symbols{t}, 'linestyle', 'none')

        scatter(d4, r_prop.(RRs{r}).(DS{t}).LONGITUDE(index), r_prop.(RRs{r}).(DS{t}).LATITUDE(index), 80, r_prop.(RRs{r}).(DS{t}).PRES_ADJUSTED(index), 'filled', 'marker', symbols{t}, 'markeredgecolor', marker_edge_color{t})
        
    end
    colorbar(d4)
    title(d4, [RRs{r} ' ' num2str(dens_bin(p)) ' ' num2str(dens_bin(p+1))])
    datetick(d1, 'x', 'YY-mm')
    datetick(d2, 'x', 'YY-mm')
    datetick(d3, 'x', 'YY-mm')
    datetick(d5, 'x', 'YY-mm')
    datetick(d6, 'x', 'YY-mm')
    datetick(d7, 'x', 'YY-mm')

    pause
end
%%
colors = brewermap(3,'Set1');
r =3;
dens_bin = [27.34 27.35];


plot_filename = [RRs{r} '_example_BGC_Change_' num2str(dens_bin(1)) '_to_' num2str(dens_bin(2))];


clf
set(gcf, 'units', 'inches')
paper_w = 9; paper_h =7;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h


d3 = subplot(3,1,1);hold on; grid on; title('O2')
d5 = subplot(3,1,2); hold on; grid on; title('NO3')
d7 = subplot(3,1,3); hold on; grid on; title('DIC')
p=1;
marker_size = 8;

l_d = NaN(2,1);
for t = 1:length(DS)
    
    index = r_prop.(RRs{r}).(DS{t}).PDENS-1000>= dens_bin(p) & r_prop.(RRs{r}).(DS{t}).PDENS-1000<= dens_bin(p+1);
    
    plot(d3, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).DOXY_ADJUSTED(index), 'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size)
    plot(d5, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).NITRATE_ADJUSTED(index), 'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size)
    l_d(t) =  plot(d7, r_prop.(RRs{r}).(DS{t}).GMT_Matlab(index),  r_prop.(RRs{r}).(DS{t}).DIC(index), 'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size);
    
end

set(d3, 'xtick', datenum(1990:5:2020,1,1))
set(d5, 'xtick', datenum(1990:5:2020,1,1))
set(d7, 'xtick', datenum(1990:5:2020,1,1))


legend(l_d, 'Float', 'GLODAP', 'location', 'northwest')
    datetick(d3, 'x', 'YY', 'keepticks')
    datetick(d5, 'x', 'YY', 'keepticks')
    datetick(d7, 'x', 'YY', 'keepticks')
ylabel(d3, '\mumol kg^-^1')
    ylabel(d5, '\mumol kg^-^1')
ylabel(d7, '\mumol kg^-^1')

print(gcf, '-dpng', '-r400', [Plot_dir '../' plot_filename '.png' ])


%% monthly means
colors = brewermap(3,'Set1');
r =3;
dens_bin = [27.34 27.35];


plot_filename = [RRs{r} '_example_BGC_Change_monthly' num2str(dens_bin(1)) '_to_' num2str(dens_bin(2))];


clf
set(gcf, 'units', 'inches')
paper_w = 9; paper_h =7;
set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]); clear paper_w paper_h


d3 = subplot(3,1,1);hold on; grid on; title('O2')
d5 = subplot(3,1,2); hold on; grid on; title('NO3')
d7 = subplot(3,1,3); hold on; grid on; title('DIC')
p=1;
marker_size = 8;

l_d = NaN(2,1);
for t = 1:length(DS)
    
%     index = r_prop.(RRs{r}).(DS{t}).PDENS-1000>= dens_bin(p) & r_prop.(RRs{r}).(DS{t}).PDENS-1000<= dens_bin(p+1);
    
    errorbar(d3, r_prop.(RRs{r}).(DS{t}).avg.GMT_Matlab,  r_prop.(RRs{r}).(DS{t}).avg.DOXY_ADJUSTED(:,1), r_prop.(RRs{r}).(DS{t}).avg.DOXY_ADJUSTED(:,2),'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size)
    errorbar(d5, r_prop.(RRs{r}).(DS{t}).avg.GMT_Matlab,  r_prop.(RRs{r}).(DS{t}).avg.NITRATE_ADJUSTED(:,1), r_prop.(RRs{r}).(DS{t}).avg.NITRATE_ADJUSTED(:,2),'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size)
    l_d(t) =  errorbar(d7, r_prop.(RRs{r}).(DS{t}).avg.GMT_Matlab,  r_prop.(RRs{r}).(DS{t}).avg.DIC(:,1), r_prop.(RRs{r}).(DS{t}).avg.DIC(:,2),'marker', symbols{t}, 'linestyle', 'none', 'markerfacecolor', colors(t+1,:), 'color', colors(t+1,:), 'markersize', marker_size);
    
end

set(d3, 'xtick', datenum(1990:5:2020,1,1))
set(d5, 'xtick', datenum(1990:5:2020,1,1))
set(d7, 'xtick', datenum(1990:5:2020,1,1))


legend(l_d, 'Float', 'GLODAP', 'location', 'northwest')
    datetick(d3, 'x', 'YY', 'keepticks')
    datetick(d5, 'x', 'YY', 'keepticks')
    datetick(d7, 'x', 'YY', 'keepticks')
ylabel(d3, '\mumol kg^-^1', 'fontsize', 16)
    ylabel(d5, '\mumol kg^-^1','fontsize', 16)
    ylabel(d7, '\mumol kg^-^1', 'fontsize',16)

print(gcf, '-dpng', '-r400', [Plot_dir '../' plot_filename '.png' ])
print(gcf, '-dpdf', '-r400', [Plot_dir '../' plot_filename '.pdf' ])


%% for each float, plot the mean o2 offset for glodap vs the mean NO3 offset
clf
subplot(1,1,1); hold on; grid on
for f = 1:length(off_SNs)
    if isfield(offsets.(off_SNs{f}), 'gdap')
        if isfield(offsets.(off_SNs{f}).gdap, 'DOXY_ADJUSTED') && isfield(offsets.(off_SNs{f}).gdap, 'NITRATE_ADJUSTED')
            plot(nanmean(offsets.(off_SNs{f}).gdap.DOXY_ADJUSTED), nanmean(offsets.(off_SNs{f}).gdap.NITRATE_ADJUSTED), '.k')
        end
    end
end