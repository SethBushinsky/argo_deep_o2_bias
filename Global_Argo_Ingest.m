%% read in float profile names from JCOMMOPS and download data
% 
% float_dir = [home_dir 'Data/ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2020_11_16/'];
% % float_dir = [home_dir 'Data/ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2022_03_10/'];
% 
% % working_dir = [home_dir 'Work/Proposals/2020_11 NOAA GOMO Float BGC Change Synthesis/Data/'];
% 
% files = {'Ocean_Ops_O2_list2', 'Ocean_Ops_pH_list2', 'Ocean_Ops_NO3_list2'};
% 
% % platform_list = 'Platforms_JCOMMOPS_KE_Region';
% 
% opts = delimitedTextImportOptions("NumVariables", 11);
% 
% % Specify range and delimiter
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% 
% % Specify column names and types
% opts.VariableNames = ["COUNTRY", "DEPLOYMENTDATE", "DEPLOYMENTLAT", "DEPLOYMENTLON", "LASTLOCATIONDATE", "MODEL", "OBSERVINGNETWORKS", "PROGRAM", "REF", "SERIALNUMBER", "STATUS"];
% opts.VariableTypes = ["char", "char", "double", "double", "char", "char", "char", "char", "double", "double", "char"];
% opts = setvaropts(opts, [1, 2, 5, 6, 7, 8, 11], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, [1, 2, 5, 6, 7, 8, 11], "EmptyFieldRule", "auto");
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Import the data
% for f = 1:length(files)
%     oceanops.(files{f}) = readtable([float_dir files{f} '.csv'], opts);
%     
%     % Convert to output type
%     oceanops.(files{f}) = table2cell(oceanops.(files{f}));
%     numIdx = cellfun(@(x) ~isnan(str2double(x)), oceanops.(files{f}));
%     oceanops.(files{f})(numIdx) = cellfun(@(x) {str2double(x)}, oceanops.(files{f})(numIdx));
% end
% % Clear temporary variables
% clear opts numIdx

%%
clear oceanops

float_dir = [home_dir 'Data/ARGO_O2_Floats/Global/GLOBAL_BGC_ARGO/2022_03_10/'];
files = {'Ocean_Ops_O2_list2', 'Ocean_Ops_pH_list2', 'Ocean_Ops_NO3_list2'};

for f = 1:length(files)
oceanops.(files{f}) = readtable([float_dir files{f} '.csv']);
end
%%  generate a list of all SNs
WMO_list_all = [];
wmo_index = 0;
for f = 1:length(files)
    for q = 1:height(oceanops.(files{f}))
%         if  strcmp(oceanops.(files{f}){q,7}, 'Argo')
            wmo_index = wmo_index+1;
            %         WMO_list_all(wmo_index,1) = oceanops.(files{f}){q,9};
            if iscell(oceanops.(files{f}).REF(q))
                WMO_list_all(wmo_index,1) =  str2double(oceanops.(files{f}).REF{q});
            else
                WMO_list_all(wmo_index,1) = oceanops.(files{f}).REF(q);
                
            end
%         end
    end
end

clear q wmo_index
WMO_list = unique(WMO_list_all);

% contains the WMO numbers of all NO3, pH, and O2 floats, globally
WMO_list = WMO_list(~isnan(WMO_list));

clear WMO_list_all f files oceanops
%% find the correct DAC location for each float

% load argo metadata list:
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [10, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["file", "profiler_type", "institution", "date_update"];
opts.VariableTypes = ["char", "double", "char", "double"];
opts = setvaropts(opts, [1, 3], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 3], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
arindexglobalmeta = readtable([float_dir 'ar_index_global_meta.txt'], opts);

% Convert to output type
arindexglobalmeta = table2cell(arindexglobalmeta);
numIdx = cellfun(@(x) ~isnan(str2double(x)), arindexglobalmeta);
arindexglobalmeta(numIdx) = cellfun(@(x) {str2double(x)}, arindexglobalmeta(numIdx));

% Clear temporary variables
clear opts numIdx

%% generate a reference list of SNs from the metadata search
clear all_argo_WMO
all_argo_WMO = NaN(length(arindexglobalmeta),1);

for m = 1:length(arindexglobalmeta)
    dir_breaks = strfind(arindexglobalmeta{m,1}, '/');
    
    all_argo_WMO(m,1) = str2double(arindexglobalmeta{m,1}(dir_breaks(1)+1:dir_breaks(2)-1));
    clear dir_breaks
    
end
clear m
%% download SProf files
% float_dir = [home_dir 'Data/ARGO_O2_Floats/Global/ALL_BGC_ARGO/2020_10_23/netcdf/'];

setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

flag_out = cell(length(WMO_list),2);
last = 1;
%%
for q = last:length(WMO_list)
    if flag_out{q,1}==0
        continue
    end
    
    float_index = all_argo_WMO==WMO_list(q,1);
    
    if sum(float_index)==0
        disp(['No match for ' num2str(WMO_list(q,1))])
        flag_out{q,2} = ['No match for ' num2str(WMO_list(q,1))];
        continue
    end
    
    dir_breaks = strfind(arindexglobalmeta{float_index,1}, '/');
    
    partial_dir = arindexglobalmeta{float_index,1}(1:dir_breaks(2));
    
    file_dir = ['ftp://ftp.ifremer.fr/ifremer/argo/dac/' partial_dir num2str(WMO_list(q,1)) '_Sprof.nc'];
    
    flag_out{q,1} = system(['wget -P ' float_dir ' ' file_dir]);
    last = last+1;
    
    flag_out{q,2} = partial_dir;
    clear partial_dir dir_breaks float_index file_dir
    
end

clear float_index q
%% save out a text file of only the floats that failed to download or were not found
clear last

failed_floats = {};
for q = 1:length(flag_out)
    if isempty(flag_out{q,1}) || flag_out{q,1}~=0 
        failed_floats{end+1,1} = flag_out{q,2};
    end
    
end
%% reading in and plotting data
clear Argo
parameter_list = {'REFERENCE_DATE_TIME', 'PLATFORM_NUMBER', 'STATION_PARAMETERS', 'DIRECTION', 'PARAMETER_DATA_MODE','JULD',  'LONGITUDE', 'LATITUDE', 'POSITION_QC', 'PARAMETER', ...
    'PRES_ADJUSTED', 'PRES_ADJUSTED_QC', 'TEMP_ADJUSTED', 'TEMP_ADJUSTED_QC', 'PSAL_ADJUSTED', 'PSAL_ADJUSTED_QC', 'DOXY_ADJUSTED', 'DOXY_ADJUSTED_QC', 'DOXY_ADJUSTED_ERROR', 'CHLA_ADJUSTED', ...
    'CHLA_ADJUSTED_QC', 'BBP700_ADJUSTED', 'BBP700_ADJUSTED_QC', 'BBP700_ADJUSTED_ERROR', 'PH_IN_SITU_TOTAL', 'PH_IN_SITU_TOTAL_ADJUSTED', 'PH_INSITU_TOTAL_ADJUSTED_QC', 'PH_IN_SITU_TOTAL_ADJUSTED_ERROR', ...
    'NITRATE_ADJUSTED', 'NITRATE_ADJUSTED_QC', 'NITRATE_ADJUSTED_ERROR'};

float_files = dir([float_dir 'Sprof/*.nc']);

for f = 1:length(float_files)
    disp(['Starting ' num2str(f) ' out of ' num2str(length(float_files))])
uscore_index = strfind(float_files(f).name, '_');

Argo.(['f' float_files(f).name(1:uscore_index(1)-1)]).WMO = str2double(float_files(f).name(1:uscore_index(1)-1));

for p1 = 1:length(parameter_list)
    try
    Argo.(['f' float_files(f).name(1:uscore_index(1)-1)]).(parameter_list{p1}) = ncread([float_dir 'Sprof/' float_files(f).name], parameter_list{p1});
    catch
    end
end
clear uscore_index p1
end

clear float_files f

SNs = fieldnames(Argo);

%% Derived quantities / adjustments to dataset %% calculate Matlab time

for f = 1:length(SNs)
    Argo.(SNs{f}).GMT_Matlab = Argo.(SNs{f}).JULD+ datenum(Argo.(SNs{f}).REFERENCE_DATE_TIME', 'YYYYmmddHHMMSS');

    float_fields = fieldnames(Argo.(SNs{f}));
    
    num_profiles = length(Argo.(SNs{f}).GMT_Matlab);
    
    % transpose fields
    for z = 1:length(float_fields)
        
        if size(Argo.(SNs{f}).(float_fields{z}), 2)==num_profiles && ~ischar(Argo.(SNs{f}).(float_fields{z}))

            Argo.(SNs{f}).(float_fields{z}) = Argo.(SNs{f}).(float_fields{z})';
        end
    end
    
    % calculate potential density
    Argo.(SNs{f}).PDENS = sw_pden(Argo.(SNs{f}).PSAL_ADJUSTED, Argo.(SNs{f}).TEMP_ADJUSTED, Argo.(SNs{f}).PRES_ADJUSTED, 1);
    
    % sets lon range to 0-360
    Argo.(SNs{f}).LONGITUDE(Argo.(SNs{f}).LONGITUDE<0) = Argo.(SNs{f}).LONGITUDE(Argo.(SNs{f}).LONGITUDE<0) + 360;
    
    if f/10==round(f/10)
        disp([num2str(f/length(SNs)*100) ' percent done'])
    end
end

%% save dataset
save([float_dir '../Argo_BGC_' datestr(now, 'YYYY-mm-dd')], 'Argo', 'SNs', 'float_dir',  '-V7.3')

%% Calculate a surface / 25 m avg.
pres_lim = 25;
temp_list = {'TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED'};

for f = 1:length(SNs)
    for t = 1:length(temp_list)
        if ~isfield(Argo.(SNs{f}), temp_list{t})
            continue
        end
        Argo.(SNs{f}).([temp_list{t} '_' num2str(pres_lim)]) = NaN(length(Argo.(SNs{f}).GMT_Matlab),1);
        for p = 1:length(Argo.(SNs{f}).GMT_Matlab)
            
            
            Argo.(SNs{f}).([temp_list{t} '_' num2str(pres_lim)])(p,1) = nanmean(Argo.(SNs{f}).(temp_list{t})(Argo.(SNs{f}).PRES_ADJUSTED(:,p)<=pres_lim,p));
            
        end
    end
    
    Argo.(SNs{f}).PDEN_25=  sw_pden(Argo.(SNs{f}).PSAL_ADJUSTED_25, Argo.(SNs{f}).TEMP_ADJUSTED_25,  12, 0);
end
%% to add:

% calculate MLD / all other derived quantities currently in your code
% this would make the dataset entirely too large I think..