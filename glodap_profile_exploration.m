gdap = load([data_dir 'Data_Products/GLODAP/GLODAPv2.2022_Merged_Master_File.mat']);
%%
gdap.PDENS = sw_pden(gdap.G2salinity,gdap.G2temperature,gdap.G2pressure, 0);
gdap.GMT_Matlab = datenum(gdap.G2year, gdap.G2month, gdap.G2day, gdap.G2hour, gdap.G2minute, zeros(size(gdap.G2minute)));


gdap_fields_no_flags = {'G2longitude', 'G2latitude', 'GMT_Matlab', 'G2year', 'G2month', 'G2cruise', 'G2station', 'G2pressure', 'G2temperature', 'PDENS'};
gdap_fields_w_flags = {'G2salinity','G2oxygen', 'G2nitrate', 'G2tco2', 'G2talk', 'G2phts25p0', 'G2aou'};

gdap.G2longitude(gdap.G2longitude<0) = gdap.G2longitude(gdap.G2longitude<0)+360;


% 
% for g = 1:length(gdap_fields_no_flags)
%     gdap_SO.(gdap_fields_no_flags{g}) = gdap.(gdap_fields_no_flags{g})(SO_gdap_index);
% end
% clear g

for g = 1:length(gdap_fields_w_flags)
%     temp_data = gdap.(gdap_fields_w_flags{g});
    gdap.(gdap_fields_w_flags{g})(gdap.([gdap_fields_w_flags{g} 'f'])~=2) = nan;
    
%     gdap_SO.(gdap_fields_w_flags{g}) =temp_data(SO_gdap_index);

end
clear g
%% 

cruise_index = gdap.G2cruise==233;
profile_range = 5:29;
% cruise_index = gdap.G2cruise==2028;
% profile_range = 2825:2830;

cruise_1 = [double(gdap.G2cruise(cruise_index)) gdap.G2station(cruise_index) gdap.G2longitude(cruise_index) gdap.G2latitude(cruise_index) ...
    gdap.G2depth(cruise_index) gdap.PDENS(cruise_index) gdap.G2oxygen(cruise_index) gdap.G2aou(cruise_index) ] ;

%

clf
d1 = subplot(2,4,1); hold on;  grid on
d2 = subplot(2,4,2); hold on; set(gca, 'ydir', 'reverse'); grid on
d3 = subplot(2,4,3); hold on; grid on
d4 = subplot(2,4,4); hold on; grid on; set(gca, 'ydir', 'reverse')
d5 = subplot(2,4,5); hold on; grid on; set(gca, 'ydir', 'reverse')
d6 = subplot(2,4,6); hold on; grid on; set(gca, 'ydir', 'reverse')
d7 = subplot(2,4,7); hold on; grid on; set(gca, 'ydir', 'reverse')
d8 = subplot(2,4,8); hold on; grid on; set(gca, 'ydir', 'reverse')



for p = 1:length(profile_range)
    profile = cruise_1(cruise_1(:,2)==profile_range(p),:);
    plot(d1, profile(:,3), profile(:,4), 'x', 'linewidth', 2)

    plot(d2, profile(:,7), profile(:,5), 'linewidth', 2)

    plot(d3, profile(:,7), profile(:,6), 'linewidth', 2)

    plot(d4, profile(:,8), profile(:,5), 'linewidth', 2)
    imagesc(d5, profile(:,3), profile(:,5), profile(:,7));

    imagesc(d7, profile(:,3), profile(:,6), profile(:,7));

end

cruise_index = gdap.G2cruise==1110;
profile_range = 144:170;

% cruise_index = gdap.G2cruise==2031;
% profile_range = 3150:3175;

cruise_1 = [double(gdap.G2cruise(cruise_index)) gdap.G2station(cruise_index) gdap.G2longitude(cruise_index) gdap.G2latitude(cruise_index) ...
    gdap.G2depth(cruise_index) gdap.PDENS(cruise_index) gdap.G2oxygen(cruise_index)  gdap.G2aou(cruise_index) ] ;

for p = 1:length(profile_range)
    profile = cruise_1(cruise_1(:,2)==profile_range(p),:);
    plot(d1, profile(:,3), profile(:,4), 'o', 'linewidth', 2)

    plot(d2, profile(:,7), profile(:,5), 'linewidth', 2, 'marker', 'o')

    plot(d3, profile(:,7), profile(:,6), 'linewidth', 2, 'marker', 'o')

    plot(d4, profile(:,8), profile(:,5), 'linewidth', 2, 'marker', 'o')
    imagesc(d6, profile(:,3), profile(:,5), profile(:,7));
    imagesc(d8, profile(:,3), profile(:,6), profile(:,7));

end
colorbar(d5)
colorbar(d6)
colorbar(d7)
colorbar(d8)
o2_lims = [190 220];
clim(d5, o2_lims)
clim(d6, o2_lims)
clim(d7, o2_lims)
clim(d8, o2_lims)