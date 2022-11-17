% script for comparing offsets between original Matlab code run for the
% proposal and new Python code
% first run Float_Crossover_Comparison_v3.m making sure to use the same
% source Argo files and Glodap
% Matlab offsets stored in "temp_Matlab_glodap_offsets.mat" so you don't
% have to re-run every time.

%% Load Python output glodap_offsets.nc

project_dir = [home_dir 'Work/Projects/2021_07_Float_BGC_QC_NOAA/code/'];

gdap_offsets.info = ncinfo([project_dir 'output/glodap_offsets.nc']);

for v = 1:length(gdap_offsets.info.Variables)
    gdap_offsets.(gdap_offsets.info.Variables(v).Name) = ncread([project_dir 'output/glodap_offsets.nc'], ...
        gdap_offsets.info.Variables(v).Name);
end

%% Search by WMO and compare offsets w/ depth

list_python_wmos = unique(gdap_offsets.main_float_wmo);

list_matlab_wmos = fieldnames(offsets);
matlab_wmo_ns = NaN(length(list_matlab_wmos),1);

% strip the 'f' from each matlab wmo and make into double
for f = 1:length(list_matlab_wmos)
   matlab_wmo_ns(f) = str2double(list_matlab_wmos{f}(2:end));
    
end

%%
c_map = brewermap(6, 'Set1');

Plot_dir = [project_dir '../plots_matlab/offset_comparison/'];
for w = 1:length(list_python_wmos)
    wmo_search = double(list_python_wmos(w));
    python_index = gdap_offsets.main_float_wmo==wmo_search;
    
    matlab_index = matlab_wmo_ns==wmo_search;
    
    if sum(matlab_index)>0 &&  ...
            ~isempty(offsets.(list_matlab_wmos{matlab_index}).gdap.GMT_Matlab_float)
        % histograms of offsets
        clf
        set(gcf, 'units', 'inches')
        paper_w = 7; paper_h =10;
        set(gcf,'PaperSize',[paper_w paper_h],'PaperPosition', [0 0 paper_w paper_h]);
        
        for c = 1:4% [1:6 8]
            if isempty(offsets.(list_matlab_wmos{matlab_index}).gdap.(comp_data{c}))
                continue
            end
            subplot(6,2,1+(c-1)*2); hold on
            title( comp_data{c}, 'interpreter', 'none')

            plot(offsets.(list_matlab_wmos{matlab_index}).gdap.([comp_data{c} '_offset']), ...
                offsets.(list_matlab_wmos{matlab_index}).gdap.PRES_ADJUSTED, 'x', 'color', c_map(2,:))
            hold on
            plot(gdap_offsets.([comp_data{c} '_offset'])(python_index), ...
                gdap_offsets.PRES_ADJUSTED(python_index), 'o', 'color', c_map(4,:))
            
            set(gca, 'ydir', 'reverse')
            
            subplot(5,2,2+(c-1)*2); hold on
            h1 = histogram(offsets.(list_matlab_wmos{matlab_index}).gdap.([comp_data{c} '_offset']));
            
            h2 = histogram(gdap_offsets.([comp_data{c} '_offset'])(python_index));
            
            h1.FaceColor = c_map(2,:);
            h2.FaceColor = c_map(4,:);

            h1.BinEdges = h2.BinEdges;
            
            matlab_leg = ['Matlab: ' num2str(nanmean(offsets.(list_matlab_wmos{matlab_index}).gdap.([comp_data{c} '_offset'])),3) ...
                ' \pm ' num2str(nanstd(offsets.(list_matlab_wmos{matlab_index}).gdap.([comp_data{c} '_offset'])),3)];
            
            python_leg = ['Python: ' num2str(nanmean(gdap_offsets.([comp_data{c} '_offset'])(python_index)),3) ...
                ' \pm ' num2str(nanstd(gdap_offsets.([comp_data{c} '_offset'])(python_index)),3)];
            
            legend([h1 h2], matlab_leg, python_leg, 'fontsize', 7)

        end
        
        
        dm1 = subplot(6,2,9); hold on
        plot(offsets.(list_matlab_wmos{matlab_index}).gdap.LONGITUDE_float, ...
                offsets.(list_matlab_wmos{matlab_index}).gdap.LATITUDE_float, '+', 'color', c_map(3,:))
          
        plot(offsets.(list_matlab_wmos{matlab_index}).gdap.LONGITUDE_gdap, ...
            offsets.(list_matlab_wmos{matlab_index}).gdap.LATITUDE_gdap, 'd', 'color', c_map(5,:))
        
        title(['M flt/gdp lat, '  num2str(nanmean(offsets.(list_matlab_wmos{matlab_index}).gdap.LATITUDE_float),4) ...
            '/ '  num2str(nanmean(offsets.(list_matlab_wmos{matlab_index}).gdap.LATITUDE_gdap),4) ...
            ' lon, '  num2str(nanmean(offsets.(list_matlab_wmos{matlab_index}).gdap.LONGITUDE_float),4) ...
            '/ '  num2str(nanmean(offsets.(list_matlab_wmos{matlab_index}).gdap.LONGITUDE_gdap),4)])
        
        legend('float', 'glodap', 'location', 'southoutside')
                    
        dm2 = subplot(6,2,10); hold on
        if nanmean(gdap_offsets.main_float_longitude(python_index)) <0
            float_long_python = gdap_offsets.main_float_longitude(python_index) + 360;
        else
            float_long_python = gdap_offsets.main_float_longitude(python_index);
        end
        
        plot(float_long_python, ...
            gdap_offsets.main_float_latitude(python_index), '+', 'color', c_map(3,:))
        
        
        plot(gdap_offsets.glodap_longitude(python_index), ...
                gdap_offsets.glodap_latitude(python_index), 'd', 'color', c_map(5,:))
            legend('float', 'glodap', 'location', 'southoutside')
            
        title(['Py flt/gdp lat, '  num2str(nanmean(gdap_offsets.main_float_latitude(python_index)),4) ...
            '/ '  num2str(nanmean(gdap_offsets.glodap_latitude(python_index)),4) ...
            ' lon, '  num2str(nanmean(float_long_python),4) ...
            '/ '  num2str(nanmean(gdap_offsets.glodap_longitude(python_index)),4)])
        
        d3 = subplot(6,2,11:12);
        
        python_glodap_Matlab_time = double(gdap_offsets.datetime(python_index) ...
            *1/60*1/24+datenum(1987,02,20));
        python_float_Matlab_time = double(gdap_offsets.main_float_juld(python_index) ...
            *1/(10^9)*1/60*1/60*1/24+datenum(2007,03,13,11,27,48));
        
        plot( python_glodap_Matlab_time, gdap_offsets.([comp_data{3} '_offset'])(python_index), 'o', 'color', c_map(4,:))
        hold on
        
        plot(offsets.(list_matlab_wmos{matlab_index}).gdap.GMT_Matlab_gdap, ...
            offsets.(list_matlab_wmos{matlab_index}).gdap.([comp_data{3} '_offset']), 'x', 'color', c_map(2,:))
        title('Glodap dates')
        datetick('x', 'YY-mmm')
        pos_map_1 = get(dm1, 'position');
        pos_map_2 = get(dm2, 'position');

        pos_map_2(3) = pos_map_1(3);
        pos_map_2(4) = pos_map_1(4);
%
        set(dm1, 'Position', pos_map_1+[0 -.03 0 .06])
        %
        set(dm2, 'Position', pos_map_2+[0 -.06 0 .06])
        
        pos_time = get(d3, 'position');
        set(d3, 'Position', pos_time+[0 -.05 0 0])

        plot_filename = ['Matlab_Python_comp_' num2str(wmo_search)];
        
        print(gcf, '-dpng', '-r400', [Plot_dir  plot_filename '_v3.png' ])

    end
end
