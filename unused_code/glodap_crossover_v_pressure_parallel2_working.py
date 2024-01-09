from multiprocessing import Process
import pressure_level_glodap_mean as pl
import xarray as xr
import os
from datetime import datetime

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

todays_date = datetime.today().strftime('%Y_%m_%d')

grouped_plot_dir = offset_dir + 'grouped_plots/' + todays_date + '/'
# grouped_plot_dir = offset_dir + 'grouped_plots/levels_mean_saved/'

if not os.path.isdir(grouped_plot_dir):
    os.mkdir(grouped_plot_dir)


def main():
    # loop through various glodap_offset plots, calculating trimmed mean for different pressure levels and storing for comparison
    glodap_offsets = []
    # glodap_offsets_filenames = ['glodap_offsets_100km_1_to_550_50m_0.05dens_0.05spice_6.nc', 'glodap_offsets_100km_400_to_2100_100m_0.005dens_0.005spice_6.nc']
    glodap_offsets_filenames = ['glodap_offsets_100km_1400_to_2100_100m_0.005dens_0.005_spice_7.nc']
    for filename in glodap_offsets_filenames:
        ds = xr.load_dataset(output_dir+filename)
        glodap_offsets.append(ds)
    pressure_levels = [1500, 2000]  # Adjust as needed

    # pressure_levels = [0, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 2000]  # Adjust as needed
    # pressure_levels = [0, 100, 200]#, 300]#, 400, 500, 750, 1000, 1250, 1500, 2000]  # Adjust as needed

    year_filt = 0
    year_plus_minus = 5
    var_list_plot = ['PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED','DOXY_ADJUSTED','NITRATE_ADJUSTED',
                        'DIC','pH_25C_TOTAL_ADJUSTED','PH_IN_SITU_TOTAL_ADJUSTED','PDENS']

    process_list = []
    for idx, gdap_offsets_n in enumerate(glodap_offsets):

        # trimmed_means = {}
        # variables with offsets that you want to trim


        for j in range(len(pressure_levels) - 1):
            # make a temporary copy of gdap_offsets_n to use for the pressure level
            gdap_offsets_n_temp = gdap_offsets_n.copy()       
            if __name__ == '__main__':
                p = Process(target=pl.pressure_level_filter, args=(argo_path,grouped_plot_dir,glodap_offsets_filenames[idx][0:-3],gdap_offsets_n_temp, \
                                                                var_list_plot,year_filt,pressure_levels[j],pressure_levels[j+1],year_plus_minus, ))
                                                            #    output_dir=grouped_plot_dir, out_filename = glodap_offsets_filenames[idx][0:-3], \
                                    #  gdap_offsets_n_temp=gdap_offsets_n_temp, var_list_plot=var_list_plot, year_filt=year_filt, pressure_level_min=pressure_levels[j],\
                                    #   pressure_level_max=pressure_levels[j+1], year_plus_minus=year_plus_minus,))
                process_list.append(p)
        
    for process in process_list:
        process.start()

    for process in process_list:
        process.join()

if __name__ == '__main__':
    main()
        