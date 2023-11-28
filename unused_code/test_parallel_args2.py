import argparse
import xarray as xr
import pressure_level_glodap_mean as pl
from multiprocessing import Process
import concurrent.futures

def main(argo_path, output_dir, grouped_plot_dir, glodap_offsets_filenames, pressure_levels, year_filt, year_plus_minus):
    # Your existing script logic goes here
    print("Executing main_pressure_process with the following arguments:")
    print(f"argo_path: {argo_path}")
    print(f"output_dir: {output_dir}")
    print(f"grouped_plot_dir: {grouped_plot_dir}")
    print(f"glodap_offsets_filenames: {glodap_offsets_filenames}")
    print(f"pressure_levels: {pressure_levels}")
    print(f"year_filt: {year_filt}")
    print(f"year_plus_minus: {year_plus_minus}")


    print('starting main pressure level processing')
    glodap_offsets = []

    for filename in glodap_offsets_filenames:
        ds = xr.load_dataset(output_dir+filename)
        glodap_offsets.append(ds)
    # pressure_levels = [0, 100, 200]#, 300]#, 400, 500, 750, 1000, 1250, 1500, 2000]  # Adjust as needed

    # year_filt = 1
    # year_plus_minus = 5
    var_list_plot = ['PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED','DOXY_ADJUSTED','NITRATE_ADJUSTED',
                        'DIC','pH_25C_TOTAL_ADJUSTED','PH_IN_SITU_TOTAL_ADJUSTED','PDENS']

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []

        for idx, gdap_offsets_n in enumerate(glodap_offsets):
            for j in range(len(pressure_levels) - 1):
                gdap_offsets_n_temp = gdap_offsets_n.copy()
                future = executor.submit(
                    pl.pressure_level_filter,
                    argo_path,
                    grouped_plot_dir,
                    glodap_offsets_filenames[idx][0:-3],
                    gdap_offsets_n_temp,
                    var_list_plot,
                    year_filt,
                    pressure_levels[j],
                    pressure_levels[j+1],
                    year_plus_minus
                )
                futures.append(future)

        # Wait for all tasks to complete
        concurrent.futures.wait(futures)
    # process_list = []
    # for idx, gdap_offsets_n in enumerate(glodap_offsets):
    #     # print(idx)
    #     # print(gdap_offsets_n)
    #     # variables with offsets that you want to trim


    #     for j in range(len(pressure_levels) - 1):
    #         # make a temporary copy of gdap_offsets_n to use for the pressure level
    #         gdap_offsets_n_temp = gdap_offsets_n.copy()       
    #         p = Process(target=pl.pressure_level_filter, args=(argo_path,grouped_plot_dir,glodap_offsets_filenames[idx][0:-3],gdap_offsets_n_temp, \
    #                                                         var_list_plot,year_filt,pressure_levels[j],pressure_levels[j+1],year_plus_minus, ))
    #         process_list.append(p)
    #         print(process_list)

    # for process in process_list:
    #     print('starting processes')
    #     process.start()

    # for process in process_list:
    #     print('joining processes')
    #     process.join()

if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Script to process pressure data")

    # Add command line arguments
    parser.add_argument("--argo_path", type=str, required=True, help="Path to Argo data")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to output directory")
    parser.add_argument("--grouped_plot_dir", type=str, required=True, help="Path to grouped plot directory")
    parser.add_argument("--glodap_offsets_filenames", nargs='+', required=True, help="List of Glodap offset filenames")
    parser.add_argument("--pressure_levels", nargs='+', type=int, required=True, help="List of pressure levels")
    parser.add_argument("--year_filt", type=int, required=True, help="Year filter")
    parser.add_argument("--year_plus_minus", type=int, required=False, help="Year plus/minus")

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(
        args.argo_path,
        args.output_dir,
        args.grouped_plot_dir,
        args.glodap_offsets_filenames,
        args.pressure_levels,
        args.year_filt,
        args.year_plus_minus
    )