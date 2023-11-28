import argo_interp_and_crossover as aiac


def process_file(file, argo_path_interpolated, offset_dir, dist, delta_dens, delta_spice, delta_press, gdap_p, p_interp, plot_profile, var_list_plot):
    result = aiac.glodap_crossover_offsets(argo_path_interpolated, offset_dir, file, dist, delta_dens, delta_spice, delta_press,
                                           gdap_p, p_interp, plot_profile, var_list_plot)
    return file, result