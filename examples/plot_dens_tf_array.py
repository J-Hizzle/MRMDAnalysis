# %%
import pathlib as plb
from mrmdanalysis.plot import plot_profile_array_in_adress_box, select_color_from_colormap
from mrmdanalysis.read import read_data_on_grid
from mrmdanalysis.util import calc_coupling_potential_from_thermodynamic_force, calc_applied_thermodynamic_force, calc_thermodynamic_force_from_coupling_potential, calc_thermodynamic_force_on_grid, cut_data_to_new_box
import numpy as np
# %%
file_basess = [["thermodynamicForce_initialGuess_rho296_T15_x10_yz20_2026_06_10",
                "thermodynamicForce_initialGuess_rho296_T15_x10_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho296_T15_x30_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho296_T15_x30_yz30_2026_06_10"],
               ["thermodynamicForce_initialGuess_rho370_T15_x10_yz20_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T15_x10_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T15_x30_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T15_x30_yz30_2026_06_10"],
               ["thermodynamicForce_initialGuess_rho370_T20_x10_yz20_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T20_x10_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T20_x30_yz30_2026_06_10",
                "thermodynamicForce_initialGuess_rho370_T20_x30_yz30_2026_06_10"]]

paths_data = [plb.Path("//media/julianhille/T7 Shield/noATsimulations_2026_06_10/rho296_T15_2026_06_10").resolve(),
              plb.Path("//media/julianhille/T7 Shield/noATsimulations_2026_06_10/rho370_T15_2026_06_10").resolve(),
              plb.Path("//media/julianhille/T7 Shield/noATsimulations_2026_06_10/rho370_T20_2026_06_10").resolve()]

iterations = [7, 7, 7, 0]
cuts_to_new_box = [False, False, True, True]

r_ref_base = 5.0
r_min_base = 0.0
r_max_base = 2.5
app_min_base = 0.0
app_max_base = 4.5

r_refs = [5.0, 5.0, 15.0, 15.0]
r_mins = [0.0, 0.0, 10.0, 10.0]
r_maxs = [2.5, 2.5, 12.5, 12.5]
app_mins = [0.0, 0.0, 10.0, 10.0]
app_maxs = [4.5, 4.5, 14.5, 14.5]

file_plt = None
format_plt = None

#out_base = "dens_tf_array_2026_06_22.svg"
#path_out = plb.Path('/srv/public/julianhille/publications/paper_noAT/').resolve()
#file_plt = path_out / '{0}'.format(out_base)
#format_plt = 'SVG'
# %%
dens_gridss = []
dens_profiless = []
tf_gridss = []
tf_profiless = []
pot_gridss = []
pot_profiless = []

for i, file_bases in enumerate(file_basess):
    dens_grids = []
    dens_profiles = []
    tf_grids = []
    tf_profiles = []
    pot_grids = []
    pot_profiles = []

    for j, file_base in enumerate(file_bases):
        iteration = iterations[j]
        cut_to_new_box = cuts_to_new_box[j]

        # get path to data
        file_dens = paths_data[i] / '{0}/{0}_dens.txt'.format(file_base)
        file_tf = paths_data[i] / '{0}/{0}_tf.txt'.format(file_base)

        # load data from file
        dens_grid, dens_profiles_temp = read_data_on_grid(file_dens)
        tf_grid, tf_profiles_temp = read_data_on_grid(file_tf)

        dens_profile = dens_profiles_temp[iteration]
        tf_profile = calc_applied_thermodynamic_force(tf_grid, tf_profiles_temp[iteration], r_refs[j], app_mins[j], app_maxs[j])

        if cut_to_new_box:
            dens_grid, dens_profile = cut_data_to_new_box(dens_grid, dens_profile, r_refs[j], app_mins[j], app_maxs[j] + 0.5, r_ref_base, app_min_base, app_max_base + 0.5)
            tf_grid, tf_profile = cut_data_to_new_box(tf_grid, tf_profile, r_refs[j], app_mins[j], app_maxs[j], r_ref_base, app_min_base, app_max_base)

        pot_grid = np.linspace(0, r_ref_base * 2, len(tf_grid) * 100)
        pot_grid += (pot_grid[1] - pot_grid[0]) / 2
        tf_interp = calc_thermodynamic_force_on_grid(pot_grid, tf_grid, tf_profile)
        pot_profile = calc_coupling_potential_from_thermodynamic_force(pot_grid, tf_interp, r_ref_base, r_min_base)

        if j != len(file_bases) - 1:
            dens_grids.append(dens_grid)
            dens_profiles.append(dens_profile)
        tf_grids.append(tf_grid)
        tf_profiles.append(tf_profile)
        pot_grids.append(pot_grid)
        pot_profiles.append(pot_profile)
    
    dens_gridss.append(dens_grids)
    dens_profiless.append(dens_profiles)
    tf_gridss.append(tf_grids)
    tf_profiless.append(tf_profiles)
    pot_gridss.append(pot_grids)
    pot_profiless.append(pot_profiles)

pos_gridsss = [dens_gridss, tf_gridss, pot_gridss]
data_profilesss = [dens_profiless, tf_profiless, pot_profiless]
# %%
legend_loc = (0.74, 0.75)

x_labelss = []
y_labelss = []
y_label_types = [r'$\rho$ in $\sigma^{-3}$', r'$F_{th}$ in $\epsilon/\sigma$', r'$\phi_{th}$ in $\epsilon$']

for i, pos_gridss in enumerate(pos_gridsss):
    x_labels = []
    y_labels = []
    for j in range(len(pos_gridss)):
        if i == len(pos_gridsss) - 1:
            x_labels.append(r'$x - x_{\mathrm{AT}/\Delta}$ in $\sigma$')
        else:
            x_labels.append(None)
        if j == 0:
            y_labels.append(y_label_types[i])
        else:
            y_labels.append(None)
    x_labelss.append(x_labels)
    y_labelss.append(y_labels)

x_lims_types = [[0.0, 5.0], [0.0, 5.0], [0.0, 5.0]]
x_limsss = [[x_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

y_lims_types = [[0.17, 0.4], [-3.6, 1.6], [-1.2, 0.2]]
y_limsss = [[y_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

data_labels = ["AdResS (minimal)",\
               "AdResS (no AT)",\
               "AdResS (complete)",\
               "initial guess"]

colors = [select_color_from_colormap(range_indicator/(len(dens_grids)) / 2, 'plasma') for range_indicator in range(len(dens_grids))]
colors.append('red')
linestyles = ['--', '-.', ':', '-']

scale = 5
# %%
plot_profile_array_in_adress_box(pos_gridsss, data_profilesss, x_limsss, y_limsss, x_labelss, y_labelss, data_labels, colors, r_ref_base, r_min_base, r_max_base, linestyles=linestyles, file_plt=file_plt, format_plt=format_plt, legend_loc=legend_loc, scale=scale)
# %%
