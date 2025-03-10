# %%
import pathlib as plb
import numpy as np
from mrmdanalysis.plot import plot_profile_in_adress_box, plot_profile_in_half_adress_box, select_color_from_colormap
from mrmdanalysis.read import read_data_on_grid
from mrmdanalysis.util import calc_coupling_potential_from_thermodynamic_force, calc_applied_thermodynamic_force
# %%
path_data = plb.Path('/home/mi/julianhille//mounted_directories/curta/project/mrmd/build_foss_2023a_CUDA_2025_02_24/examples/TracerThermodynamicForce').resolve()
#path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce').resolve()

dens_bin_widths = [0.1, 0.125, 0.15, 0.2] 

file_bases = ['thermoForce_2025_02_26_' + str(dens_bin).replace(".", "") for dens_bin in dens_bin_widths]

files_dens = [path_data / '{0}_dens.txt'.format(file_base) for file_base in file_bases]
files_tf = [path_data / '{0}_tf.txt'.format(file_base) for file_base in file_bases]

file_plt_dens = None
format_plt_dens = None
file_plt_tf = None
format_plt_tf = None
file_plt_pot = None
format_plt_pot = None

out_base = "dens_res_comparison"

path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_02_28').resolve()

file_plt_dens = path_out / '{0}.png'.format(out_base)
format_plt_dens = 'PNG'
file_plt_tf = path_out / '{0}_tf.png'.format(out_base)
format_plt_tf = 'PNG'
file_plt_pot = path_out / '{0}_pot.png'.format(out_base)
format_plt_pot = 'PNG'

r_ref = 22.5
r_min = 3.0
r_max = 5.5
app_min = 3.0
app_max = 8.0
# %%
dens_valss = []
dens_grids = []
for file_dens in files_dens:
    dens_grid, dens_vals = read_data_on_grid(file_dens)
    dens_valss.append(dens_vals[-1])
    dens_grids.append(dens_grid)

tf_valss = []
tf_grids = []
pot_valss = []
pot_grids = []
for file_tf in files_tf:
    tf_grid, tf_vals = read_data_on_grid(file_tf)
    tf_vals = calc_applied_thermodynamic_force(tf_grid, tf_vals[-1], r_ref, app_min, app_max)
    print(len(tf_vals))
    pot_grid = tf_grid
    pot_vals = calc_coupling_potential_from_thermodynamic_force(pot_grid, tf_vals, r_ref, r_min)

    tf_valss.append(tf_vals)
    tf_grids.append(tf_grid)
    pot_grids.append(pot_grid)
    pot_valss.append(pot_vals)
# %%
colors = [select_color_from_colormap(range_indicator/len(dens_grids), 'plasma_r') for range_indicator in range(len(dens_grids))]
x_label = r'$x$ in $\sigma$'

y_lims_dens = [np.min(dens_valss[0]) - 0.05, np.max(dens_valss[0]) + 0.05]
y_label_dens = r'$\rho$ in $\sigma^{-3}$'
data_labels_dens = [r'bin width = ' + str(bin_width) for bin_width in dens_bin_widths]
fig_title_dens = r'MRMD density profiles'

y_lims_tf = [np.min(tf_valss[0]) - 1, 4.0]
x_lims_tf = [app_min - 1.0, app_max]
y_label_tf = r'$F_{th}$ in $\epsilon/\sigma$'
data_labels_tf = [r'bin width = ' + str(bin_width) for bin_width in dens_bin_widths]
fig_title_tf = r'MRMD thermodynamic forces'

y_lims_pot = [np.min(pot_valss[0]) - 0.2, np.max(pot_valss[0]) + 0.2]
x_lims_pot = [app_min - 1.0, app_max]
y_label_pot = r'$\phi_{th}$ in $\epsilon$'
data_labels_pot = [r'bin width = ' + str(bin_width) for bin_width in dens_bin_widths]
fig_title_pot = r'MRMD coupling potentials'
# %%
plot_profile_in_adress_box(dens_grids, dens_valss, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens, file_plt_dens, format_plt_dens)
# %%
plot_profile_in_half_adress_box(tf_grids, tf_valss, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, x_lims_tf, y_lims_tf, file_plt_tf, format_plt_tf)
# %%
plot_profile_in_half_adress_box(pot_grids, pot_valss, r_ref, r_min, r_max, x_label, y_label_pot, data_labels_pot, fig_title_pot, colors, x_lims_pot, y_lims_pot, file_plt_pot, format_plt_pot)

# %%
