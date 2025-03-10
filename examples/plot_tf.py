# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt

from mrmdanalysis.plot import plot_profile_in_adress_box
from mrmdanalysis.read import read_data_on_grid
import os

from mrmdanalysis.util import calc_coupling_potential_from_thermodynamic_force
# %%
path_data = plb.Path('/home/mi/julianhille//mounted_directories/curta/project/mrmd/build_foss_2023a_CUDA_2025_02_24/examples/TracerThermodynamicForce').resolve()
#path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce').resolve()

dens_bin = 0.15

file_base = 'thermoForce_2025_02_26_{0}'.format(str(dens_bin).replace(".", ""))

file_dens = path_data / '{0}_dens.txt'.format(file_base)
file_tf = path_data / '{0}_tf.txt'.format(file_base)

file_plt_dens = None
file_plt_tf = None
file_plt_pot = None

path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_02_28').resolve()
file_plt_dens = path_out / '{0}_itAll_dens.png'.format(file_base)
format_plt_dens = 'PNG'

file_plt_tf = path_out / '{0}_tf.png'.format(file_base)
format_plt_tf = 'PNG'

file_plt_pot = path_out / '{0}_pot.png'.format(file_base)
format_plt_pot = 'PNG'
# %%
pos_x_grid_dens, dens_vals = read_data_on_grid(file_dens)
pos_x_grid_tf, tf_vals = read_data_on_grid(file_tf)
# %%
r_ref = 22.5
r_min = 3.0
r_max = 5.5

first_step = 0
last_step = np.shape(dens_vals)[0]
last_step = 41
step_size = 10
iterations_to_plot = [i for i in range(first_step, last_step, step_size)]
x_label = r'$x$ in $\sigma$'
colors = [plt.cm.plasma_r((i - first_step)/(last_step - first_step)) for i in iterations_to_plot]

dens_vals_to_plot = dens_vals[iterations_to_plot]
y_lims_dens = [np.min(dens_vals_to_plot) - 0.05, np.max(dens_vals_to_plot) + 0.05]
y_label_dens = r'$\rho$ in $\sigma^{-3}$'
data_labels_dens = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_dens = r'dens'

tf_vals_to_plot = tf_vals[iterations_to_plot]
y_lims_tf = [np.min(tf_vals_to_plot) - 2, np.max(tf_vals_to_plot) + 2]
y_label_tf = r'$F_{th}$ in $\epsilon/\sigma$'
data_labels_tf = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_tf = r'thermodynamic force'

pot_vals_to_plot = [calc_coupling_potential_from_thermodynamic_force(pos_x_grid_tf, tf_vals[i], r_ref, r_min) for i in iterations_to_plot]
y_lims_pot = [np.min(pot_vals_to_plot) - 0.2, np.max(pot_vals_to_plot) + 0.2]
y_label_pot = r'$\phi_{th}$ in $\epsilon$'
data_labels_pot = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_pot = r'coupling potential'
# %%
if file_plt_dens:
    if os.path.isfile(file_plt_dens):
            print('dens file already exists! Aborting.')
    else:
        plot_profile_in_adress_box(pos_x_grid_dens, dens_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens, file_plt_dens, format_plt_dens)
else:
    plot_profile_in_adress_box(pos_x_grid_dens, dens_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens)
# %%
if file_plt_tf: 
    if os.path.isfile(file_plt_tf):
            print('tf file already exists! Aborting.')
    else:
        plot_profile_in_adress_box(pos_x_grid_tf, tf_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, y_lims_tf, file_plt_tf, format_plt_tf)
else:
    plot_profile_in_adress_box(pos_x_grid_tf, tf_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, y_lims_tf)
# %%
if file_plt_pot: 
    if os.path.isfile(file_plt_pot):
            print('pot file already exists! Aborting.')
    else:
        plot_profile_in_adress_box(pos_x_grid_tf, pot_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_pot, data_labels_pot, fig_title_pot, colors, y_lims_pot, file_plt_pot, format_plt_pot)
else:
    plot_profile_in_adress_box(pos_x_grid_tf, pot_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_pot, data_labels_pot, fig_title_pot, colors, y_lims_pot)
# %%
