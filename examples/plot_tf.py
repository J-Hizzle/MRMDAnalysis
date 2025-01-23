# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt

from mrmdanalysis.plot import plot_profile_in_adress_box
from mrmdanalysis.read import read_data_on_grid
# %%
path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce').resolve()
file_dens = path_data / 'densityProfile.txt'
file_tf = path_data / 'thermodynamicForce.txt'
# %%
pos_x_grid_dens, dens_vals = read_data_on_grid(file_dens)
# %%
step_size = 20
iterations_to_plot = [i for i in range(0, np.shape(dens_vals)[0], step_size)]
r_ref = 0.0
r_min = 5.0
r_max = 7.5
x_label = r'$x - x_{ref}$ in $nm$'
colors = [plt.cm.plasma_r(i/np.max(iterations_to_plot)) for i in iterations_to_plot]
# %%
file_plt_dens = path_data / 'mrmd_dens_2025_01_23.png'
format_plt_dens = 'PNG'

dens_vals_to_plot = dens_vals[iterations_to_plot]
y_lims_dens = [np.min(dens_vals_to_plot) - 20, np.max(dens_vals_to_plot) + 20]
y_label_dens = r'$\rho$ in [Density]'
data_labels_dens = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_dens = r'dens'
# %%
plot_profile_in_adress_box(pos_x_grid_dens, dens_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens, file_plt_dens, format_plt_dens)
# %%
pos_x_grid_tf, tf_vals = read_data_on_grid(file_tf)
# %%
file_plt_tf = path_data / 'mrmd_tf_2025_01_23.png'
format_plt_tf = 'PNG'

tf_vals_to_plot = tf_vals[iterations_to_plot]
y_lims_tf = [np.min(tf_vals_to_plot) - 2, np.max(tf_vals_to_plot) + 2]
y_label_tf = r'$F_{th}$ in [Force]'
data_labels_tf = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_tf = r'thermodynamic force'
# %%
plot_profile_in_adress_box(pos_x_grid_tf, tf_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, y_lims_tf, file_plt_tf, format_plt_tf)
# %%
