# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt

from mrmdanalysis.plot import plot_profile_in_adress_box
from mrmdanalysis.read import read_data_on_grid
import os
# %%
path_data = plb.Path('/home/julianhille/project/mrmd/build_foss_2023a_CUDA_2025_02_24/examples/TracerThermodynamicForce/simulations').resolve()

dens_bin = 0.2

file_base = 'noForce_2025_03_04_{0}'.format(str(dens_bin).replace(".", ""))

file_dens = path_data / '{0}_dens.txt'.format(file_base)
file_tf = path_data / '{0}_tf.txt'.format(file_base)

file_plt_dens = None

path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_02_28').resolve()
file_plt_dens = path_out / '{0}_itAll_dens.png'.format(file_base)
format_plt_dens = 'PNG'
# %%
pos_x_grid_dens, dens_vals = read_data_on_grid(file_dens)
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
# %%
if file_plt_dens:
    if os.path.isfile(file_plt_dens):
            print('dens file already exists! Aborting.')
    else:
        plot_profile_in_adress_box(pos_x_grid_dens, dens_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens, file_plt_dens, format_plt_dens)
else:
    plot_profile_in_adress_box(pos_x_grid_dens, dens_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_dens, data_labels_dens, fig_title_dens, colors, y_lims_dens)
