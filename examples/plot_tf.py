# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt

from mrmdanalysis.plot import plot_profile_in_half_adress_box, plot_profile_in_adress_box
from mrmdanalysis.read import read_data_on_grid

from mrmdanalysis.util import calc_coupling_potential_from_thermodynamic_force, calc_applied_thermodynamic_force
# %%
#file_base = "test"
dir_base = "tracerProduction_23828869_2025_08_15"
file_base = "thermoForce_NoAT_n15000_23827578_2025_08_15_final"

#path_data = plb.Path("/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce").resolve()
path_data = plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/{}'.format(dir_base)).resolve()
#path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/DeltaDeltaInterface').resolve()

#dens_bin = 0.15
#file_base = 'thermoForce_2025_02_26_{0}'.format(str(dens_bin).replace(".", ""))


file_tf = path_data / '{0}_tf.txt'.format(file_base)

file_plt_tf = None
format_plt_tf = None
file_plt_pot = None
format_plt_pot = None

#path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_06_10').resolve()
#file_plt_tf = path_out / '{0}_tf.png'.format(file_base)
#format_plt_tf = 'PNG'
#
#file_plt_pot = path_out / '{0}_pot.png'.format(file_base)
#format_plt_pot = 'PNG'
# %%
pos_x_grid_tf, tf_vals = read_data_on_grid(file_tf)
# %%
tf_grids = [pos_x_grid_tf for i in range(len(tf_vals))]
pot_grids = tf_grids
binWidth = pos_x_grid_tf[1] - pos_x_grid_tf[0]

r_ref = 22.5
r_min = 3.0
r_max = 5.5
app_min = 3.0
app_max = 8.0

first_step = 0
last_step = np.shape(tf_vals)[0]
#last_step = 11
step_size = 1
iterations_to_plot = [i for i in range(first_step, last_step, step_size)]

x_label = r'$x$ in $\sigma$'
colors = ["blue"]
loc='lower right'

tf_vals_to_plot = [calc_applied_thermodynamic_force(tf_grids[i], tf_vals[i], r_ref, app_min, app_max) for i in iterations_to_plot]
#tf_vals_to_plot = [tf_vals[i] for i in iterations_to_plot]
y_lims_tf = [-10.0, 6.0]
x_lims_tf = [r_min - 0.5, r_max + 4.5]
y_label_tf = r'$F_{th}$ in $\epsilon/\sigma$'
data_labels_tf = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_tf = r'thermodynamic force ($\Delta x = {0:.2f} \sigma$)'.format(binWidth)

pot_vals_to_plot = [calc_coupling_potential_from_thermodynamic_force(pos_x_grid_tf, tf_vals[i], r_ref, r_min) for i in iterations_to_plot]
y_lims_pot = [-1.3, 0.1]
x_lims_pot = x_lims_tf
y_label_pot = r'$\phi_{th}$ in $\epsilon$'
data_labels_pot = [r'iteration = ' + str(i) for i in iterations_to_plot]
fig_title_pot = r'coupling potential ($\Delta x = {0:.2f} \sigma$)'.format(binWidth)
# %%
plot_profile_in_half_adress_box(tf_grids, tf_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, x_lims_tf, y_lims_tf, file_plt_tf, format_plt_tf, loc=loc)
# %%
plot_profile_in_adress_box(tf_grids, tf_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_tf, data_labels_tf, fig_title_tf, colors, y_lims_tf, file_plt_tf, format_plt_tf)
# %%
plot_profile_in_half_adress_box(pot_grids, pot_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_pot, data_labels_pot, fig_title_pot, colors, x_lims_pot, y_lims_pot, file_plt_pot, format_plt_pot, loc=loc)
# %%
plot_profile_in_adress_box(pot_grids, pot_vals_to_plot, r_ref, r_min, r_max, x_label, y_label_pot, data_labels_pot, fig_title_pot, colors, y_lims_pot, file_plt_pot, format_plt_pot)
# %%
