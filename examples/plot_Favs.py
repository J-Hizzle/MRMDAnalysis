# %%
'''
Calculate average force through interface numerically inspired by Liouville-type hierarchy paper.
'''
import numpy as np
from scipy.integrate import cumulative_trapezoid

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

import pathlib as plb

from mrmdanalysis.interactions import calc_average_force_through_interface_for_particles_around_x_planes
# %%
path_data = plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mosceto/testing/ghost_face_lj_fa_prod_26_03_2024').resolve()
file_top = path_data / 'topol.tpr'
file_trj = path_data / 'traj_comp.xtc'
format_top = 'TPR'
format_trj = 'XTC'

#file_out = None
file_out = path_data / 'interface_force_test_2025_08_21.txt'

r_min = 3.0
r_max = 5.5
r_cut = 2.5
epsilon = 1
sigma = 1
r_ref = 22.5

densBinWidth = 0.2
densShift = densBinWidth/2.0

forceBinWidth = 0.01
forceShift = forceBinWidth/2.0

x_plane_grid = np.arange(r_ref + r_min, r_ref + r_max, densBinWidth) + densShift
x_margin = 0.05

atom_name = 'Ar'
frame_lims = [0, -1]
frame_freq = 100

file_plt = path_data / 'interaction_correction_numerical_test_2025_08_21.png'
# %%
average_force_vals = calc_average_force_through_interface_for_particles_around_x_planes(x_plane_grid, x_margin, r_ref, r_max, r_cut, epsilon, sigma, str(file_top), str(file_trj), format_top, format_trj, atom_name, frame_lims, frame_freq, file_out)
# %%
_, average_force_vals_x, average_force_vals_y, average_force_vals_z = np.loadtxt(file_out, dtype=float).T
# %%
average_pot_vals_x = -cumulative_trapezoid(average_force_vals_x, x_plane_grid, initial=0)
# %%
ylims = [np.min(average_force_vals_x) - (np.max(average_force_vals_x) - np.min(average_force_vals_x))/10, np.max(average_force_vals_x) + (np.max(average_force_vals_x) - np.min(average_force_vals_x))/10]
xlims = [2.5, 8.0]
alpha_val = 0.2
textheight = ylims[1] - (ylims[1] - ylims[0])/20

plt.plot((x_plane_grid - r_ref), average_force_vals_x, label='average force values x', color='green', linewidth=0.75)


plt.axvline(-r_min, color='black', alpha=2*alpha_val)
plt.axvline(r_min, color='black', alpha=2*alpha_val)
plt.axvline(-r_max, color='black', alpha=2*alpha_val)
plt.axvline(r_max, color='black', alpha=2*alpha_val)

plt.text(xlims[0] + (r_min - xlims[0])/2, textheight, r'AT', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.text(r_min + (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
#plt.text(-r_min - (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.text(xlims[1] - (xlims[1] - r_max)/2, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
#plt.text(-r_max - r_min, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.ylim(ylims[0], ylims[1])
plt.xlim(xlims[0], xlims[1])
#plt.yscale('log')

plt.xlabel(r'$x$ in $\sigma$')
plt.ylabel(r'$V$ in $\epsilon$')
plt.title(r'Force profile through $\Delta$/TR interfaces ($x_{margin} = ' + str(x_margin) + r'\sigma $)')
plt.legend(framealpha=0.0, fontsize=12, loc='lower right')
#plt.savefig(fname=file_plt, dpi=150)
# %%
ylims = [-1, 0.25]
xlims = [2.5, 8.0]
alpha_val = 0.2
textheight = ylims[1] - (ylims[1] - ylims[0])/20

plt.plot((x_plane_grid - r_ref), average_pot_vals_x, label='average pot values x', color='green', linewidth=0.75)


plt.axvline(-r_min, color='black', alpha=2*alpha_val)
plt.axvline(r_min, color='black', alpha=2*alpha_val)
plt.axvline(-r_max, color='black', alpha=2*alpha_val)
plt.axvline(r_max, color='black', alpha=2*alpha_val)

plt.text(xlims[0] + (r_min - xlims[0])/2, textheight, r'AT', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.text(r_min + (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
#plt.text(-r_min - (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.text(xlims[1] - (xlims[1] - r_max)/2, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
#plt.text(-r_max - r_min, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
plt.ylim(ylims[0], ylims[1])
plt.xlim(xlims[0], xlims[1])
#plt.yscale('log')

plt.xlabel(r'$x$ in $\sigma$')
plt.ylabel(r'$V$ in $\epsilon$')
plt.title(r'Potential profile through $\Delta$/TR interfaces ($x_{margin} = ' + str(x_margin) + r'\sigma $)')
plt.legend(framealpha=0.0, fontsize=12, loc='lower right')
#plt.savefig(fname=file_plt, dpi=150)
# %%
