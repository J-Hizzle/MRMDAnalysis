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

from mrmdanalysis.interactions import calc_average_force_through_interface_for_particles_around_x_planes, calc_average_force_vector_per_particle_in_sample_group
from mrmdanalysis.plot import select_color_from_colormap, plot_profile_at_interface
# %%
#file_bases = ["tracerProduction_23828863_2025_08_15"]

file_bases = ["tracerProduction_23828863_2025_08_15",\
              "tracerProduction_23843205_2025_08_19",\
              "tracerProduction_23848988_2025_08_21",\
              "tracerProduction_23828869_2025_08_15",\
              "atomisticProduction_23854725_2025_08_23",\
              "atomisticProduction_23849820_2025_08_21"]
files_trj = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/{0}/{0}.h5md'.format(file_base)).resolve() for file_base in file_bases]
files_top = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/{0}/{0}.gro'.format(file_base)).resolve() for file_base in file_bases]
format_trj = 'H5MD'
format_top = 'GRO'

fig_title_Fav = r"Average force through AT/$\Delta$-interface"
xlabel = r'$x - x_{AT/\Delta}$ in $\sigma$'
ylabel = r'$F_{av}$ in $\epsilon/\sigma$'

labels = [r"AdResS ($N = 370$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 3330$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 9990$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 15000$, $\gamma = 20$, $d_{TR} = 34$)", r"Full-AT ($N = 9990$, $\gamma = 20$)", r"Full-AT ($N = 9990$, $\gamma = 1$)"]
markers = ["^", "o", "v", "s", "D", "P"]
colors = [select_color_from_colormap(range_indicator/len(file_bases), 'plasma') for range_indicator in range(len(file_bases))]
fontsize = 8

file_plt_Fav = None

out_base = "average_interface_forces_2025_08_26"
path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_08_26').resolve()
files_out = [path_out / '{0}_data.txt'.format(file_base) for file_base in file_bases]
file_plt_Fav = path_out / '{0}.png'.format(out_base)

x_inters =  [5.0, 5.0, 25.0,  25.5, 25.0, 25.0]

cuts_to_new_box = [False, False, True, True, True, True]

r_cut = 2.5
epsilon = 1.0
sigma = 1.0

densBinWidth = 0.2
densShift = densBinWidth/2.0

#forceBinWidth = 0.01
#forceShift = forceBinWidth/2.0

x_margin = 0.1

atom_name = 'Ar'
frame_lims = [0, -1]
frame_freq = 1

ylims = [-1.5, 1.5]
xlims = [-2.5, 0.2]
fontsize = 8
# %%
for i, file_trj in enumerate(files_trj):
    file_top = files_top[i]
    file_out = files_out[i]

    x_inter = x_inters[i]

    x_plane_grid = np.arange(x_inter - r_cut, x_inter, densBinWidth) + densShift
    average_force_vals = calc_average_force_through_interface_for_particles_around_x_planes(x_plane_grid, x_margin, x_inter, r_cut, epsilon, sigma, str(file_top), str(file_trj), format_top, format_trj, atom_name, frame_lims, frame_freq, file_out)
# %%
x_plane_grids = []
average_force_valss_x = []

for i, file_trj in enumerate(files_trj):
    print("file number = ", i)
    x_plane_grid, average_force_vals_x, average_force_vals_y, average_force_vals_z = np.loadtxt(files_out[i], dtype=float).T

    x_plane_grids.append(x_plane_grid - x_inters[i])
    average_force_valss_x.append(average_force_vals_x)
# %%
plot_profile_at_interface(x_plane_grids, average_force_valss_x, 0.0, x_label=xlabel, y_label=ylabel, data_labels=labels, fig_title=fig_title_Fav, colors=colors, x_lims=xlims, y_lims=ylims, fontsize=fontsize, file_plt=file_plt_Fav)
# %%
