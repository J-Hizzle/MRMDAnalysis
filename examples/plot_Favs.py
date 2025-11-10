# %%
'''
Calculate average force through interface numerically inspired by Liouville-type hierarchy paper.
'''
import numpy as np

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

import pathlib as plb

from mrmdanalysis.interactions import calc_average_force_through_interface_for_particles_around_x_planes, calc_average_force_vector_per_particle_in_sample_group
from mrmdanalysis.plot import select_color_from_colormap, plot_profile_at_interface
# %%
#file_bases = ["tracerProduction_23945719_2025_09_09",\
#              "tracerProduction_23945717_2025_09_09",\
#              "tracerProduction_23945715_2025_09_09",\
#              "tracerProduction_23945712_2025_09_09",\
#              "atomisticProduction_23945334_2025_09_09"]
#file_bases = ["tracerProduction_23945720_2025_09_09",\
#              "tracerProduction_23945721_2025_09_09",\
#              "tracerProduction_23945722_2025_09_09",\
#              "tracerProduction_23945723_2025_09_09",\
#              "atomisticProduction_23945366_2025_09_09"]
#file_bases = ["tracerProduction_23949517_2025_09_10",\
#              "tracerProduction_23949521_2025_09_10",\
#              "tracerProduction_23949523_2025_09_10",\
#              "tracerProduction_23956781_2025_09_12",\
#              "atomisticProduction_23945374_2025_09_09"]
#file_bases = ["tracerProduction_23949525_2025_09_10",\
#              "tracerProduction_23949527_2025_09_10",\
#              "tracerProduction_23949529_2025_09_10",\
#              "tracerProduction_23956779_2025_09_12",\
#              "atomisticProduction_23945376_2025_09_09"]
file_bases = ["tracerProduction_n370_T20_24019679_2025_09_23",\
              "tracerProduction_n1480_T20_24019714_2025_09_23",\
              "tracerProduction_n3330_T20_24019715_2025_09_23",\
              "tracerProduction_n9990_T20_24019724_2025_09_23",\
              "atomisticProduction_24019672_2025_09_23"]

subdomain_dimensionss = [(10.0, 10.0, 10.0),\
                        (10.0, 20.0, 20.0),\
                        (10.0, 30.0, 30.0),\
                        (30.0, 30.0, 30.0),
                        (30.0, 30.0, 30.0)]

density = 0.370
#density = 0.296
#density = 0.247
#density = 0.198

#path_parent = plb.Path("/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations").resolve()
path_parent = plb.Path("/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations").resolve()

#files_trj = [path_parent / 'rho{0:04d}_2025_09_03/{1}/{1}.h5md'.format(int(density * 1000), file_base) for file_base in file_bases]
#files_top = [path_parent / 'rho{0:04d}_2025_09_03/{1}/{1}.gro'.format(int(density * 1000), file_base) for file_base in file_bases]
files_trj = [path_parent / 'rho{0:04d}_T20_2025_09_19/{1}/{1}.h5md'.format(int(density * 1000), file_base) for file_base in file_bases]
files_top = [path_parent / 'rho{0:04d}_T20_2025_09_19/{1}/{1}.gro'.format(int(density * 1000), file_base) for file_base in file_bases]

format_trj = 'H5MD'
format_top = 'GRO'

fig_title_Fav = "Average force through AT/${2}$-interface (${1} = {0:1.3f}$)".format(density, r"\rho", r"\Delta")
xlabel = r'$x - x_{AT/\Delta}$ in $\sigma$'
ylabel = r'$F_{av}$ in $\epsilon/\sigma$'

labels = ["$x = {0}, y = {1}, z = {2}$".format(subdomain_dimensions[0], subdomain_dimensions[1], subdomain_dimensions[2]) for subdomain_dimensions in subdomain_dimensionss]
labels[-1] = "$x = {0}, y = {1}, z = {2}$, Full-AT".format(subdomain_dimensionss[-1][0], subdomain_dimensionss[-1][1], subdomain_dimensionss[-1][2])

colors = [select_color_from_colormap(range_indicator/len(file_bases), 'plasma') for range_indicator in range(len(file_bases))]
fontsize = 8

file_plt_Fav = None

paths_out = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho{0:04d}_T20_2025_09_19/{1}'.format(int(density * 1000), file_base)).resolve() for file_base in file_bases]
files_out = [paths_out[i] / '{0}_Fav_data.txt'.format(file_base) for i, file_base in enumerate(file_bases)]

path_plt = plb.Path('/srv/public/julianhille/project/MRMDAnalysis/data').resolve()
plt_base = "average_interface_forces_rho{0:04d}_T20_2025_09_25".format(int(density * 1000))
file_plt_Fav = path_plt / '{0}.png'.format(plt_base)

x_inters =  [5.0, 5.0, 5.0, 25.0, 25.0]

r_cut = 2.5
epsilon = 1.0
sigma = 1.0

densBinWidth = 0.2
densShift = densBinWidth/2.0

#forceBinWidth = 0.01
#forceShift = forceBinWidth/2.0

x_margin = densBinWidth/2.0

atom_name = 'Ar'
frame_lims = [0, -1]
frame_freq = 1

ylims = [-1.5, 1.5]
xlims = [-2.5, 0.2]
fontsize = 8
# %%
#for i, file_trj in enumerate(files_trj):
#    file_top = files_top[i]
#    file_out = files_out[i]
#
#    x_inter = x_inters[i]
#
#    x_plane_grid = np.arange(x_inter - r_cut, x_inter, densBinWidth) + densShift
#    average_force_vals = calc_average_force_through_interface_for_particles_around_x_planes(x_plane_grid, x_margin, x_inter, r_cut, epsilon, sigma, str(file_top), str(file_trj), format_top, format_trj, atom_name, frame_lims, frame_freq, file_out)
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
