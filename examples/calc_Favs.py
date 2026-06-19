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

from mrmdanalysis.interactions import calc_average_force_through_interface_for_particles_around_x_planes
# %%
file_basess = [["tracerProduction_rho296_T15_x10_yz20_2026_06_17",
                "tracerProduction_rho296_T15_x10_yz30_2026_06_17",
                "tracerProduction_rho296_T15_x30_yz30_2026_06_17",
                "atomisticProduction_rho296_T15_x30_yz30_2026_06_17"],
               ["tracerProduction_rho370_T15_x10_yz20_2026_06_17",
                "tracerProduction_rho370_T15_x10_yz30_2026_06_17",
                "tracerProduction_rho370_T15_x30_yz30_2026_06_17",
                "atomisticProduction_rho370_T15_x30_yz30_2026_06_17"],
               ["tracerProduction_rho370_T20_x10_yz20_2026_06_17",
                "tracerProduction_rho370_T20_x10_yz30_2026_06_17",
                "tracerProduction_rho370_T20_x30_yz30_2026_06_17",
                "atomisticProduction_rho370_T20_x30_yz30_2026_06_17"]]

subdomain_dimensionss = [(10.0, 10.0, 10.0),\
                        (10.0, 20.0, 20.0),\
                        (10.0, 30.0, 30.0),\
                        (30.0, 30.0, 30.0),
                        (30.0, 30.0, 30.0)]

paths_parent = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho296_T15_2026_06_10').resolve(),\
                plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho370_T15_2026_06_10').resolve(),\
                plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho370_T20_2026_06_10').resolve()]

format_trj = 'H5MD'
format_top = 'GRO'

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

file_out = None
# %%
for i, file_bases in enumerate(file_basess):
    path_parent = paths_parent[i]

    x_inter = x_inters[i]

    x_plane_grid = np.arange(x_inter - r_cut, x_inter, densBinWidth) + densShift

    for file_base in file_bases:
        file_trj = path_parent / file_base / '{0}.h5md'.format(file_base)
        file_top = path_parent / file_base / '{0}_final.gro'.format(file_base)
        file_out = path_parent / file_base / '{0}_Fav_data.txt'.format(file_base)
        average_force_vals = calc_average_force_through_interface_for_particles_around_x_planes(x_plane_grid, x_margin, x_inter, r_cut, epsilon, sigma, str(file_top), str(file_trj), format_top, format_trj, atom_name, frame_lims, frame_freq, file_out)
# %%
