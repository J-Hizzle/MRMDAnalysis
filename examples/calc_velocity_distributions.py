# %%
import pathlib as plb

import numpy as np

from mrmdanalysis.statistics import count_velocities
# %%
file_basess = [["tracerProduction_rho296_T15_x10_yz20_2026_06_11",
                "tracerProduction_rho296_T15_x10_yz30_2026_06_11",
                "tracerProduction_rho296_T15_x30_yz30_2026_06_11",
                "atomisticProduction_rho296_T15_x30_yz30_2026_06_11"],
               ["tracerProduction_rho370_T15_x10_yz20_2026_06_11",
                "tracerProduction_rho370_T15_x10_yz30_2026_06_11",
                "tracerProduction_rho370_T15_x30_yz30_2026_06_11",
                "atomisticProduction_rho370_T15_x30_yz30_2026_06_12"],
               ["tracerProduction_rho370_T20_x10_yz20_2026_06_11",
                "tracerProduction_rho370_T20_x10_yz30_2026_06_11",
                "tracerProduction_rho370_T20_x30_yz30_2026_06_11",
                "atomisticProduction_rho370_T20_x30_yz30_2026_06_12"]]

paths_parent = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho296_T15_2026_06_10').resolve(),\
                plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho370_T15_2026_06_10').resolve(),\
                plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/rho370_T20_2026_06_10').resolve()]

format_trj = 'H5MD'
format_top = 'GRO'

atom_group_1_name = 'Ar'
atom_group_2_name = 'Ar'

r_refs = [5.0, 5.0, 15.0, 15.0]
r_mins = [0.0, 0.0, 10.0, 10.0]
r_maxs = [2.5, 2.5, 12.5, 12.5]
app_mins = [0.0, 0.0, 10.0, 10.0]
app_maxs = [4.5, 4.5, 14.5, 14.5]

frame_lims = [0, -1]
frame_freq = 1

x_margin = 0.1

velocity_bin_grid = np.linspace(0, 10, 100)

path_out = plb.Path('/srv/public/julianhille/project/MRMDAnalysis/data').resolve()
# %%
'''
Calculate velocity distributions
'''
for i, file_bases in enumerate(file_basess):
    path_parent = paths_parent[i]

    for j, file_base in enumerate(file_bases):
        savename = path_parent / file_base / '{0}_velocityDistribution.txt'.format(file_base)
        if savename.is_file():
            print("File {0} already exists. Skipping calculation.".format(savename))
            continue

        file_top = path_parent / '{0}/{0}_final.gro'.format(file_base)
        file_trj = path_parent / '{0}/{0}.h5md'.format(file_base)
    
        atom_selection = 'name Ar and prop x > {0} and prop x < {1}'.format(r_refs[0] + r_mins[0] - x_margin, r_refs[0] + r_mins[0] + x_margin)
        velocity_distribution = count_velocities(velocity_bin_grid, file_top, file_trj, atom_selection, format_top, format_trj, frame_freq, frame_lims)

        write_array = np.column_stack((velocity_bin_grid, velocity_distribution)).T

        np.savetxt(savename, write_array, fmt='%05.6e')

# %%
