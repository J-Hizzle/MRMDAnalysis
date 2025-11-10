# %%
import pathlib as plb
import MDAnalysis as mda

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.plot import select_color_from_colormap
from mrmdanalysis.statistics import count_velocities
# %%
file_basess = [["tracerProduction_23945721_2025_09_09",\
                "tracerProduction_23945722_2025_09_09",\
                "tracerProduction_23945723_2025_09_09",\
                "atomisticProduction_23945366_2025_09_09"],\
               ["tracerProduction_23945717_2025_09_09",\
                "tracerProduction_23945715_2025_09_09",\
                "tracerProduction_23945712_2025_09_09",\
                "atomisticProduction_23945334_2025_09_09"],\
               ["tracerProduction_n1480_T20_24019714_2025_09_23",\
                "tracerProduction_n3330_T20_24019715_2025_09_23",\
                "tracerProduction_n9990_T20_24019724_2025_09_23",\
                "atomisticProduction_24019672_2025_09_23"]]

paths_parent = [plb.Path('/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho0296_2025_09_03').resolve(),\
                plb.Path('/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho0370_2025_09_03').resolve(),\
                plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho0370_T20_2025_09_19').resolve()]

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
        file_top = path_parent / '{0}/{0}.gro'.format(file_base)
        file_trj = path_parent / '{0}/{0}.h5md'.format(file_base)
    
        atom_selection = 'name Ar and prop x > {0} and prop x < {1}'.format(r_refs[0] + r_mins[0] - x_margin, r_refs[0] + r_mins[0] + x_margin)
        velocity_distribution = count_velocities(velocity_bin_grid, file_top, file_trj, atom_selection, format_top, format_trj, frame_freq, frame_lims)

        write_array = np.column_stack((velocity_bin_grid, velocity_distribution)).T
        savename = path_out / 'velocity_distribution_{0}_{1}_2025_09_30.txt'.format(i, j)

        np.savetxt(savename, write_array, fmt='%05.6e')
# %%
velocity_bin_gridss = []
velocity_distributionss = []
for i, file_bases in enumerate(file_basess):
    velocity_bin_grids = []
    velocity_distributions = []

    path_parent = paths_parent[i]

    for j, file_base in enumerate(file_bases):
        file_read = path_out / 'velocity_distribution_{0}_{1}_2025_09_30.txt'.format(i, j)
        velocity_bin_grid, velocity_distribution = np.loadtxt(file_read, dtype=float, delimiter=" ")

        velocity_bin_grids.append(velocity_bin_grid)
        velocity_distributions.append(velocity_distribution)
    
    velocity_bin_gridss.append(velocity_bin_grids)
    velocity_distributionss.append(velocity_distributions)
# %%
pos_gridsss = [velocity_bin_gridss]
data_profilesss = [velocity_distributionss]
# %%
legend_loc = (0.82, 0.55)

x_labelss = []
y_labelss = []
y_label_types = [r'$f_{\mathrm{M}}(v; x_{\mathrm{AT}/\Delta})$ in $\tau/\sigma$']

for i, pos_gridss in enumerate(pos_gridsss):
    x_labels = []
    y_labels = []
    for j in range(len(pos_gridss)):
        if i == len(pos_gridsss) - 1:
            x_labels.append(r'$v$ in $\sigma/\tau$')
        else:
            x_labels.append(None)
        if j == 0:
            y_labels.append(y_label_types[i])
        else:
            y_labels.append(None)
    x_labelss.append(x_labels)
    y_labelss.append(y_labels)

x_lims_types = [[0.0, 10.0]]
x_limsss = [[x_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

y_lims_types = [[0.0, 0.05]]
y_limsss = [[y_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

data_labels = ["AdResS (minimal)",\
               "AdResS (no AT)",\
               "AdResS (full)",\
               "Reference (full)"]

axis_titles = [r'$\rho = 0.296 \sigma^{-3}$, $T = 1.5 \epsilon/k_{\mathrm{B}}$',\
               r'$\rho = 0.370 \sigma^{-3}$, $T = 1.5 \epsilon/k_{\mathrm{B}}$',\
               r'$\rho = 0.370 \sigma^{-3}$, $T = 2.0 \epsilon/k_{\mathrm{B}}$']

colors = [select_color_from_colormap(range_indicator/(len(velocity_distributionss[0]) - 2) / 2, 'plasma') for range_indicator in range(len(velocity_distributionss[0]))]
colors[-1] = 'black'
linestyles = ['--', '-.', ':', '-', '-'] 

scale = 5
# %%
fig, ax = plt.subplots(len(pos_gridsss), len(pos_gridsss[0]), figsize=(len(pos_gridsss[0]) * scale, len(pos_gridsss) * scale), constrained_layout=True, squeeze=False, )

for row_index in range(0, len(pos_gridsss)):  
    for column_index in range(0, len(pos_gridsss[0])):
        pos_grids = pos_gridsss[row_index][column_index]
        data_profiles = data_profilesss[row_index][column_index]
        x_label = x_labelss[row_index][column_index]      
        y_label = y_labelss[row_index][column_index]
        x_lims = x_limsss[row_index][column_index]
        y_lims = y_limsss[row_index][column_index]
        if axis_titles and row_index == 0:
            axis_title = axis_titles[column_index]
        else:
            axis_title = None

        plt.sca(ax[row_index, column_index])
        
        for pos_grid, data_profile, data_label, color, linestyle in zip(pos_grids, data_profiles, data_labels, colors, linestyles):
            plt.plot(pos_grid, data_profile/np.sum(data_profile), label=data_label, color=color, linewidth=0.75, linestyle=linestyle)
        
        plt.xlim(x_lims)
        plt.ylim(y_lims)
        plt.title(axis_title)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        if row_index < len(pos_gridsss) - 1:
            plt.tick_params(axis='x', labelbottom=False)
        #if row_index == 0:
        #    plt.tick_params(axis='x', labeltop=True)
        if column_index > 0:
            plt.tick_params(axis='y', labelleft=False)
            
        #if column_index == len(pos_gridsss[0]) - 1:
        #    plt.tick_params(axis='y', labelright=True)

    handles = [mlines.Line2D([], [], color=colors[i], label=data_labels[i], linestyle=linestyles[i]) for i in range(len(data_labels))]
    fig.legend(handles=handles, loc=legend_loc, framealpha=0.0)

# %%
