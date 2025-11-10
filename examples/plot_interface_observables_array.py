# %%
import pathlib as plb
import MDAnalysis as mda
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
plt.rc('text.latex', preamble=r'\usepackage{siunitx}')
import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": r"\usepackage{siunitx}",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

from mrmdanalysis.plot import select_color_from_colormap, plot_profile_at_interface, plot_PofNs
from mrmdanalysis.read import read_mrmd_out
from mrmdanalysis.statistics import count_velocities
# %%
file_basess_velocity_distribution = [["tracerProduction_23945721_2025_09_09",\
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
file_basess_flux_distribution = [["tracerProduction_23945721_2025_09_09/slurm-23945721.out",\
                "tracerProduction_23945722_2025_09_09/slurm-23945722.out",\
                "tracerProduction_23945723_2025_09_09/slurm-23945723.out",\
                "atomisticProduction_23945366_2025_09_09/slurm-23945366.out"],\
                ["tracerProduction_23945717_2025_09_09/slurm-23945717.out",\
                "tracerProduction_23945715_2025_09_09/slurm-23945715.out",\
                "tracerProduction_23945712_2025_09_09/slurm-23945712.out",\
                "atomisticProduction_23945334_2025_09_09/slurm-23945334.out"],\
                ["tracerProduction_n1480_T20_24019714_2025_09_23/slurm-24019714.out",\
                "tracerProduction_n3330_T20_24019715_2025_09_23/slurm-24019715.out",\
                "tracerProduction_n9990_T20_24019724_2025_09_23/slurm-24019724.out",\
                "atomisticProduction_24019672_2025_09_23/slurm-24019672.out"]]
file_basess_Fav = [["tracerProduction_23945721_2025_09_09",\
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
                plb.Path('/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho0370_T20_2025_09_19').resolve()]

format_trj = 'H5MD'
format_top = 'GRO'

atom_group_1_name = 'Ar'
atom_group_2_name = 'Ar'

r_ref_base = 5.0
r_min_base = 0.0
r_max_base = 2.5
app_min_base = 0.0
app_max_base = 4.5

r_refs = [5.0, 5.0, 15.0, 15.0]
r_mins = [0.0, 0.0, 10.0, 10.0]
r_maxs = [2.5, 2.5, 12.5, 12.5]
app_mins = [0.0, 0.0, 10.0, 10.0]
app_maxs = [4.5, 4.5, 14.5, 14.5]

frame_lims = [0, -1]
frame_freq = 1

x_margin = 0.1

velocity_bin_grid = np.linspace(0, 10, 100)

subdomain_dimensionss = [(10.0, 20.0, 20.0),\
                         (10.0, 30.0, 30.0),\
                         (30.0, 30.0, 30.0),
                         (30.0, 30.0, 30.0)]

interface_areas = [subdomain_dimensions[1] * subdomain_dimensions[2] for subdomain_dimensions in subdomain_dimensionss] # in sigma^2
time_interval = 10000 * 0.002 # in tau
conversion_factors = [np.sqrt(interface_area)/(interface_area * time_interval) for interface_area in interface_areas]

x_inters =  [5.0, 5.0, 25.0, 25.0]

path_read = plb.Path('/srv/public/julianhille/project/MRMDAnalysis/data').resolve()

file_plt = None
format_plt = None

out_base = "interface_observables_array_2025_11_05.PGF"
#path_out = plb.Path('/srv/public/julianhille/project/MRMDAnalysis/data').resolve()
path_out = plb.Path('/srv/public/julianhille/publications/paper_noAT/').resolve()
file_plt = path_out / '{0}'.format(out_base)
format_plt = 'PGF'
# %%
'''
Read velocity distributions
'''
velocity_bin_gridss = []
velocity_distributionss = []
for i, file_bases in enumerate(file_basess_velocity_distribution):
    velocity_bin_grids = []
    velocity_distributions = []

    path_parent = paths_parent[i]

    for j, file_base in enumerate(file_bases):
        file_read = path_read / 'velocity_distribution_{0}_{1}_2025_09_30.txt'.format(i, j)
        velocity_bin_grid, velocity_distribution = np.loadtxt(file_read, dtype=float, delimiter=" ")

        velocity_bin_grids.append(velocity_bin_grid)
        velocity_distributions.append(velocity_distribution/(np.sum(velocity_distribution) * (velocity_bin_grid[1] - velocity_bin_grid[0])))
    
    velocity_bin_gridss.append(velocity_bin_grids)
    velocity_distributionss.append(velocity_distributions)
# %%
'''
Read particle flux distributions
'''
N_gridss = []
flux_distributionss = []

for i, file_bases in enumerate(file_basess_flux_distribution):
    flux_distributions = []
    N_grids = []

    path_parent = paths_parent[i]

    for j, file_base in enumerate(file_bases):
        file_read = path_parent / file_base
        fluxes = read_mrmd_out(file_read, usecols=(9))
        N_grid = np.arange(np.min(fluxes), np.max(fluxes), 2)
        fluxHist = np.histogram(fluxes, N_grid, density=True)

        N_grids.append((N_grid[:-1] + 0.5 * (N_grid[1] - N_grid[0])) * conversion_factors[j])
        flux_distributions.append(fluxHist[0]/conversion_factors[j])

    N_gridss.append(N_grids)
    flux_distributionss.append(flux_distributions)
# %%
'''
Read average interface forces 
'''
x_gridss = []
average_force_profiless = []

for i, file_bases in enumerate(file_basess_Fav):
    x_grids = []
    average_force_profiles = []

    path_parent = paths_parent[i]

    for j, file_base in enumerate(file_bases):
        file_read = path_parent / '{0}/{0}_Fav_data.txt'.format(file_base)
        x_grid, average_force_profile, _, __ = np.loadtxt(file_read, dtype=float).T

        x_grids.append(x_grid - x_inters[j])
        average_force_profiles.append(average_force_profile)
    
    x_gridss.append(x_grids)
    average_force_profiless.append(average_force_profiles)
# %%
pos_gridsss = [velocity_bin_gridss, N_gridss, x_gridss]
data_profilesss = [velocity_distributionss, flux_distributionss, average_force_profiless]
# %%
legend_loc = (0.8, 0.85)

x_labelss = []
y_labelss = []
y_label_types = [r'$f_{\mathrm{M}}(p; x_{\mathrm{AT}/\Delta})$ in $\tau/(M\sigma)$', r'$P(J)/\sqrt{A}$ in $\sigma \tau$', r'$F_{\mathrm{av}}^x(x_{\mathrm{AT}/\Delta} - x)$ in $\epsilon/\sigma$']
x_label_types = [r'$p$ in $M\sigma/\tau$', r'$J \cdot \sqrt{A}$ in $1/(\sigma \tau)$', r'$x - x_{\mathrm{AT}/\Delta}$ in $\sigma$']
for i, pos_gridss in enumerate(pos_gridsss):
    x_labels = []
    y_labels = []
    for j in range(len(pos_gridss)):
        x_labels.append(x_label_types[i])
        if j == 0:
            y_labels.append(y_label_types[i])
        else:
            y_labels.append(None)
    x_labelss.append(x_labels)
    y_labelss.append(y_labels)

x_lims_types = [[0.0, 10.0], [-0.15, 0.15], [-2.5, 2.5]]
x_limsss = [[x_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

y_lims_types = [[0.0, 0.5], [0.0, 16.0], [-2.0, 1.6]]
y_limsss = [[y_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

data_labels = ["AdResS (minimal)",\
               "AdResS (no AT)",\
               "AdResS (complete)",\
               "Full-AT (complete)"]

colors = [select_color_from_colormap(range_indicator/3 / 2, 'plasma') for range_indicator in range(3)]
colors.append('black')
linestyles = ['--', '-.', ':', '-']

axis_titles = [r'$\rho = \SI{0.296}{\sigma^{-3}}$, $T = \SI{1.5}{\epsilon/k_{\mathrm{B}}}$',\
               r'$\rho = \SI{0.370}{\sigma^{-3}}$, $T = \SI{1.5}{\epsilon/k_{\mathrm{B}}}$',\
               r'$\rho = \SI{0.370}{\sigma^{-3}}$, $T = \SI{2.0}{\epsilon/k_{\mathrm{B}}}$']
#axis_titles = ["a", "b", "c"]
scale = 5
# %%
fig, ax = plt.subplots(len(pos_gridsss), len(pos_gridsss[0]), figsize=(len(pos_gridsss[0]) * scale, len(pos_gridsss) * scale), constrained_layout=True, squeeze=False, )

for column_index in range(0, len(pos_gridsss[i])):
    pos_grids = pos_gridsss[0][column_index]
    data_profiles = data_profilesss[0][column_index]
    x_label = x_labelss[0][column_index]      
    y_label = y_labelss[0][column_index]
    x_lims = x_limsss[0][column_index]
    y_lims = y_limsss[0][column_index]
    if axis_titles and 0 == 0:
        axis_title = axis_titles[column_index]
    else:
        axis_title = None

    plt.sca(ax[0, column_index])
    
    for pos_grid, data_profile, data_label, color, linestyle in zip(pos_grids, data_profiles, data_labels, colors, linestyles):
        plt.plot(pos_grid, data_profile, label=data_label, color=color, linewidth=0.75, linestyle=linestyle)
    
    plt.xlim(x_lims)
    plt.ylim(y_lims)
    plt.title(axis_title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    #if row_index < len(pos_gridsss) - 1:
    #    plt.tick_params(axis='x', labelbottom=False)
    #if row_index == 0:
    #    plt.tick_params(axis='x', labeltop=True)
    if column_index > 0:
        plt.tick_params(axis='y', labelleft=False)
        
    #if column_index == len(pos_gridsss[0]) - 1:
    #    plt.tick_params(axis='y', labelright=True)

handles = [mlines.Line2D([], [], color=colors[i], label=data_labels[i], linestyle=linestyles[i]) for i in range(len(data_labels))]

fig.legend(handles=handles, loc=legend_loc, framealpha=0.0)
for column_index in range(0, len(pos_gridsss[1])):
    pos_grids = pos_gridsss[1][column_index]
    data_profiles = data_profilesss[1][column_index]
    x_label = x_labelss[1][column_index]      
    y_label = y_labelss[1][column_index]
    x_lims = x_limsss[1][column_index]
    y_lims = y_limsss[1][column_index]
    if axis_titles and 1 == 0:
        axis_title = axis_titles[column_index]
    else:
        axis_title = None

    plt.sca(ax[1, column_index])
    
    plot_PofNs(pos_gridsss[1][column_index], data_profilesss[1][column_index], fig_title=None, colors=colors, x_lims=x_lims, y_lims=y_lims, data_labels=data_labels, legend=None, xlabel=x_label, ylabel=y_label, linestyles=linestyles)
    
    plt.xlim(x_lims)
    plt.ylim(y_lims)
    plt.title(axis_title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    #if row_index < len(pos_gridsss) - 1:
    #    plt.tick_params(axis='x', labelbottom=False)
    #if row_index == 0:
    #    plt.tick_params(axis='x', labeltop=True)
    if column_index > 0:
        plt.tick_params(axis='y', labelleft=False)
        
    #if column_index == len(pos_gridsss[0]) - 1:
    #    plt.tick_params(axis='y', labelright=True)

handles = [mlines.Line2D([], [], color=colors[i], label=data_labels[i], linestyle=linestyles[i]) for i in range(len(data_labels))]
fig.legend(handles=handles, loc=legend_loc, framealpha=0.0)

for column_index in range(0, len(pos_gridsss[2])):
    pos_grids = pos_gridsss[2][column_index]
    data_profiles = data_profilesss[2][column_index]
    x_label = x_labelss[2][column_index]      
    y_label = y_labelss[2][column_index]
    x_lims = x_limsss[2][column_index]
    y_lims = y_limsss[2][column_index]
    if axis_titles and 2 == 0:
        axis_title = axis_titles[column_index]
    else:
        axis_title = None

    plt.sca(ax[2, column_index])
    
    plot_profile_at_interface(pos_grids, data_profiles, 0.0, x_label=x_label, y_label=y_label, data_labels=data_labels, fig_title=None, colors=colors, x_lims=x_lims, y_lims=y_lims, legend=None, label_left='AT')

    #if row_index < len(pos_gridsss) - 1:
    #    plt.tick_params(axis='x', labelbottom=False)
    #if row_index == 0:
    #    plt.tick_params(axis='x', labeltop=True)
    if column_index > 0:
        plt.tick_params(axis='y', labelleft=False)
        
    #if column_index == len(pos_gridsss[0]) - 1:
    #    plt.tick_params(axis='y', labelright=True)

handles = [mlines.Line2D([], [], color=colors[i], label=data_labels[i], linestyle=linestyles[i]) for i in range(len(data_labels))]
fig.legend(handles=handles, loc=legend_loc, framealpha=0.0)


if file_plt:
    if os.path.isfile(file_plt):
        print('file already exists! Aborting.')
    else:
        plt.savefig(fname=file_plt, dpi=150, format=format_plt)

# %%
