# %%
import pathlib as plb
import numpy as np
from mrmdanalysis.plot import plot_profile_array_in_adress_box, select_color_from_colormap
from mrmdanalysis.util import calc_coupling_potential_from_thermodynamic_force, calc_applied_thermodynamic_force

from ghostface.LJ import Pot_lj_shifted_cap, Rdf_lj_goldman, Rdf_lj_morsali, Rdf_lj_read
from ghostface.interaction_correction.interaction_correction import calc_av_force

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
plt.rc('text.latex', preamble=r'\usepackage{siunitx}')
# %%
file_basess = [["equilibrateLangevin_rho296_T15_x10_yz20_2026_06_03",\
                "equilibrateLangevin_rho296_T15_x10_yz30_2026_06_03",\
                "equilibrateLangevin_rho296_T15_x30_yz30_2026_06_03"],\
               ["equilibrateLangevin_rho370_T15_x10_yz20_2026_06_03",\
                "equilibrateLangevin_rho370_T15_x10_yz30_2026_06_03",\
                "equilibrateLangevin_rho370_T15_x30_yz30_2026_06_03"],\
               ["equilibrateLangevin_rho370_T20_x10_yz20_2026_06_03",\
                "equilibrateLangevin_rho370_T20_x10_yz30_2026_06_03",\
                "equilibrateLangevin_rho370_T20_x30_yz30_2026_06_03"]]

file_bases_parametrized = ["goldman", "morsali"]

densities = [0.296,\
             0.370,\
             0.370]

temperatures = [1.5,\
                1.5,\
                2.0]

path_parent = plb.Path("/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/noATsimulations_2026_06_10/equilibrateLangevin_2026_06_03").resolve()
#path_parent = plb.Path("/media/julianhille/T7 Shield/backups_noATPaper_2026_04_29/equilibrateLangevin").resolve()

format_trj = 'H5MD'
format_top = 'GRO'

atom_group_1_name = 'Ar'
atom_group_2_name = 'Ar'

densBinWidth = 0.01
densShift = densBinWidth/2.0

forceGridSpacing = 0.01
forceShift = forceGridSpacing/2.0

r_cut = 2.5
r_cap = 0.82417464
pot_lj_shifted_cap = Pot_lj_shifted_cap(r_cut, r_cap)

r_ref_base = 5.0
r_min_base = 0.0
r_max_base = 2.5
app_min_base = 0.0
app_max_base = 4.5

app_min = 0.0
app_max = 5.0
r_ref = 5.0
x_delta_TR_initialGuess = 2.5
x_grid_initialGuess = np.arange(0.0, app_max, densBinWidth) + densShift

rdf_range = [0.01, 2.5]
nbins = 100

t_init = 0
t_fin = -1
t_step = 1

r_grid = np.linspace(0.01, 2.5, 10000)

file_plt = None
format_plt = None

path_out = plb.Path("/srv/public/julianhille/publications/paper_noAT").resolve()
out_base = "tf_apriori_array_2026_06_30.svg"
file_plt = path_out / '{0}'.format(out_base)
format_plt = 'SVG'
# %%
'''
Calculate the a priori thermodynamic forces
'''
rdfss = []
for i, file_bases in enumerate(file_basess):
    rdfs = []
    rdfs_eval = []
    for j, file_base in enumerate(file_bases):
        file_top = path_parent / '{0}.gro'.format(file_base)
        file_trj = path_parent / '{0}.h5md'.format(file_base)
    
        rdf_lj_read = Rdf_lj_read(file_top, file_trj, format_top, format_trj, atom_group_1_name, atom_group_2_name, nbins, rdf_range, t_init, t_fin, t_step)
        color = select_color_from_colormap(j/len(file_bases), 'plasma') 
        rdfs.append(rdf_lj_read)
    
    rdf_lj_goldman = Rdf_lj_goldman(dens=densities[i], temp=temperatures[i])
    rdf_lj_morsali = Rdf_lj_morsali(dens=densities[i], temp=temperatures[i])

    rdfs.append(rdf_lj_goldman)
    rdfs.append(rdf_lj_morsali)
    rdfss.append(rdfs)

rdfss_eval = []
for i, rdfs in enumerate(rdfss):
    rdfs_eval = []
    for j, rdf in enumerate(rdfs):
        rdf_eval = rdf.rdf_lj(r_grid)
        rdfs_eval.append(rdf_eval)
    rdfss_eval.append(rdfs_eval)

for i, rdfs in enumerate(rdfss):
    for j, rdf in enumerate(rdfs):
        force_vals_av_initialGuess, _ = calc_av_force(x_grid_initialGuess[(x_grid_initialGuess >= app_min) & (x_grid_initialGuess <= app_max)], x_delta_TR_initialGuess, r_cut, densities[i], rdf.rdf_lj, pot_lj_shifted_cap.force_lj)

        # fill in force values from initial guess
        force_vals_initialGuess = np.zeros_like(x_grid_initialGuess)
        force_vals_initialGuess[(x_grid_initialGuess >= app_min) & (x_grid_initialGuess <= app_max)] = force_vals_av_initialGuess
        # symmetrize
        x_grid_initTotal, indices = np.unique(np.concatenate([r_ref - np.flip(x_grid_initialGuess), r_ref + x_grid_initialGuess]), return_index=True)
        force_vals_initTotal = np.concatenate([-np.flip(force_vals_initialGuess), force_vals_initialGuess])[indices]
        #applied_force_vals_av = calc_applied_thermodynamic_force(x_grid_initTotal, force_vals_initTotal, 5.0, 0.0, 2.5)
        #applied_pot_vals_av = calc_coupling_potential_from_thermodynamic_force(x_grid_initTotal, applied_force_vals_av, 0.0, 2.5)
        applied_force_vals_av = calc_applied_thermodynamic_force(x_grid_initTotal, force_vals_initTotal, r_ref, app_min, x_delta_TR_initialGuess)
        applied_pot_vals_av = calc_coupling_potential_from_thermodynamic_force(x_grid_initTotal, applied_force_vals_av, r_ref, app_min)
        # interpolate 
        x_grid_interp = np.arange(0.0, 2.0 * r_ref, forceGridSpacing) + forceShift
        applied_force_vals_av_interp = np.interp(x_grid_interp, x_grid_initTotal, applied_force_vals_av)

        write_array = np.column_stack((x_grid_interp, applied_force_vals_av_interp)).T
        if j > len(file_basess[i]) - 1:
            savename = path_parent / '{0}_rho{1:.3}_T{2}_initialGuess.txt'.format(file_bases_parametrized[j - len(file_basess[i])], str(densities[i]).lstrip('0.'), str(temperatures[i]).replace('.', ''))
        else:
            savename = path_parent / '{0}_initialGuess.txt'.format(file_basess[i][j])

        np.savetxt(savename, write_array, fmt='%05.6e')
# %%
x_gridss_interp = []
tf_profiless = []
for i in range(len(file_basess)):
    x_grids_interp = []
    tf_profiles = []
    for j in range(len(file_basess[0]) + 2):
        if j > len(file_basess[i]) - 1:
            readname = path_parent / '{0}_rho{1:.3}_T{2}_initialGuess.txt'.format(file_bases_parametrized[j - len(file_basess[i])], str(densities[i]).lstrip('0.'), str(temperatures[i]).replace('.', ''))
        else:
            readname = path_parent / '{0}_initialGuess.txt'.format(file_basess[i][j])

        x_grid_interp, applied_force_vals_av_interp = np.loadtxt(readname, dtype=float, delimiter=" ")
        x_grids_interp.append(x_grid_interp)
        tf_profiles.append(applied_force_vals_av_interp)

    x_gridss_interp.append(x_grids_interp)
    tf_profiless.append(tf_profiles)
# %%
pos_gridsss = [x_gridss_interp]
data_profilesss = [tf_profiless]
# %%
legend_loc = (0.85, 0.15)

x_labelss = []
y_labelss = []
y_label_types = [r'$F_{\mathrm{av}}^x$ in $\epsilon/\sigma$']

for i, pos_gridss in enumerate(pos_gridsss):
    x_labels = []
    y_labels = []
    for j in range(len(pos_gridss)):
        if i == len(pos_gridsss) - 1:
            x_labels.append(r'$x - x_{\mathrm{AT}/\Delta}$ in $\sigma$')
        else:
            x_labels.append(None)
        if j == 0:
            y_labels.append(y_label_types[i])
        else:
            y_labels.append(None)
    x_labelss.append(x_labels)
    y_labelss.append(y_labels)

x_lims_types = [[0.0, 5.0]]
x_limsss = [[x_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

y_lims_types = [[-2.0, 1.6]]
y_limsss = [[y_lims_types[i] for j in range(len(pos_gridss))] for i, pos_gridss in enumerate(pos_gridsss)]

data_labels = ["EQ (minimal)",\
               "EQ (no AT)",\
               "EQ (complete)",\
               "Goldman",\
               "Morsali"]

axis_titles = [r'$\rho = \SI{0.296}{\sigma^{-3}}$, $T = \SI{1.5}{\epsilon/k_{\mathrm{B}}}$',\
               r'$\rho = \SI{0.370}{\sigma^{-3}}$, $T = \SI{1.5}{\epsilon/k_{\mathrm{B}}}$',\
               r'$\rho = \SI{0.370}{\sigma^{-3}}$, $T = \SI{2.0}{\epsilon/k_{\mathrm{B}}}$']

colors = [select_color_from_colormap(range_indicator/(len(tf_profiles) - 2) / 2, 'plasma') for range_indicator in range(len(tf_profiles))]
colors[-2] = 'red'
colors[-1] = 'green'
linestyles = ['--', '-.', ':', '-', '-'] 


scale = 5
# %%
plot_profile_array_in_adress_box(pos_gridsss, data_profilesss, x_limsss, y_limsss, x_labelss, y_labelss, data_labels, colors, r_ref_base, r_min_base, r_max_base, linestyles=linestyles, file_plt=file_plt, format_plt=format_plt, legend_loc=legend_loc, scale=scale, axis_titles=axis_titles)
# %%
