# %%
import numpy as np
import pathlib as plb

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.plot import plot_PofNs, plot_PofNs_fit, plot_Noft
from mrmdanalysis.statistics import count_particle_numbers, create_particle_number_hist
# %%
path_data = plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations').resolve()
file_bases = ["equilibrateLangevin_x45_y30_z30_n15000_23108150_2025_06_04", "tracerProduction_23145081_2025_06_10"]
files_top = [path_data / '{0}/{0}.gro'.format(file_base) for file_base in file_bases]
#files_top = [path_data / '{0}/{0}.gro'.format(file_base) for file_base in file_bases]
files_trj = [path_data / '{0}/{0}.h5md'.format(file_base) for file_base in file_bases]
r_refs = [22.5, 5.0]
width = 5.0
atom_selections = ['name Ar and prop x > {0} and prop x < {1}'.format(r_ref - width/2, r_ref + width/2) for r_ref in r_refs]
formats_top = ['GRO' for file_base in file_bases]
formats_trj = ['H5MD' for file_base in file_bases]
labels = [r'full-atomistic', r'minimal AdResS']
colors = ['blue', 'red']
markers = ['x', '+']
frame_freqs = [1 for file_base in file_bases]
frame_limss = [[0, -1] for file_base in file_bases]
initial_guesses = [np.asarray([0.01, 1200, 100]), np.asarray([0.01, 1200, 100])]
fit_xlims = [-150.0, 150.0]

N_grid = np.arange(1250.0, 2000.0, 10.0)
step_time = 2 # in tau
total_steps = 10000000 # total simulated steps
total_time = step_time * total_steps# in tau
save_step = 10000
save_time = step_time * save_step # in tau

time_grids = [np.arange(0, total_time, save_time) for file_base in file_bases]
run_mod = 10

hist_title = r"Particle number distributions (width $= " + str(width) + r"\tau$)"
fit_title = "Shifted and Gaussian-fitted P(N)'s (width $= " + str(width) + r"\tau$)"
Noft_title = "Number of particles in AT over time (width $= " + str(width) + r"\tau$)"

hist_out = None
fit_out = None
Noft_out = None

path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_06_10').resolve()
hist_out = path_out / 'PofN_w{0}_2025_06_10.png'.format(width)
fit_out = path_out / 'PofN_fit_w{0}_2025_06_10.png'.format(width)
Noft_out = path_out / 'Noft_w{0}_2025_06_10.png'.format(width)
# %%
# calculate numbers of particles and construct histograms
particle_numberss = []
P_valss = []

for i, file_top in enumerate(files_top):
    print('file number =', i)
    file_trj = files_trj[i]
    format_top = formats_top[i]
    format_trj = formats_trj[i]
    label = labels[i]
    color = colors[i]
    marker = markers[i]
    frame_freq = frame_freqs[i]
    frame_lims = frame_limss[i]
    atom_selection = atom_selections[i]

    particle_numbers = count_particle_numbers(file_top, file_trj, atom_selection, format_top, format_trj, frame_freq, frame_lims)

    P_vals = create_particle_number_hist(particle_numbers, N_grid)

    particle_numberss.append(np.copy(particle_numbers))
    P_valss.append(np.copy(P_vals))
# %%
plot_PofNs(N_grid, P_valss, labels, colors, markers, hist_title, file_out=hist_out)
# %%
plot_PofNs_fit(N_grid, P_valss, initial_guesses, fit_title, fit_xlims, labels, colors, markers, file_out=fit_out)
# %%
plot_Noft(time_grids, particle_numberss, labels, colors, markers, Noft_title, file_out=Noft_out, run_mod=run_mod)
# %%
