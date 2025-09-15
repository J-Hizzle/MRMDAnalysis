# %%
import pathlib as plb
import glob

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.plot import plot_PofNs, plot_PofNs_fit, select_color_from_colormap
from mrmdanalysis.read import read_mrmd_out
# %%
file_bases = ["tracerProduction_23945719_2025_09_09",\
              "tracerProduction_23945717_2025_09_09",\
              "tracerProduction_23945715_2025_09_09",\
              "tracerProduction_23945712_2025_09_09",\
              "atomisticProduction_23945334_2025_09_09"]

density = 0.370
#density = 0.198

files_obs = [next(plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho{0:04d}_2025_09_03/{1}/'.format(int(density * 1000), file_base)).resolve().glob("*.out")) for file_base in file_bases]

conversion_factors = [3, 1.5, 1, 1, 1, 1]

initial_guesses = [np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100])]
fit_xlims = [-50.0, 50.0]

fig_title_flux = "Boundary particle flux distributions"
fig_title_fit = "Gaussian-fit boundary particle flux distributions"
xlabel = r'$\delta N$ in $1$'
ylabel = r'$P(\delta N)$ in $1$'

labels = [r"AdResS ($N = 370$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 3330$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 9990$, $\gamma = 20$, $d_{TR} = 5$)", r"AdResS ($N = 15000$, $\gamma = 20$, $d_{TR} = 34$)", r"Full-AT ($N = 9990$, $\gamma = 20$)", r"Full-AT ($N = 9990$, $\gamma = 1$)"]
markers = ["^", "o", "v", "s", "D", "P"]
colors = [select_color_from_colormap(range_indicator/len(file_bases), 'plasma') for range_indicator in range(len(file_bases))]
fontsize = 8

file_plt_flux = None

out_base = "boundary_flux_distributions_test_2025_09_12"
path_out = plb.Path('/srv/public/julianhille/project/MRMDAnalysis/data').resolve()
file_plt_flux = path_out / '{0}.png'.format(out_base)
file_plt_fit = path_out / '{0}_fit.png'.format(out_base)
# %%
flux_distributions = []
N_grids = []

for i, file_obs in enumerate(files_obs):
    print("file number = ", i)
    fluxes = read_mrmd_out(file_obs, usecols=(9))
    N_grid = np.arange(np.min(fluxes), np.max(fluxes), 1) * conversion_factors[i]
    fluxHist = np.histogram(fluxes * conversion_factors[i], N_grid, density=True)

    N_grids.append(N_grid)
    flux_distributions.append(fluxHist[0])
# %%
plot_PofNs(N_grids, flux_distributions, labels, colors, fig_title_flux, xlabel=xlabel, ylabel=ylabel, file_out=file_plt_flux, fontsize=fontsize)
# %%
plot_PofNs_fit(N_grids, flux_distributions, initial_guesses, fig_title_fit, fit_xlims, labels, colors, markers, xlabel=xlabel, ylabel=ylabel, file_out=file_plt_fit, fontsize=fontsize)
# %%
