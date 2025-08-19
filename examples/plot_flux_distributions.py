# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.plot import plot_PofNs_fit, select_color_from_colormap
# %%
file_bases = ["tracerProduction_23828863_2025_08_15",\
              "tracerProduction_23843205_2025_08_19",\
              "tracerProduction_23828869_2025_08_15"]
files_obs = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/{0}/{0}_observables.txt'.format(file_base)).resolve() for file_base in file_bases]

conversion_factors = [3, 1, 1]

initial_guesses = [np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100])]
fit_xlims = [-50.0, 50.0]

fig_title = "Boundary particle flux distributions"
xlabel = r'$\delta N$ in $1$'
ylabel = r'$P(\delta N)$ in $1$'

labels = ["$N = 370$", "$N = 3330$", "$N = 15000$"]
markers = ["^", "o", "s"]
colors = [select_color_from_colormap(range_indicator/len(file_bases), 'plasma') for range_indicator in range(len(file_bases))]

file_plt_flux = None

out_base = "boundary_flux_distributions"
path_out = plb.Path('/srv/public/julianhille/presentations/reports_2025/report_2025_08_26').resolve()
file_plt_flux = path_out / '{0}.png'.format(out_base)
# %%
flux_distributions = []
N_grids = []

for i, file_obs in enumerate(files_obs):
    with open(file_obs, 'r') as file:
        lines = map(lambda x: x.replace("│", ""), file.readlines())
        fluxes = np.loadtxt(lines, skiprows=2, usecols=(8), dtype=int)
    N_grid = np.arange(np.min(fluxes), np.max(fluxes), 1) * conversion_factors[i]
    fluxHist = np.histogram(fluxes * conversion_factors[i], N_grid, density=True)

    N_grids.append(N_grid)
    flux_distributions.append(fluxHist[0])
# %%
plot_PofNs_fit(N_grids, flux_distributions, initial_guesses, fig_title, fit_xlims, labels, colors, markers, xlabel=xlabel, ylabel=ylabel, file_out=file_plt_flux)
# %%
