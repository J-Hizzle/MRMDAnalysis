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
file_bases = ["tracerProduction_23945717_2025_09_09",\
              "tracerProduction_23945715_2025_09_09",\
              "tracerProduction_23945712_2025_09_09",\
              "atomisticProduction_23945334_2025_09_09"]
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

#subdomain_dimensionss = [(10.0, 10.0, 10.0),\
#                        (10.0, 20.0, 20.0),\
#                        (10.0, 30.0, 30.0),\
#                        (30.0, 30.0, 30.0),
#                        (30.0, 30.0, 30.0)]
subdomain_dimensionss = [(10.0, 20.0, 20.0),\
                        (10.0, 30.0, 30.0),\
                        (30.0, 30.0, 30.0),
                        (30.0, 30.0, 30.0)]

density = 0.370
#density = 0.296
#density = 0.247
#density = 0.198

files_obs = [next(plb.Path('/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/rho{0:04d}_2025_09_03/{1}/'.format(int(density * 1000), file_base)).resolve().glob("*.out")) for file_base in file_bases]
interface_areas = [subdomain_dimensions[1] for subdomain_dimensions in subdomain_dimensionss] # in sigma^2
time_interval = 1 # in tau 
conversion_factors = [10/(interface_area * time_interval) for interface_area in interface_areas]

initial_guesses = [np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100]), np.asarray([0.01, 0, 100])]
fit_xlims = [-50.0, 50.0]

fig_title_flux = "Boundary particle flux distributions"
fig_title_fit = "Gaussian-fit boundary particle flux distributions"
xlabel = r'$\delta N$ in $1$'
ylabel = r'$P(\delta N)$ in $1$'

labels = ["$x = {0}, y = {1}, z = {2}$".format(subdomain_dimensions[0], subdomain_dimensions[1], subdomain_dimensions[2]) for subdomain_dimensions in subdomain_dimensionss]
labels[-1] = "$x = {0}, y = {1}, z = {2}$, Full-AT".format(subdomain_dimensionss[-1][0], subdomain_dimensionss[-1][1], subdomain_dimensionss[-1][2])

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
    N_grid = np.arange(np.min(fluxes) * conversion_factors[i], np.max(fluxes) * conversion_factors[i], 1)
    fluxHist = np.histogram(fluxes * conversion_factors[i], N_grid, density=True)

    N_grids.append(N_grid)
    flux_distributions.append(fluxHist[0])
# %%
plot_PofNs(N_grids, flux_distributions, labels, colors, fig_title_flux, xlabel=xlabel, ylabel=ylabel, file_out=file_plt_flux, fontsize=fontsize)
# %%
plot_PofNs_fit(N_grids, flux_distributions, initial_guesses, fig_title_fit, fit_xlims, labels, colors, markers, xlabel=xlabel, ylabel=ylabel, file_out=file_plt_fit, fontsize=fontsize)
# %%
