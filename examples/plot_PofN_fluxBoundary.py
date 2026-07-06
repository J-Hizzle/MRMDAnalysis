# %%
import numpy as np
import pathlib as plb
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
plt.rc('text.latex', preamble=r'\usepackage{siunitx}')

from mrmdanalysis.plot import plot_PofNs

from mrmdanalysis.read import read_mrmd_out
# %%
file_read_NVE = plb.Path("/srv/public/julianhille/project/mrmd/build/examples/06_OpenBoundary/fluxBoundary_NVE.out").resolve()
data_NVE = read_mrmd_out(file_read_NVE, usecols=(7))

file_read_NVT = plb.Path("/srv/public/julianhille/project/mrmd/build/examples/06_OpenBoundary/fluxBoundary.out").resolve()
data_NVT = read_mrmd_out(file_read_NVT, usecols=(7))
# %%
N_grid_NVE = np.arange(np.min(data_NVE), np.max(data_NVE), 1)
hist_NVE = np.histogram(data_NVE, N_grid_NVE, density=True)

N_grid_NVT = np.arange(np.min(data_NVT), np.max(data_NVT), 1)
hist_NVT = np.histogram(data_NVT, N_grid_NVT, density=True)
#  %%
plt.plot(N_grid_NVE[:-1], hist_NVE[0], drawstyle='steps-mid', label='NVE')
plt.plot(N_grid_NVT[:-1], hist_NVT[0], drawstyle='steps-mid', label='NVT')
plt.legend()
plt.xlabel(r'$N$')
plt.ylabel(r'$P(N)$')
plt.title(r'2-sided Kinetic Chemostat connected to ideal gas')
plt.savefig('data/PofN_fluxBoundary_2026_07_06.png', dpi=300, format='png')
# %%
