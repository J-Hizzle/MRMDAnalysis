# %%
import pathlib as plb

import numpy as np
import matplotlib.pyplot as plt
# %%
file_bases = ["tracerProduction_23828863_2025_08_15",\
              "tracerProduction_23843205_2025_08_19",\
              "tracerProduction_23828869_2025_08_15"]
files_obs = [plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations/{0}/{0}_observables.txt'.format(file_base)).resolve() for file_base in file_bases]

labels = ["$N = 370$", "$N = 3330$", "$N = 15000$"]
conversion_factors = [3, 1, 1]
# %%
for i, file_obs in enumerate(files_obs):
    with open(file_obs, 'r') as file:
        lines = map(lambda x: x.replace("│", ""), file.readlines())
        fluxes = np.loadtxt(lines, skiprows=2, usecols=(8), dtype=int)
    N_grid = np.arange(np.min(fluxes), np.max(fluxes), 1) * conversion_factors[i]
    fluxHist = np.histogram(fluxes * conversion_factors[i], N_grid, density=True)
    print(i, N_grid * conversion_factors[i])
    plt.plot(N_grid[:-1], fluxHist[0], label=labels[i])

plt.legend()
# %%
