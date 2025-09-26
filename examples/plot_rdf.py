# %%
import pathlib as plb

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.structure import calc_rdf_in_slab
# %%
density = 0.370

file_base = "tracerProduction_n370_T20_24019679_2025_09_23"

#path_parent = plb.Path("/media/julianhille/T7 Shield/backups_unpacked/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations").resolve()
path_parent = plb.Path("/home/mi/julianhille/mounted_directories/curta/project/mrmd_simulations").resolve()

file_trj = path_parent / 'rho{0:04d}_T20_2025_09_19/{1}/{1}.h5md'.format(int(density * 1000), file_base)
file_top = path_parent / 'rho{0:04d}_T20_2025_09_19/{1}/{1}.gro'.format(int(density * 1000), file_base)

format_trj = 'H5MD'
format_top = 'GRO'

atom_group_1_name = 'Ar'
atom_group_2_name = 'Ar'
x_min = 00
x_max = 100
nbins = 1000
range = [0.5, 15.0]
t_init = 0
t_fin = -1
t_step = 1
# %%
rdf = calc_rdf_in_slab(file_top, file_trj, format_top, format_trj, atom_group_1_name, atom_group_2_name, x_min, x_max, nbins, range, t_init, t_fin, t_step)
# %%
title = '{0}-{1} RDFs in slab '.format(atom_group_1_name, atom_group_2_name) + r'$x_{min} = ' + '{0}'.format(x_min) + r' \AA$ to $x_{max} = ' + '{0}'.format(x_max) + r' \AA$'

plt.plot(rdf.results.bins, rdf.results.rdf)

plt.xlabel(r'$r$ in $\AA$')
plt.ylabel(r'$g(r)$ in a.u.')
plt.legend()
plt.xlim(0.0, 2.5)
plt.title(title)
#plt.savefig(out_file, format='png', dpi=300)
plt.show()
# %%
