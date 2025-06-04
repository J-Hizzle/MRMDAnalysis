# %%
import pathlib as plb

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.structure import calc_rdf_in_slab
# %%
path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce').resolve()
#path_data = plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mosceto/testing/ghost_face_lj_fa_prod_26_03_2024').resolve()
file_top = path_data / 'equilibrateLangevin.gro'
file_trj = path_data / 'equilibrateLangevin.h5md'
format_top = 'GRO'
format_trj = 'H5MD'
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
