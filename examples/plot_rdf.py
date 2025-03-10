# %%
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

import pathlib as plb

import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
#path_data = plb.Path('/srv/public/julianhille/project/mrmd/build/examples/TracerThermodynamicForce').resolve()
path_data = plb.Path('/home/mi/julianhille/mounted_directories/curta/project/mosceto/testing/ghost_face_lj_fa_prod_26_03_2024').resolve()
file_top = path_data / 'confout.gro'
file_trj = path_data / 'confout.gro'
format_top = 'GRO'
format_trj = 'GRO'
atom_group_1_name = 'Ar'
atom_group_2_name = 'Ar'
x_min = 00
x_max = 450
nbins = 200
range = [0.5, 15.0]
t_init = 0
t_fin = 10
t_step = 1
# %%
def calc_rdf_in_slab(file_top : str, file_trj : str, format_top : str, format_trj : str, atom_group_1_name : str, atom_group_2_name : str, x_min : float, x_max : float, nbins : int, range : list, t_init : int = 0, t_fin : int = -1, t_step : int = 1):
    universe = mda.Universe(file_top, file_trj, topology_format=format_top, format=format_trj, convert_units=False)
    atom_group_1 = universe.select_atoms('name {0} and prop x > {1} and prop x < {2}'.format(atom_group_1_name, x_min, x_max), updating=True)
    atom_group_2 = universe.select_atoms('name {0} and prop x > {1} and prop x < {2}'.format(atom_group_2_name, x_min, x_max), updating=True)
    rdf = InterRDF(atom_group_1, atom_group_2, nbins=nbins, range=range)
    rdf.run(t_init, t_fin, t_step, verbose=True)
    return rdf
# %%
rdf = calc_rdf_in_slab(file_top, file_trj, format_top, format_trj, atom_group_1_name, atom_group_2_name, x_min, x_max, nbins, range, t_init, t_fin, t_step)
# %%
title = '{0}-{1} RDFs in slab '.format(atom_group_1_name, atom_group_2_name) + r'$x_{min} = ' + '{0}'.format(x_min) + r' \AA$ to $x_{max} = ' + '{0}'.format(x_max) + r' \AA$'

plt.plot(rdf.results.bins, rdf.results.rdf)

plt.xlabel(r'$r$ in $\AA$')
plt.ylabel(r'$g(r)$ in a.u.')
plt.legend()
plt.title(title)
#plt.savefig(out_file, format='png', dpi=300)
plt.show()
# %%
