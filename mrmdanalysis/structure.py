# %%
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
# %%
def calc_rdf_in_slab(file_top : str, file_trj : str, format_top : str, format_trj : str, atom_group_1_name : str, atom_group_2_name : str, x_min : float, x_max : float, nbins : int, range : list, t_init : int = 0, t_fin : int = -1, t_step : int = 1):
    universe = mda.Universe(file_top, file_trj, topology_format=format_top, format=format_trj, convert_units=False)
    atom_group_1 = universe.select_atoms('name {0} and prop x > {1} and prop x < {2}'.format(atom_group_1_name, x_min, x_max), updating=True)
    atom_group_2 = universe.select_atoms('name {0} and prop x > {1} and prop x < {2}'.format(atom_group_2_name, x_min, x_max), updating=True)
    rdf = InterRDF(atom_group_1, atom_group_2, nbins=nbins, range=range)
    rdf.run(t_init, t_fin, t_step, verbose=True)
    return rdf
# %%