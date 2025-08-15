# %%
import numpy as np
import MDAnalysis as mda
# %%
def create_particle_number_hist(particle_numbers, N_grid):
    hist, _ = np.histogram(particle_numbers, N_grid, density=False)

    return hist
# %%
def count_particle_numbers(file_top, file_trj, atom_selection, format_top, format_trj, frame_freq, frame_lims):
    '''
    Parameters:
    -----------
        file_top : str
            Topology binary file.
        file_trj : str
            Trajectory file.
        atom_selection: str
            Atom selection including the type and the boundaries of the open system in x-direction.
        format_top : str (Options : ['TPR', 'GRO'])
        format_trj : str (Options : ['XTC', 'TRJ', 'H5MD'])
            File format for the trajectory file.
        frame_freq : int
            Frame frequency - every 'frame_freq'-th frame is used in the calculation of the particle numbers.
            (Useful if too many frames in trajectory)
        frame_lims : list (dtype : int)
            Limitation for the frames to use in the plot.
            (Useful if one needs to match different simulations to each other)
        label : str
            Label for this plot.

    Returns:
    --------
        particle_numbers : np.ndarray (dtype=int)
            
    '''
    # load files into mda universe
    universe = mda.Universe(file_top, file_trj, format=format_trj, topology_format=format_top, convert_units=False)
    atom_group = universe.select_atoms(atom_selection, updating=True)

    # cut trajectory according to user input
    traj_cut = universe.trajectory[frame_lims[0]:frame_lims[1]][::frame_freq]

    particle_numbers = []

    for i in traj_cut:
        particle_numbers.append(len(atom_group))

    particle_numbers = np.asarray(particle_numbers)

    return particle_numbers
# %%
def gaussian(x, A, x0, sigma): 
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
# %%