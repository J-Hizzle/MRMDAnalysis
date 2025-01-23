# %%
import numpy as np
# %%
def read_data_on_grid(filename):
    '''
    Read a .txt file created by an MRMD thermodynamic force iteration run
    and return its contents in numpy arrays.

    Parameters
    ----------
    filename : str
        Name of the density file to be read.

    Returns
    -------
    pos_x_grid : np.ndarray (dtype=float)
        1D-grid on which the density is given.
    dens_vals : np.ndarray (dtype=float)
        Array containing the densities read from the file.
    '''
    file_vals = np.loadtxt(filename, dtype=float, delimiter=" ")
    pos_x_grid = file_vals[0, :] 
    data_vals = file_vals[1:, :]

    return pos_x_grid, data_vals
# %%