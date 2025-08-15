# %%
import numpy as np
import os
from scipy.optimize import curve_fit
import matplotlib as mpl
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title

from mrmdanalysis.statistics import gaussian
# %%
def plot_profile_in_adress_box(pos_grids, data_valss, r_ref, r_min, r_max, x_label, y_label, data_labels, fig_title, colors, y_lims, file_plt=None, format_plt=None, linestyles=None):
    '''
    Plot a data profile (e.g. density or thermodynamic force) in an AdResS simulation box.
    
    Parameters
    ----------
    pos_grids : np.2darray (dtype=float, shape=(number of datasets, number of gridpoints))
        Grids on which the data points are given.
    data_valss : np.2darray (dtype=float, shape=(number of datasets, number of gridpoints))
        Array of data values to be plotted on the grid. 
    r_ref : float
        Center point of the AT region.
    r_min : float
        Distance from AT center to the AT/Delta interface.
    r_max : float 
        Distance from AT center to the Delta/TR interface. 
    x_label : str
        Label for the x-axis.
    y_label : str
        Label for the y-axis.
    data_labels : np.1darray (dtype=str)
        Labels for the datasets.
    fig_title : str
        Title of the figure.
    colors : np.1darray (dtype=str)
        Colors for the datasets.
    y_lims : list (dtype=float)
        Limits for the y-axis.
    file_plt : str
        Name of the file to print plot to.
    format_plt : str
        Format in which the plot is supposed to be saved. 
    linestyles : list (dtype=str)
        List of linestyles for the datasets. If None, all datasets will be plotted with a solid line.
    '''
    
    if linestyles is None:
        linestyles = ['solid'] * len(data_valss)

    alpha_val = 0.2
    textheight = y_lims[1] - (y_lims[1] - y_lims[0])/20

    for i in range(len(data_valss)):
        plt.plot(pos_grids[i] - r_ref, data_valss[i], label=data_labels[i], color=colors[i], linewidth=0.75, linestyle=linestyles[i])

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.title(fig_title)
    plt.legend(framealpha=0.0, fontsize=12)

    plt.axvline(-r_min, color='black', alpha=2*alpha_val)
    plt.axvline(r_min, color='black', alpha=2*alpha_val)
    plt.axvline(-r_max, color='black', alpha=2*alpha_val)
    plt.axvline(r_max, color='black', alpha=2*alpha_val)

    plt.text(0.0, textheight, r'AT', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(r_min + (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(-r_min - (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(r_max + r_min, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(-r_max - r_min, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.ylim(y_lims[0], y_lims[1])

    if file_plt:
        if os.path.isfile(file_plt):
            print('file already exists! Aborting.')
        else:
            plt.savefig(fname=file_plt, dpi=150, format=format_plt)
# %%
def plot_profile_in_half_adress_box(pos_grids, data_valss, r_ref, r_min, r_max, x_label, y_label, data_labels, fig_title, colors, x_lims, y_lims, file_plt=None, format_plt=None, loc='best', linestyles=None):
    if linestyles is None:
        linestyles = ['solid'] * len(data_valss)

    alpha_val = 0.2
    textheight = y_lims[1] - (y_lims[1] - y_lims[0])/20

    for i in range(len(data_valss)):
        r_grid = np.array_split(pos_grids[i] - r_ref, 2)[-1]
        data_vals = np.array_split(data_valss[i], 2)[-1]
        plt.plot(r_grid, data_vals, label=data_labels[i], color=colors[i], linewidth=0.75, linestyle=linestyles[i])

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    plt.title(fig_title)
    plt.legend(framealpha=0.0, fontsize=12, loc=loc)

    plt.axvline(r_min, color='black', alpha=2*alpha_val)
    plt.axvline(r_max, color='black', alpha=2*alpha_val)

    plt.text(x_lims[0] + (r_min - x_lims[0])/2, textheight, r'AT', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(r_min + (r_max - r_min)/2, textheight, r'$\Delta$', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.text(r_max + (x_lims[1] - r_max)/2, textheight, r'TR', fontsize=18, horizontalalignment='center', verticalalignment='center')
    plt.xlim(x_lims[0], x_lims[1])
    plt.ylim(y_lims[0], y_lims[1])

    if file_plt:
        if os.path.isfile(file_plt):
            print('file already exists! Aborting.')
        else:
            plt.savefig(fname=file_plt, dpi=150, format=format_plt)
# %%
def plot_PofNs(N_grid, P_valss, labels, colors, markers, fig_title, file_out=None, **kwargs):
    '''
    Calculate and plot the probability density of finding a certain number of particles within an open system 
    in an AdResS or full-atomistic MD simulation done in Gromacs.

    Parameters:
    -----------
        N_grid : np.ndarray (dtype : float)
            Grid of particle numbers for which to plot the probability density. 
        P_valss : list (dtype : np.ndarray (dtype : int))
            Histograms of particle number frequencies.
        **kwargs : Keyword arguments
            Keyword arguments directly passed on to plt.scatter.
    '''
    for i, P_vals in enumerate(P_valss):
        label = labels[i]
        color = colors[i]
        marker = markers[i]

        plt.scatter(N_grid[:-1], P_vals, label=label, color=color, marker=marker, **kwargs)
    
    plt.legend(framealpha=0)
    plt.title(fig_title)
    plt.xlabel(r'$N$ in $1$')
    plt.ylabel(r'$P(N)$ in $1$')

    if file_out:
        plt.savefig(fname=file_out, format='png', dpi=150)
    
    plt.show()
# %%
def plot_PofNs_fit(N_grid, P_valss, initial_guesses, fig_title, xlims, labels, colors, markers, file_out=None, **kwargs):
    '''
    Calculate and plot Gaussian-fitted P(N)'s shifted to zero.

    Parameters:
    -----------
        N_grid : np.ndarray (dtype : float)
            Grid of particle numbers for which to plot the probability density. 
        P_valss : list (dtype : np.ndarray (dtype : int))
            Histograms of particle number frequencies.
        initial_guesses : np.ndarray()
        **kwargs : Keyword arguments
            Keyword arguments directly passed on to plt.scatter.
    '''
    for i, P_vals in enumerate(P_valss):
        initial_guess = initial_guesses[i]
        label = labels[i]
        color = colors[i]
        marker = markers[i]

        # final iteration of tf
        fit_params, _ = curve_fit(gaussian, N_grid[:-1], P_vals, p0=initial_guess) 
        fit_function = gaussian(N_grid[:-1], fit_params[0], fit_params[1], fit_params[2])

        # shift values
        N_max = N_grid[np.where(fit_function == np.max(fit_function))[0][0]]
        N_grid_shifted = N_grid[:-1] - N_max

        plt.scatter(N_grid_shifted, P_vals, color=color, marker=marker, alpha=0.2, **kwargs)
        plt.plot(N_grid_shifted, fit_function, label=label, color=color)

    plt.legend(framealpha=0)
    plt.title(fig_title)
    plt.xlim(xlims[0], xlims[1])
    plt.xlabel(r'$N$ - $\langle N \rangle$ in $1$')
    plt.ylabel(r'$P(N)$ in $1$')

    if file_out:
        plt.savefig(fname=file_out, format='png', dpi=150)
    
    plt.show()
# %%
def plot_Noft(time_grids, particle_numberss, labels, colors, markers, fig_title, file_out=None, run_mod=None):
    '''
    Plot the number of particles over the time of a simulation.

    Parameters: 
    -----------
        time_grid : np.ndarray (dtype : float)
            Grid displaying the simulation time.
        particle_numberss : list (dtype : np.ndarray (dtype : float))
            list of particle number arrays over time for different simulations.
        run_mod : int
            calculate and plot the running mean of every 'run_mod' elements.
    '''
    for i, time_grid in enumerate(time_grids):
        particle_numbers = particle_numberss[i]
        label = labels[i]
        color = colors[i]
        marker = markers[i]

        # do averaging of each 'run_mod'-steps
        running_mean = np.mean(particle_numbers.reshape(-1, run_mod), axis=1)
        time_grid_avrg = time_grid[::run_mod][:] + run_mod * (time_grid[1] - time_grid[0]) / 2 # plot within half the interval
        
        plt.plot(time_grid, particle_numbers, alpha=0.2, color=color)   
        plt.plot(time_grid_avrg, running_mean, label=label, color=color, marker=marker)

    plt.xlabel(r'$t$ in ${\tau}$')
    plt.ylabel(r'$N$ in $1$')
    plt.legend(framealpha=0)
    plt.title(fig_title)
    
    if file_out: 
        plt.savefig(fname=file_out, format='png', dpi=150)
    plt.show()
# %%
def select_color_from_colormap(range_indicator, colormap):
    '''
    indicator : float
        Value between 0 and 1 that indicates where in the range of the colormap the chosen color is supposed to be.
    colormap : str
        Name of the colormap to be picked from.
    '''
    cmap = mpl.colormaps[colormap]

    color = cmap(range_indicator)

    return color