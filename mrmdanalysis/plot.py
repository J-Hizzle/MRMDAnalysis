# %%
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'notebook', 'grid'])
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=22, labelsize=22)     # fontsize of the axes title
# %%
def plot_profile_in_adress_box(pos_grid, data_vals, r_ref, r_min, r_max, x_label, y_label, data_labels, fig_title, colors, y_lims, file_plt=None, format_plt=None):
    '''
    Plot a data profile (e.g. density or thermodynamic force) in an AdResS simulation box.
    
    Parameters
    ----------
    pos_grid : np.1darray (dtype=float)
        Grid on which the data points are given.
    data_vals : np.2darray (dtype=float, shape=(number of datasets, number of gridpoints))
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
    '''
    
    alpha_val = 0.2
    textheight = y_lims[1] - (y_lims[1] - y_lims[0])/20

    for i in range(len(data_vals)):
        plt.plot(pos_grid - r_ref, data_vals[i], label=data_labels[i], color=colors[i], linewidth=0.75)

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
        plt.savefig(fname=file_plt, dpi=150, format=format_plt)
