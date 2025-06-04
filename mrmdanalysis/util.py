# %%
import numpy as np
from scipy import integrate
# %%
def calc_coupling_potential_from_thermodynamic_force(x_grid, force_vals, r_ref, r_min):
    pot_vals = -integrate.cumulative_trapezoid(force_vals, x_grid, initial=0)
    pot_vals -= pot_vals[np.where((x_grid >= r_ref + r_min))][0]
    return pot_vals
# %%
def calc_applied_thermodynamic_force(x_grid, force_vals, r_ref, app_min, app_max):
    r_grid = np.abs(x_grid - r_ref)
    force_vals[np.where((r_grid <= app_min) | (r_grid >= app_max))] = 0.0
    return force_vals
# %%
