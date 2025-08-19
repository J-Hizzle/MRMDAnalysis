# %%
import numpy as np
from scipy import integrate
# %%
def calc_coupling_potential_from_thermodynamic_force(x_grid, force_vals, r_ref, r_min):
    pot_vals = -integrate.cumulative_simpson(force_vals, x=x_grid, initial=0)
    pot_vals -= pot_vals[np.where((x_grid >= r_ref + r_min))][0]
    return pot_vals
# %%
def calc_applied_thermodynamic_force(x_grid, force_vals, r_ref, app_min, app_max):
    r_grid = np.abs(x_grid - r_ref)
    force_vals[np.where((r_grid <= app_min) | (r_grid >= app_max))] = 0.0
    return force_vals
# %%
def calc_thermodynamic_force_from_coupling_potential(x_grid, pot_vals):
    force_vals = -np.gradient(pot_vals, x_grid)
    return force_vals
# %%
def cut_data_to_new_box(x_grid, data_vals, r_ref, app_min, app_max, r_ref_new, app_min_new, app_max_new):
    x_grid_cut = np.arange(0, r_ref_new * 2, x_grid[1] - x_grid[0])
    data_vals_cut = np.zeros_like(x_grid_cut)
    print((np.where((x_grid_cut >= r_ref_new - app_max_new) & (x_grid_cut <= r_ref_new - app_min_new))))
    print((np.where((x_grid >= r_ref - app_max) & (x_grid <= r_ref - app_min))))
    data_vals_cut[np.where((x_grid_cut >= r_ref_new - app_max_new) & (x_grid_cut <= r_ref_new - app_min_new))] = data_vals[np.where((x_grid >= r_ref - app_max) & (x_grid <= r_ref - app_min))]
    data_vals_cut[np.where((x_grid_cut >= r_ref_new + app_min_new) & (x_grid_cut <= r_ref_new + app_max_new))] = data_vals[np.where((x_grid >= r_ref + app_min) & (x_grid <= r_ref + app_max))]
    return x_grid_cut, data_vals_cut