# %%
import numpy as np
import os 
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis as mda
# %%
def calc_average_force_through_interface_for_particles_around_x_planes(x_plane_grid, x_margin, x_interface, r_cut, epsilon, sigma, file_tpr, file_trj, format_tpr, format_trj, atom_name, frame_lims, frame_freq, file_out=None):
    average_force_vals = np.zeros((len(x_plane_grid), 3))
    
    for i, x_plane in enumerate(x_plane_grid):
        print('plane = {0:2d} of {1}'.format(i + 1, len(x_plane_grid)), flush=True, end='\r')
        average_force_vals[i] = calc_average_force_vector_through_interface_for_particles_around_x_plane(x_plane, x_margin, x_interface, r_cut, epsilon, sigma, file_tpr, file_trj, format_tpr, format_trj, atom_name, frame_lims, frame_freq)

        if file_out:
            f_out = open(file_out, "a")
            f_out.write("{0:10.4f} {1:10.4f} {2:10.4f} {3:10.4f}\n".format(x_plane, average_force_vals[i, 0], average_force_vals[i, 1], average_force_vals[i, 2]))
            f_out.close()

    return average_force_vals
# %%
def calc_average_force_vector_through_interface_for_particles_around_x_plane(x_plane, x_margin, x_interface, r_cut, epsilon, sigma, file_tpr, file_trj, format_tpr, format_trj, atom_name, frame_lims, frame_freq):
    # load files into mda universe
    universe = mda.Universe(file_tpr, file_trj, format=format_trj, topology_format=format_tpr, convert_units=False)

    # select atoms within bin
    cis_selection = 'name {0} and prop x >= {1} and prop x <= {2} and prop x <= {3}'.format(atom_name, (x_plane - x_margin), (x_plane + x_margin), x_interface)
    cis_group = universe.select_atoms(cis_selection, updating=True)

    trans_selection = 'name {0} and prop x >= {1}'.format(atom_name, x_interface)
    trans_group = universe.select_atoms(trans_selection, updating=True)

    # cut trajectory according to user input
    traj_cut = universe.trajectory[frame_lims[0]:frame_lims[1]][::frame_freq]
    
    average_forces_x_plane = np.zeros((len(traj_cut), 3))

    for i, frame in enumerate(traj_cut):
        if len(cis_group.positions) == 0 or len(trans_group.positions) == 0:
            average_forces_x_plane[i] = np.array([np.nan, np.nan, np.nan])
            continue
    
        cis_positions = cis_group.positions
        trans_positions = trans_group.positions
        average_forces_x_plane[i] = calc_average_force_vector_per_particle_in_sample_group(cis_positions, trans_positions, r_cut, epsilon, sigma, force_cap=None, capping_scheme=None) 

    average_force_x_plane = np.nanmean(average_forces_x_plane, axis=0)

    return average_force_x_plane
# %%
def calc_average_force_vector_per_particle_in_sample_group(sample_positions, reference_positions, r_cut, epsilon, sigma, force_cap, capping_scheme):
    force_vectors = calc_force_vectors_between_atom_groups(sample_positions, reference_positions, r_cut, epsilon, sigma, force_cap, capping_scheme)
    average_force_vector = np.sum(force_vectors, axis=0) / len(sample_positions)

    return average_force_vector
# %%
def calc_force_vectors_between_atom_groups(sample_positions, reference_positions, r_cut, epsilon, sigma, force_cap=None, capping_scheme='pair-wise'):    
    # calculate distances at this frame
    distances_frame = distance_array(sample_positions, reference_positions)

    # find indices of atoms that have distances smaller than cutoff
    pair_indices = np.where((distances_frame <= r_cut) & (distances_frame > 0.0))

    # compile into new distances array with only cut-off distances
    distances_cut = distances_frame[pair_indices]

    # calculate forces
    forces = calc_forces_from_distances(distances_cut, epsilon, sigma)

    # calculate distance vectors and dimensional distances between corresponding particles
    distance_vectors = sample_positions[pair_indices[0]] - reference_positions[pair_indices[1]]

    force_vectors = project_forces_onto_distance_vectors(forces, distance_vectors, distances_cut)

    if force_cap:
        # cap particle-wise or pair-wise forces
        if capping_scheme == 'pair-wise':
            pass
        elif capping_scheme == 'particle-wise':
            # add up all forces for the same inserted test particle
            force_vectors = add_force_vectors_for_each_particle(force_vectors, pair_indices[0], len(sample_positions))
            
        force_vectors = cap_force_vectors(force_vectors, force_cap)

    return force_vectors
# %%
def calc_forces_from_distances(distances, epsilon, sigma):
    forces = (24 * epsilon/distances) * (2 * (sigma/distances)**12 - (sigma/distances)**6) # lennard-jones force
    return forces
# %%
def project_forces_onto_distance_vectors(forces, distance_vectors, distances):
    distances_x, distances_y, distances_z = np.transpose(distance_vectors)

    # calculate forces between the particles 
    forces_x = project_forces_onto_directional_distances(forces, distances_x, distances)
    forces_y = project_forces_onto_directional_distances(forces, distances_y, distances)
    forces_z = project_forces_onto_directional_distances(forces, distances_z, distances)

    # format forces similar to distance vectors
    force_vectors = np.stack((forces_x, forces_y, forces_z), axis=1)

    return force_vectors
# %%
def project_forces_onto_directional_distances(forces, distances_directional, distances):
    forces = distances_directional/distances * forces
    return forces 
# %%
def add_force_vectors_for_each_particle(force_vectors, particle_indices, num_particles):
    force_vectors_per_particle = np.zeros((num_particles, 3))

    for i in range(num_particles):
        force_vectors_per_particle[i] = np.sum(force_vectors[np.where(particle_indices == i)], axis=0)
    
    return force_vectors_per_particle
# %%
def cap_force_vectors(force_vectors, force_cap):
    force_vectors[force_vectors > force_cap] = force_cap
    force_vectors[force_vectors < -force_cap] = -force_cap

    return force_vectors