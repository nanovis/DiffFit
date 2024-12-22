import os
import numpy as np
import torch
from datetime import datetime
from pathlib import Path
import operator
from DiffAtomComp import (mrc_to_npy,
                          filter_volume,
                          random_sample_indices,
                          numpy2tensor,
                          linear_norm_tensor,
                          conv_volume,
                          read_file_and_get_coordinates,
                          center_atom_coords_list,
                          generate_random_quaternions,
                          transform_coords,
                          add_conv_density)

from scipy.spatial.transform import Rotation as R
from math import pi


def parse_precision(precision_str):
    """
    Parse precision string into torch.dtype.
    """
    precision_map = {
        "float32": torch.float32,
        "float64": torch.float64,
        "float16": torch.float16
    }
    return precision_map.get(precision_str.lower(), torch.float32)


def process_volume(target_vol_path,
                   target_surface_threshold,
                   N_shifts,
                   negative_space_value,
                   conv_loops,
                   conv_kernel_sizes,
                   conv_weights,
                   Gaussian_mode,
                   device,
                   precision):
    """
    Load and process the target volume.
    """
    target_no_negative, target_steps, target_origin = mrc_to_npy(target_vol_path)
    target_no_negative, eligible_volume, _ = filter_volume(target_no_negative, target_surface_threshold)

    sampled_indices = random_sample_indices(eligible_volume, N_shifts)
    sampled_coords = np.array([np.array(idx) * np.array(target_steps) for idx in sampled_indices])
    sampled_coords = sampled_coords[:, [2, 1, 0]] + target_origin  # Convert to [x, y, z] and shift

    target_no_negative, target_dim = numpy2tensor(target_no_negative, device, precision)
    # target as [1, 1, z, y, x]
    # target_dim as [z, y, x]
    target_size = np.array(list(map(operator.mul, target_dim, target_steps)))  # in [z, y, x]
    # coordinates is in [x, y, z]
    # target_size is in [z, y, x]
    target_size_x_y_z = [target_size[2], target_size[1], target_size[0]]
    target_size_x_y_z_tensor = torch.tensor(target_size_x_y_z, device=device, dtype=precision)
    target_origin_tensor = torch.tensor(target_origin, device=device, dtype=precision)

    target_no_negative = linear_norm_tensor(target_no_negative)
    # negative space in target volume
    eligible_volume_tensor = torch.tensor(eligible_volume, device=device, dtype=torch.bool).unsqueeze_(0).unsqueeze_(0)
    target = target_no_negative.clone()
    target[~eligible_volume_tensor] = negative_space_value  # Placeholder for negative space value

    # ======= create convoluted target volumes

    if len(conv_weights) != conv_loops:
        raise ValueError("Length of conv_weights does not match conv_loops! ")

    target_gaussian_conv_list = conv_volume(target_no_negative, device, conv_loops, conv_kernel_sizes,
                                            negative_space_value, kernel_type="Gaussian", mode=Gaussian_mode)

    return target_gaussian_conv_list, target, target_no_negative, target_size.mean(), target_size_x_y_z_tensor, target_origin_tensor, sampled_coords


def prepare_atoms(structure_path, fit_atom_mode):
    """
    Read and prepare atom coordinates from files.
    """
    atom_coords_list = [read_file_and_get_coordinates(structure_path, fit_atom_mode)]
    mol_centers = [np.mean(coords, axis=0) for coords in atom_coords_list]
    atom_coords_list = center_atom_coords_list(atom_coords_list, mol_centers)

    mol_num_atoms = [len(coords) for coords in atom_coords_list]

    return atom_coords_list, mol_centers, mol_num_atoms


def initialize_tensors(N_quaternions, N_shifts, num_molecules, sampled_coords, n_iters, device, precision):
    """
    Initialize quaternion and shift tensors for optimization.

    :param N_quaternions: Number of quaternion rotations.
    :param N_shifts: Number of sampled shifts.
    :param num_molecules: Number of molecules to process.
    :param sampled_coords: Coordinates sampled for shifts.
    :param n_iters: The number of iterations.
    :param device: Device for tensor computation (e.g., 'cuda' or 'cpu').
    :param precision: Torch precision (e.g., torch.float32).
    :return: Initialized quaternion and shift tensors as torch Tensors.
    """
    e_quaternions = generate_random_quaternions(N_quaternions * N_shifts).reshape([N_quaternions, N_shifts, 4])
    e_quaternions = np.repeat(e_quaternions[np.newaxis, :, :, :], num_molecules, axis=0)
    e_quaternions = torch.tensor(e_quaternions, device=device, dtype=precision)

    e_shifts = np.tile(
        sampled_coords.reshape(1, 1, N_shifts, 1, 3),
        (num_molecules, N_quaternions, 1, 1, 1)
    )
    e_shifts = torch.tensor(e_shifts, device=device, dtype=precision)

    e_sqd_log = torch.zeros([num_molecules, N_quaternions, N_shifts, int(n_iters / 10) + 2, 9], device=device,
                            dtype=precision)
    # [x, y, z, w, -x, -y, -z, occupied_density_sum]

    return e_quaternions, e_shifts, e_sqd_log


def optimize_fitting(target,
                     target_gaussian_conv_list,
                     atom_coords_list,
                     sampled_coords,
                     e_quaternions,
                     e_shifts,
                     e_sqd_log,
                     target_size_mean,
                     target_size_x_y_z_tensor,
                     target_origin_tensor,
                     N_quaternions,
                     N_shifts,
                     conv_loops,
                     conv_weights,
                     num_molecules,
                     n_iters,
                     learning_rate,
                     device,
                     precision,
                     out_dir,
                     out_dir_exist_ok=True):
    """
    Perform optimization for fitting structures into the target volume.
    """
    timer_start = datetime.now()

    e_shifts = e_shifts.detach().clone().requires_grad_(True)
    e_quaternions = e_quaternions.detach().clone().requires_grad_(True)

    # Training loop
    log_every = 10

    with torch.no_grad():
        e_sqd_log[:, :, :, 0, 0:3] = e_shifts.squeeze(-2)
        e_sqd_log[:, :, :, 0, 3:7] = e_quaternions

    log_idx = 0
    os.makedirs(out_dir, exist_ok=out_dir_exist_ok)

    # Create the optimizer with different learning rates
    optimizer = torch.optim.Adam([
        {'params': [e_shifts], 'lr': target_size_mean * learning_rate},
        {'params': [e_quaternions], 'lr': learning_rate}
    ])

    atom_coords_torch_list = [torch.tensor(atom_coords, device=device, dtype=precision) for atom_coords in
                              atom_coords_list]

    for epoch in range(n_iters):
        # Forward pass

        first_layer_positive_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device,
                                                       dtype=precision)
        in_contour_percentage = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device, dtype=precision)
        occupied_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device, dtype=precision)

        for mol_idx in range(num_molecules):
            # sampled_coords = atom_coords_torch_list[mol_idx][torch.randint(0, atom_coords_torch_list[mol_idx].shape[0], (500,), device=device)]
            grid = transform_coords(atom_coords_torch_list[mol_idx],
                                    e_quaternions[mol_idx:mol_idx + 1],
                                    e_shifts[mol_idx:mol_idx + 1],
                                    target_size_x_y_z_tensor, target_origin_tensor, device)
            render = torch.nn.functional.grid_sample(target, grid, 'bilinear', 'border', align_corners=True)
            occupied_density_sum[mol_idx] = torch.sum(render, dim=-1).squeeze()
            add_conv_density(conv_loops, target_gaussian_conv_list, conv_weights, grid, occupied_density_sum[mol_idx])
            occupied_density_sum[mol_idx] /= len(atom_coords_list[mol_idx])

            with torch.no_grad():
                positive_mask = render > 0
                in_contour_percentage[mol_idx] = positive_mask.to(precision).mean(dim=-1)
                first_layer_positive_density_sum[mol_idx] = torch.sum(render * positive_mask, dim=-1).squeeze()

        # loss
        loss = -torch.sum(occupied_density_sum)
        # gradients
        loss.backward()

        # update weights
        optimizer.step()
        optimizer.zero_grad()

        # log
        if (epoch - 1) % log_every == (log_every - 1):
            with torch.no_grad():
                log_idx += 1
                e_sqd_log[:, :, :, log_idx, 0:3] = e_shifts.squeeze(-2)
                e_sqd_log[:, :, :, log_idx, 3:7] = e_quaternions
                e_sqd_log[:, :, :, log_idx, 7] = first_layer_positive_density_sum
                e_sqd_log[:, :, :, log_idx, 8] = in_contour_percentage

                with open(f"{out_dir}/log.log", "a") as log_file:
                    log_file.write(f"Epoch: {epoch + 1:05d}, "
                                   f"loss = {loss:.4f}\n")

    timer_stop = datetime.now()

    with open(f"{out_dir}/log.log", "a") as log_file:
        log_file.write(f"Time elapsed: {timer_stop - timer_start}\n\n")

    # convert quaternion to ChimeraX, Houdini, scipy system and normalize it

    e_sqd_ChimeraX_q = torch.cat([-e_sqd_log[..., 4:7], e_sqd_log[..., 3].unsqueeze(-1)], dim=-1)
    e_sqd_log[:, :, :, :, 3:7] = e_sqd_ChimeraX_q

    q_norms = torch.linalg.vector_norm(e_sqd_log[:, :, :, :, 3:7], dim=-1, keepdim=True)
    e_sqd_log[:, :, :, :, 3:7] /= q_norms

    return e_sqd_log


def cluster_and_sort_sqd_fast(e_sqd_log, shift_tolerance: float = 3.0, angle_tolerance: float = 6.0,
                              sort_column_idx: int = 7,
                              in_contour_threshold: float = 0.5,
                              save_log=False,
                              log_path="",
                              max_fits=10000,
                              max_clusters=100):
    """
    Cluster the fitting results in sqd table by thresholding on shift and quaternion
    Return the sorted cluster representatives

    How it works briefly:
    1. From all iterations, get the iteration with the highest correlation, or the metric at the sort_column_idx
    2. For each molecule:
        2.1. cluster the shift using half shift_tolerance as radius in Birch clustering algorithm
        2.2. convert the quaternion by applying to [0, 1, 0] and [1, 0, 0] to form 6 dim coords
        2.3. cluster the 6 dim coords using half secant calculated from half angle_tolerance as radius in Birch clustering algorithm
        2.4. combine two clusters to form unique clusters
        2.5. select a representative from each cluster as the one with the highest correlation, or the metric at the sort_column_idx
        2.5. record [mol_idx, max_idx, iter_idx, cluster size, correlation] for each cluster's representative
    3. sort the cluster table in descending order by correlation, or the metric at the sort_column_idx

    @param e_sqd_log: fitting results in sqd table
    @param shift_tolerance: shift tolerance in Angstrom
    @param angle_tolerance: angle tolerance in degrees
    @param sort_column_idx: the column to sort, 9-th column is the correlation
    @return: cluster representative table sorted in descending order
    """
    from chimerax.geometry import Place
    from DiffFit_bins import DiffFit_Binned_Transforms

    N_mol, N_record, N_iter, N_metric = e_sqd_log.shape

    sort_column_metric = e_sqd_log[:, :, 1:22, sort_column_idx]  # remove the 0 iteration, which is before optimization
    max_sort_column_metric_idx = np.argmax(sort_column_metric, axis=-1) + 1  # add back 0 iteration

    # Generate meshgrid for the dimensions you're not indexing through
    dims_0, dims_1 = np.meshgrid(
        np.arange(e_sqd_log.shape[0]),
        np.arange(e_sqd_log.shape[1]),
        indexing='ij'
    )

    # Use the generated meshgrid and max_sort_column_metric_idx to index into e_sqd_log
    sqd_highest_corr_np = e_sqd_log[dims_0, dims_1, max_sort_column_metric_idx]

    timer_start = datetime.now()
    if save_log:
        with open(log_path, "a") as log_file:
            log_file.write(f"DiffFit fit_res filtering starts: {timer_start}\n")

    fit_res_filtered = []
    fit_res_filtered_indices = []
    in_contour_col_idx = 8

    for mol_idx in range(N_mol):
        sqd_highest_corr_np_mol = sqd_highest_corr_np[mol_idx]

        # Fetch the columns of interest
        in_contour_percentage_column = sqd_highest_corr_np_mol[:, in_contour_col_idx]

        # Create masks for the filtering conditions
        in_contour_mask = in_contour_percentage_column >= in_contour_threshold

        # Apply the mask to filter the original array and also retrieve the indices
        filtered_indices = np.where(in_contour_mask)  # Get the indices of the filtered rows
        filtered_array = sqd_highest_corr_np_mol[filtered_indices]

        sorted_indices = np.argsort(filtered_array[:, sort_column_idx])[::-1]
        top_indices = sorted_indices[:min(max_fits, len(filtered_array))]
        filtered_array_top = filtered_array[top_indices]
        filtered_indices_top = filtered_indices[0][top_indices]

        fit_res_filtered.append(filtered_array_top)
        fit_res_filtered_indices.append(filtered_indices_top)

    if save_log:
        with open(log_path, "a") as log_file:
            log_file.write(f"DiffFit fit_res filtering time elapsed: {datetime.now() - timer_start}\n"
                           f"-------\n")

    sqd_clusters = []
    for mol_idx in range(N_mol):
        sqd_clusters_mol = []
        mol_shift = fit_res_filtered[mol_idx][:, :3]
        mol_q = fit_res_filtered[mol_idx][:, 3:7]

        if save_log:
            with open(log_path, "a") as log_file:
                log_file.write(f"Clustering {len(mol_shift)} fits for mol_idx: {mol_idx}\n")
        timer_start = datetime.now()

        T = []
        for i in range(len(mol_shift)):
            shift = mol_shift[i]
            quat = mol_q[i]
            R_matrix = R.from_quat(quat).as_matrix()

            T_matrix = np.zeros([3, 4])
            T_matrix[:, :3] = R_matrix
            T_matrix[:, 3] = shift

            transformation = Place(matrix=T_matrix)
            T.append(transformation)

        if save_log:
            with open(log_path, "a") as log_file:
                log_file.write(f"Convert to matrix time: {datetime.now() - timer_start}\n")
        timer_start = datetime.now()

        b = DiffFit_Binned_Transforms(angle_tolerance * pi / 180, shift_tolerance)
        mol_transform_label = []
        unique_id = 0
        T_ID_dict = {}
        for i in range(len(mol_shift)):
            ptf = T[i]
            in_cluster = b.any_close_transform(ptf)
            if in_cluster is None:
                b.add_transform(ptf)
                mol_transform_label.append(unique_id)
                T_ID_dict[id(ptf)] = unique_id
                unique_id = unique_id + 1
            else:
                mol_transform_label.append(T_ID_dict[id(in_cluster)])
                T_ID_dict[id(ptf)] = T_ID_dict[id(in_cluster)]

            if save_log and (i + 1) % 10000 == 0:
                with open(log_path, "a") as log_file:
                    log_file.write(f"Clustered {i+1} fits: {datetime.now()}\n")

        if save_log:
            with open(log_path, "a") as log_file:
                log_file.write(f"ChimeraX bin clustering: {datetime.now() - timer_start}\n")

        unique_labels, indices, counts = np.unique(mol_transform_label, axis=0, return_inverse=True, return_counts=True)

        for cluster_idx in range(len(unique_labels)):
            sqd_idx = np.argwhere(indices == cluster_idx).reshape([-1])
            max_idx_in_filtered = sqd_idx[np.argsort(-fit_res_filtered[mol_idx][sqd_idx, sort_column_idx])[0]]
            max_idx = fit_res_filtered_indices[mol_idx][max_idx_in_filtered]

            # [mol_idx, max_idx (in e_sqd_log), iter_idx (giving the largest sort_column),
            #  cluster size, sort_metric]
            sqd_clusters_mol.append([mol_idx, max_idx, max_sort_column_metric_idx[mol_idx, max_idx],
                                     counts[cluster_idx],
                                     fit_res_filtered[mol_idx][max_idx_in_filtered, sort_column_idx]])

        # ======= Filter cluster by the density, keep max_clusters=100 clusters for each mol
        if len(sqd_clusters_mol) > 0:
            sqd_clusters_mol = np.array(sqd_clusters_mol)
            sqd_clusters_mol = sqd_clusters_mol[np.argsort(-sqd_clusters_mol[:, -1])]
            sqd_clusters_mol = sqd_clusters_mol[:min(max_clusters, len(sqd_clusters_mol)), :]

            sqd_clusters.append(sqd_clusters_mol)

    sqd_clusters = np.vstack(sqd_clusters)

    if len(sqd_clusters) == 0:
        return None

    # e_sqd_clusters_ordered = sqd_clusters[np.argsort(-sqd_clusters[:, -1])]

    return sqd_clusters