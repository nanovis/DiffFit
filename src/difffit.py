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

    return target_gaussian_conv_list, target, target_no_negative, target_size_x_y_z_tensor, target_origin_tensor, sampled_coords

def prepare_atoms(structure_path, fit_atom_mode):
    """
    Read and prepare atom coordinates from files.
    """
    atom_coords_list = [read_file_and_get_coordinates(structure_path, fit_atom_mode)]
    mol_centers = [np.mean(coords, axis=0) for coords in atom_coords_list]
    atom_coords_list = center_atom_coords_list(atom_coords_list, mol_centers)

    mol_num_atoms = [len(coords) for coords in atom_coords_list]

    return atom_coords_list, mol_centers, mol_num_atoms


def optimize_fitting(target,
                     target_gaussian_conv_list,
                     atom_coords_list,
                     sampled_coords,
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

    e_quaternions = generate_random_quaternions(N_quaternions * N_shifts).reshape([N_quaternions, N_shifts, 4])
    e_quaternions = np.repeat(e_quaternions[np.newaxis, :, :, :], num_molecules, axis=0)

    e_shifts = np.tile(
        sampled_coords.reshape(1, 1, N_shifts, 1, 3),
        (num_molecules, N_quaternions, 1, 1, 1)
    )

    e_shifts = torch.tensor(e_shifts, device=device, dtype=precision).detach().requires_grad_(True)
    e_quaternions = torch.tensor(e_quaternions, device=device, dtype=precision).detach().requires_grad_(True)


    # Training loop
    log_every = 10

    e_sqd_log = torch.zeros([num_molecules, N_quaternions, N_shifts, int(n_iters / 10) + 2, 9], device=device,
                            dtype=precision)
    # [x, y, z, w, -x, -y, -z, occupied_density_sum]

    with torch.no_grad():
        e_sqd_log[:, :, :, 0, 0:3] = e_shifts.squeeze(-2)
        e_sqd_log[:, :, :, 0, 3:7] = e_quaternions

    log_idx = 0
    os.makedirs(out_dir, exist_ok=out_dir_exist_ok)

    # Create the optimizer with different learning rates
    optimizer = torch.optim.Adam([
        {'params': [e_shifts], 'lr': target_size_x_y_z_tensor.mean() * learning_rate},
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


def save_results(out_dir, target_vol_path, target_surface_threshold, structures_dir, mol_num_atoms, e_sqd_log):
    """
    Save the results to the output directory.
    """
    mol_paths = [str(Path(structures_dir) / file) for file in sorted(os.listdir(structures_dir)) if file.endswith(('.pdb', '.cif'))]

    os.makedirs(out_dir, exist_ok=True)
    np.savez_compressed(f"{out_dir}/fit_res.npz",
                        target_vol_path=target_vol_path,
                        target_surface_threshold=target_surface_threshold,
                        mol_paths=mol_paths,
                        mol_num_atoms=mol_num_atoms,
                        opt_res=e_sqd_log.detach().cpu().numpy())


def difffit(target_vol_path,
            target_surface_threshold,
            structures_dir,
            fit_atom_mode="Backbone",
            Gaussian_mode="Gaussian with negative (shrink)",
            N_shifts=10,
            N_quaternions=100,
            negative_space_value=-0.5,
            learning_rate=0.01,
            n_iters=101,
            out_dir="out",
            out_dir_exist_ok=True,
            conv_loops=3,
            conv_kernel_sizes=(5, 5, 5),
            conv_weights=(1.0, 1.0, 1.0),
            device="cuda",
            precision=torch.float32):
    """
    Main function to fit structures to the target volume.
    """
    timer_start = datetime.now()

    target, target_no_negative, target_dim, target_size, target_origin, sampled_coords = process_volume(
        target_vol_path, target_surface_threshold, device, precision)

    atom_coords_list, mol_centers, mol_num_atoms = prepare_atoms(structures_dir, fit_atom_mode)

    e_sqd_log = optimize_fitting(target,
                                 atom_coords_list, sampled_coords, target_size, target_origin,
                                 conv_loops, conv_kernel_sizes, conv_weights,
                                 len(atom_coords_list), n_iters, learning_rate,
                                 device, precision)

    save_results(out_dir, target_vol_path, target_surface_threshold, structures_dir, mol_num_atoms, e_sqd_log)

    print(f"Time elapsed: {datetime.now() - timer_start}")
    return target_vol_path, target_surface_threshold, structures_dir, mol_num_atoms, e_sqd_log
