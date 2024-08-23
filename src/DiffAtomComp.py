import argparse
import os, sys
from datetime import datetime

import numpy as np
import torch
from typing import Tuple, List
import torch.nn.functional as F

from Bio.PDB import MMCIFParser, PDBParser
import mrcfile
import operator

import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from scipy.ndimage import label, center_of_mass

from .common import generate_random_quaternions

from scipy.spatial.transform import Rotation as R

from sklearn.cluster import Birch
import math

from math import pi
from chimerax.geometry import bins
from chimerax.geometry import Place


# Ignore PDBConstructionWarning for unrecognized 'END' record
warnings.filterwarnings("ignore", message="Ignoring unrecognized record 'END'", category=PDBConstructionWarning)


def q2_unit_coord(Q):
    rotations = [R.from_quat(q) for q in Q]

    up = np.array([0, 1, 0])
    right = np.array([1, 0, 0])

    rotated_up = np.array([rot.apply(up) for rot in rotations])
    rotated_right = np.array([rot.apply(right) for rot in rotations])

    return np.concatenate((rotated_up, rotated_right), axis=-1)


def cluster_and_sort_sqd_fast(e_sqd_log, mol_centers, shift_tolerance: float = 3.0, angle_tolerance: float = 6.0,
                              sort_column_idx: int = 7, in_contour_threshold = 0.5, correlation_threshold = 0.5):
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
    @param mol_centers: molecule atom coords centers
    @param shift_tolerance: shift tolerance in Angstrom
    @param angle_tolerance: angle tolerance in degrees
    @param sort_column_idx: the column to sort, 9-th column is the correlation
    @return: cluster representative table sorted in descending order
    """

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

    fit_res_filtered = []
    fit_res_filtered_indices = []
    in_contour_col_idx = 11
    correlation_col_idx = 9
    for mol_idx in range(N_mol):
        sqd_highest_corr_np_mol = sqd_highest_corr_np[mol_idx]

        # Fetch the columns of interest
        in_contour_percentage_column = sqd_highest_corr_np_mol[:, in_contour_col_idx]
        correlation_column = sqd_highest_corr_np_mol[:, correlation_col_idx]

        # Create masks for the filtering conditions
        in_contour_mask = in_contour_percentage_column >= in_contour_threshold
        correlation_mask = correlation_column >= correlation_threshold

        # Combine the masks to get a final filter
        combined_mask = in_contour_mask & correlation_mask

        # Apply the mask to filter the original array and also retrieve the indices
        filtered_indices = np.where(combined_mask)  # Get the indices of the filtered rows
        filtered_array = sqd_highest_corr_np_mol[filtered_indices]

        fit_res_filtered.append(filtered_array)
        fit_res_filtered_indices.append(filtered_indices[0])


    sqd_clusters = []
    for mol_idx in range(N_mol):
        mol_shift = fit_res_filtered[mol_idx][:, :3]
        mol_q = fit_res_filtered[mol_idx][:, 3:7]

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

        b = bins.Binned_Transforms(angle_tolerance * pi / 180, shift_tolerance, mol_centers[mol_idx])
        mol_transform_label = []
        unique_id = 0
        T_ID_dict = {}
        for i in range(len(mol_shift)):
            ptf = T[i]
            close = b.close_transforms(ptf)
            if len(close) == 0:
                b.add_transform(ptf)
                mol_transform_label.append(unique_id)
                T_ID_dict[id(ptf)] = unique_id
                unique_id = unique_id + 1
            else:
                mol_transform_label.append(T_ID_dict[id(close[0])])
                T_ID_dict[id(ptf)] = T_ID_dict[id(close[0])]

        unique_labels, indices, counts = np.unique(mol_transform_label, axis=0, return_inverse=True, return_counts=True)

        for cluster_idx in range(len(unique_labels)):
            sqd_idx = np.argwhere(indices == cluster_idx).reshape([-1])
            max_idx_in_filtered = sqd_idx[np.argsort(-fit_res_filtered[mol_idx][sqd_idx, sort_column_idx])[0]]
            max_idx = fit_res_filtered_indices[mol_idx][max_idx_in_filtered]

            # [mol_idx, max_idx (in e_sqd_log), iter_idx (giving the largest sort_column),
            #  cluster size, sort_metric]
            sqd_clusters.append([mol_idx, max_idx, max_sort_column_metric_idx[mol_idx, max_idx],
                                 counts[cluster_idx], fit_res_filtered[mol_idx][max_idx_in_filtered, sort_column_idx]])

    sqd_clusters = np.array(sqd_clusters)
    e_sqd_clusters_ordered = sqd_clusters[np.argsort(-sqd_clusters[:, -1])]

    return e_sqd_clusters_ordered


def get_grid3D(w, h, d, device):
    # using from torch sample
    grid_size = [w, h, d]
    return torch.from_numpy(np.indices(grid_size).reshape((len(grid_size), -1)).T).type(torch.FloatTensor).to(device)


def generate_3d_gaussian_filter(size, std, device='cpu'):
    """Generate a 3D Gaussian filter using PyTorch.

    Args:
        size (int): The size of the filter (cubic).
        std (float): The standard deviation of the Gaussian distribution.
        device (str): The device to store the tensor ('cpu' or 'cuda').

    Returns:
        torch.Tensor: A 3D Gaussian filter tensor.
    """
    # Create a grid of (x,y,z) coordinates
    axis = torch.arange(size, dtype=torch.float32, device=device) - (size - 1) / 2.0
    x, y, z = torch.meshgrid(axis, axis, axis, indexing='ij')

    # Calculate the Gaussian distribution
    gauss = torch.exp(-(x ** 2 + y ** 2 + z ** 2) / (2 * std ** 2))

    # Normalize to make the sum of all elements equal to 1
    gauss /= gauss.sum()

    return gauss


def generate_3d_laplacian_filter(size, device='cpu'):
    """Generate a 3D Laplacian filter using PyTorch.

    Args:
        size (int): The size of the filter should be odd to have a central pixel (e.g., 3).
        device (str): The device to store the tensor ('cpu' or 'cuda').

    Returns:
        torch.Tensor: A 3D Laplacian filter tensor.
    """
    # Ensure size is odd to have a center pixel
    if size % 2 == 0:
        raise ValueError("Size must be odd.")

    # Initialize the filter with zeros
    laplacian = torch.zeros((size, size, size), dtype=torch.float32, device=device)

    # Set the central value
    center = size // 2
    laplacian[center, center, center] = -6

    # Set the immediate neighbors to 1
    for d in range(1, center + 1):
        laplacian[center + d, center, center] = 1
        laplacian[center - d, center, center] = 1
        laplacian[center, center + d, center] = 1
        laplacian[center, center - d, center] = 1
        laplacian[center, center, center + d] = 1
        laplacian[center, center, center - d] = 1

    return laplacian


def compute_padding(kernel_size: Tuple[int, int, int]) -> List[int]:
    """Computes padding tuple."""
    # Adapted from kornia\filters\filter.py
    # 4 ints:  (padding_left, padding_right,padding_top,padding_bottom)
    # https://pytorch.org/docs/stable/nn.html#torch.nn.functional.pad
    assert len(kernel_size) == 3, kernel_size
    computed = [k // 2 for k in kernel_size]
    # for even kernels we need to do asymetric padding :(
    return [computed[2] - 1 if kernel_size[0] % 2 == 0 else computed[2],
            computed[2],
            computed[1] - 1 if kernel_size[1] % 2 == 0 else computed[1],
            computed[1],
            computed[0] - 1 if kernel_size[2] % 2 == 0 else computed[0],
            computed[0]]


def quaternion_to_matrix(quaternion):
    """Convert a quaternion into a 3x3 rotation matrix.
    :param quaternion: w, x, y, z
    """
    q = quaternion / torch.norm(quaternion)
    q0, q1, q2, q3 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return torch.stack([
        1 - 2 * (q2 ** 2 + q3 ** 2), 2 * (q1 * q2 - q3 * q0), 2 * (q1 * q3 + q2 * q0),
        2 * (q1 * q2 + q3 * q0), 1 - 2 * (q1 ** 2 + q3 ** 2), 2 * (q2 * q3 - q1 * q0),
        2 * (q1 * q3 - q2 * q0), 2 * (q2 * q3 + q1 * q0), 1 - 2 * (q1 ** 2 + q2 ** 2)
    ], dim=-1).view(-1, 3, 3)


def loss_between_volumes(volume1, volume2):
    """Compute the loss between two volumes
    :param volume1:
    :param volume2:
    """
    return torch.sum((volume1 - volume2) ** 2)
    # return -torch.sum(volume1 * volume2)


def numpy2tensor(np_array, device):
    target = torch.tensor(np_array, device=device).float()
    target = target.to(device)
    target_dim = target.shape
    target = target.unsqueeze(dim=0).unsqueeze(dim=0)
    return target, target_dim


def linear_norm_tensor(tensor, norm_sep: float = 0.0):
    # Normalize the tensor "target"
    # Values above the separator are normalized to [0, 1]
    # Values below the separator are normalized to [-1, 0]

    # First, create masks for the two ranges
    mask_pos = tensor > norm_sep
    mask_neg = tensor < norm_sep

    # Calculate the max and min for positive and negative values separately
    max_val_pos = tensor.max()
    min_val_pos = norm_sep

    max_val_neg = norm_sep
    min_val_neg = tensor.min()

    # Apply linear normalization for values above the separator
    tensor[mask_pos] = (tensor[mask_pos] - min_val_pos) / (max_val_pos - min_val_pos)

    # Apply linear normalization for values below the separator
    tensor[mask_neg] = (tensor[mask_neg] - min_val_neg) / (max_val_neg - min_val_neg) * -1

    # Check if any values are exactly the separator, and set them to 0
    tensor[tensor == norm_sep] = 0

    return tensor


def mrc_to_npy(mrc_filename):
    with mrcfile.open(mrc_filename, mode='r') as mrc:
        data = mrc.data
        steps = mrc.voxel_size.tolist()
        origin = mrc.header.origin.tolist()
    return data, steps, origin


def mrc_folder_to_npy_list(mrc_folder):
    sim_map_list = []

    for file_name in os.listdir(mrc_folder):
        full_path = os.path.join(mrc_folder, file_name)
        # Check if the current path is a file and not a directory
        if os.path.isfile(full_path):
            data, steps, origin = mrc_to_npy(full_path)
            sim_map_list.append((data, steps, origin))

    return sim_map_list


def normalize_coordinates_to_map_origin(coordinates, box_size, box_origin=(0.0, 0.0, 0.0)):
    # coordinates is in [x, y, z]
    # box_size is in [z, y, x]
    box_size_x_y_z = [box_size[2], box_size[1], box_size[0]]

    coordinates = coordinates - box_origin

    # Normalize coordinates to [-1, 1]
    normalized_coordinates = 2 * coordinates / box_size_x_y_z
    normalized_coordinates -= 1.0

    return normalized_coordinates


def normalize_coordinates_to_map_origin_torch(coordinates, box_size_x_y_z_tensor, box_origin_tensor):
    # Normalize coordinates to [-1, 1]
    coordinates = coordinates - box_origin_tensor

    normalized_coordinates = 2 * coordinates / box_size_x_y_z_tensor
    normalized_coordinates -= 1.0

    return normalized_coordinates


def sample_sim_map(atom_coords_list, sim_map_list, num_molecules, device):
    elements_sim_density_list = []
    for mol_idx in range(num_molecules):
        # coordinates is in [x, y, z]
        # target_size is in [z, y, x]
        sim_map, sim_map_dim = numpy2tensor(sim_map_list[mol_idx][0], device)
        sim_map_steps = sim_map_list[mol_idx][1]
        sim_map_size = np.array(list(map(operator.mul, sim_map_dim, sim_map_steps)))  # in [z, y, x]
        sim_map_size_x_y_z = [sim_map_size[2], sim_map_size[1], sim_map_size[0]]
        sim_map_size_x_y_z_tensor = torch.tensor(sim_map_size_x_y_z, device=device).float()
        sim_map_origin_tensor = torch.tensor(sim_map_list[mol_idx][2], device=device).float()

        atom_coords = torch.tensor(atom_coords_list[mol_idx], device=device).float().unsqueeze(0).unsqueeze(
            0).unsqueeze(0)
        grid = normalize_coordinates_to_map_origin_torch(atom_coords,
                                                         sim_map_size_x_y_z_tensor,
                                                         sim_map_origin_tensor)

        render = torch.nn.functional.grid_sample(sim_map, grid, 'bilinear', 'border', align_corners=True)
        element_density_flat = render.flatten()
        elements_sim_density_list.append(element_density_flat)

    return elements_sim_density_list


def transform_to_angstrom_space(ndc_shift, box_size, box_origin, atom_center_in_angstrom):
    # coordinates is in [x, y, z]
    # box_size is in [z, y, x]
    box_size_x_y_z = [box_size[2], box_size[1], box_size[0]]
    angstrom_space_shift = (ndc_shift + 1.0) * box_size_x_y_z / 2 + box_origin - atom_center_in_angstrom
    return angstrom_space_shift


def read_file_and_get_coordinates(file_path):
    # Determine file extension
    file_extension = os.path.splitext(file_path)[1].lower()

    # Initialize parser based on file extension
    if file_extension == '.cif':
        parser = MMCIFParser()
    elif file_extension == '.pdb':
        parser = PDBParser()
    else:
        raise ValueError("Unsupported file format. Please provide a .mmcif or .pdb file.")

    # Parse the structure
    structure_id = os.path.basename(file_path).split('.')[0]  # Use file name as structure ID
    structure = parser.get_structure(structure_id, file_path)

    # Initialize a list to hold all atom coordinates
    all_atom_coordinates = []

    # Iterate through each model, chain, and residue in the structure to get atoms
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # Append the atom's coordinates to the list
                    all_atom_coordinates.append(atom.get_coord())

    # Convert the list of coordinates to a numpy array
    coordinates_array = np.array(all_atom_coordinates)

    # Return the numpy array of coordinates
    return coordinates_array


def pad_and_convert_to_tensor(atom_coords_list, device):
    # Find the maximum length among all arrays
    max_length = max(len(coords) for coords in atom_coords_list)

    padded_list = []
    for coords in atom_coords_list:
        # Calculate how much padding is needed
        padding_size = max_length - len(coords)
        # Pad the current array if necessary
        if padding_size > 0:
            padded_coords = np.pad(coords, ((0, padding_size), (0, 0)), mode='constant', constant_values=0)
        else:
            padded_coords = coords
        # Convert the padded array to a torch tensor and add to the list
        padded_list.append(torch.tensor(padded_coords, device=device).float())

    # Stack all tensors into a higher-dimensional tensor
    atom_coords_tensor = torch.stack(padded_list)

    return atom_coords_tensor


def quaternion_to_matrix_batch(quaternions):
    """Convert batches of quaternions into batches of 3x3 rotation matrices for tensors of shape
    [num_quaternions, num_shifts, 4] containing quaternions (w, x, y, z)
    """
    # Normalize the quaternions
    q_norm = torch.norm(quaternions, dim=-1, keepdim=True)
    q = quaternions / q_norm
    q0, q1, q2, q3 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]

    # Calculate the components of the rotation matrices
    m00 = 1 - 2 * (q2 ** 2 + q3 ** 2)
    m01 = 2 * (q1 * q2 - q3 * q0)
    m02 = 2 * (q1 * q3 + q2 * q0)
    m10 = 2 * (q1 * q2 + q3 * q0)
    m11 = 1 - 2 * (q1 ** 2 + q3 ** 2)
    m12 = 2 * (q2 * q3 - q1 * q0)
    m20 = 2 * (q1 * q3 - q2 * q0)
    m21 = 2 * (q2 * q3 + q1 * q0)
    m22 = 1 - 2 * (q1 ** 2 + q2 ** 2)

    # Stack the components into rotation matrices
    rotation_matrices = torch.stack([m00, m01, m02, m10, m11, m12, m20, m21, m22], dim=-1)
    rotation_matrices = rotation_matrices.view(-1, quaternions.size(1), 3, 3)
    return rotation_matrices


def filter_volume(volume_np, threshold, min_island_size):
    # Step 1: Threshold the volume to create a binary mask
    binary_volume = volume_np > threshold

    # Step 2: Apply the binary mask to the original volume
    filtered_volume = np.zeros_like(volume_np)
    filtered_volume[binary_volume] = volume_np[binary_volume]

    return filtered_volume, binary_volume, None


def random_sample_indices(binary_volume, sample_size):
    # Find indices of eligible voxels
    eligible_indices = np.where(binary_volume)
    eligible_indices_list = list(zip(*eligible_indices))

    # Step 4: Sample from eligible voxels
    if len(eligible_indices_list) > sample_size:
        sampled_indices = np.random.choice(len(eligible_indices_list), size=sample_size, replace=False)
        sampled_indices = [eligible_indices_list[i] for i in sampled_indices]
    else:
        sampled_indices = eligible_indices_list  # If less than desired sample size, take all

    return sampled_indices


def transform_coords(atom_coords, e_quaternions, e_shifts, target_size_x_y_z_tensor, target_origin_tensor, device):
    atom_coords = torch.tensor(atom_coords, device=device).float()

    e_rotation_matrices = quaternion_to_matrix_batch(e_quaternions)

    transformed_coords = torch.matmul(atom_coords, e_rotation_matrices)

    transformed_coords = transformed_coords.unsqueeze(0) + e_shifts.unsqueeze(2).unsqueeze(0)

    atom_coords_normalized_to_target = normalize_coordinates_to_map_origin_torch(transformed_coords,
                                                                                 target_size_x_y_z_tensor,
                                                                                 target_origin_tensor)

    # np.save("data3D/DomainFitExample1_4.0A/atom_coords_normalized_to_target.npy", atom_coords_normalized_to_target)

    return atom_coords_normalized_to_target


def conv_volume(volume, device, conv_loops, kernel_sizes,
                negative_space_value, kernel_type="Gaussian"):
    volume_conv_list = [None] * (conv_loops + 1)
    volume_conv_list[0] = volume
    for conv_idx in range(1, conv_loops + 1):

        if kernel_type == "Gaussian":
            filter = generate_3d_gaussian_filter(kernel_sizes[conv_idx - 1], 1.0, device)
        elif kernel_type == "Laplacian":
            filter = generate_3d_laplacian_filter(kernel_sizes[conv_idx - 1], device)

        filter_padding = compute_padding(filter.shape)
        filter.unsqueeze_(0).unsqueeze_(0)

        volume_pad = F.pad(volume_conv_list[conv_idx - 1], filter_padding, mode="reflect")
        volume_conv = F.conv3d(volume_pad, filter)

        eligible_volume_tensor = volume_conv > 0.0
        volume_conv[~eligible_volume_tensor] = negative_space_value

        volume_conv_list[conv_idx] = volume_conv

    return volume_conv_list


def add_conv_density(conv_loops, conv_list, conv_weights, grid, occupied_density_sum_mol):
    if conv_loops > 0:
        for conv_idx in range(1, conv_loops + 1):
            render_conv = torch.nn.functional.grid_sample(conv_list[conv_idx], grid, 'bilinear',
                                                          'border',
                                                          align_corners=True)
            occupied_density_sum_mol += torch.sum(render_conv * conv_weights[conv_idx - 1], dim=-1).squeeze()


def read_all_files_to_atom_coords_list(structures_dir):
    atom_coords_list = []
    # List all files in the given directory
    for file_name in os.listdir(structures_dir):
        full_path = os.path.join(structures_dir, file_name)
        # Check if the current path is a file and not a directory
        if os.path.isfile(full_path):
            # Read the atom coordinates from the file
            atom_coords = read_file_and_get_coordinates(full_path)
            # Append the coordinates to the list
            atom_coords_list.append(atom_coords)
    return atom_coords_list


def rotate_centers(mol_centers, e_quaternions):
    Q = e_quaternions[:, [1, 2, 3, 0]]
    Q[:, 0:3] *= -1  # convert to scipy system, which is the same as ChimeraX and Houdini system
    # suspect that DiffFit system has something consistently negative (Right-hand vs. Left-hand probably)

    rotations = [R.from_quat(q) for q in Q]

    rotated_centers_list = []
    for atom_center in mol_centers:
        rotated_centers = np.array([rot.apply(atom_center) for rot in rotations])
        rotated_centers_list.append(rotated_centers)

    return rotated_centers_list


def calculate_metrics(render, elements_sim_density):
    # Mask to filter elements in render that are greater than zero
    mask = render > 0

    # Apply the mask to the render and elements_sim_density tensors
    render_filtered = render * mask
    elements_sim_density_filtered = elements_sim_density * mask
    mask_sum = mask.float().sum(dim=-1, keepdim=True)
    in_contour_percentage = mask.float().mean(dim=-1)

    # Calculation of correlation
    # First, normalize the inputs to have zero mean and unit variance, as Pearson's correlation requires
    render_mean = render_filtered.sum(dim=-1, keepdim=True) / mask_sum
    elements_sim_density_mean = elements_sim_density_filtered.sum(dim=-1, keepdim=True) / mask_sum

    render_std = torch.sqrt(((render_filtered - render_mean * mask) ** 2).sum(dim=-1, keepdim=True) / mask_sum)
    elements_sim_density_std = torch.sqrt(
        ((elements_sim_density_filtered - elements_sim_density_mean * mask) ** 2).sum(dim=-1, keepdim=True) / mask_sum)

    render_normalized = (render_filtered - render_mean * mask) / render_std
    elements_sim_density_normalized = (elements_sim_density_filtered - elements_sim_density_mean * mask) / elements_sim_density_std

    # Now, both tensors can be multiplied directly, and then we sum over the last dimension to compute the dot product
    # The denominator for the Pearson correlation coefficient simplifies to the product of stds, times the length of the vectors,
    # because we normalized the inputs. This results in 1 for each pair, so the dot product gives us the correlation directly.
    cam = (render_normalized * elements_sim_density_normalized).sum(dim=-1) / mask_sum.squeeze(-1)

    overlap = (render_filtered * elements_sim_density_filtered).sum(dim=-1)
    overlap_mean = overlap / mask_sum.squeeze(-1)

    render_norm = torch.linalg.vector_norm(render_filtered, dim=-1)
    elements_sim_density_norm = torch.linalg.vector_norm(elements_sim_density_filtered, dim=-1)

    correlation = overlap / (render_norm * elements_sim_density_norm)

    return torch.nan_to_num(torch.stack((overlap_mean, correlation, cam, in_contour_percentage), dim=-1))


def diff_fit(volume_list: list,
             volume_steps: list,
             volume_origin: list,
             min_island_size: int,
             mol_coords: list,
             mol_sim_maps: list,
             N_shifts: int = 10,
             N_quaternions: int = 100,
             negative_space_value: float = -0.5,
             learning_rate: float = 0.01,
             n_iters: int = 201,
             save_results: bool = False,
             out_dir: str = "DiffFit_out",
             out_dir_exist_ok: bool = False,
             device: str = "cpu"
             ):
    timer_start = datetime.now()

    if save_results:
        os.makedirs(out_dir, exist_ok=out_dir_exist_ok)
        with open(f"{out_dir}/log.log", "a") as log_file:
            log_file.write(f"Wall clock time: {datetime.now()}\n")

    # ======= load target volume to fit into
    target_no_negative = volume_list[0]
    target_steps = volume_steps
    target_origin = volume_origin
    # target steps is in [z, y, x]
    # target_origin is in [x, y, z]

    target_no_negative, eligible_volume, cluster_center_indices = filter_volume(target_no_negative, 0, min_island_size)

    sampled_indices = random_sample_indices(eligible_volume, N_shifts)
    # convert sampled_indices to angstrom space coords
    sampled_coords = np.array([np.array(idx) * np.array(target_steps) for idx in sampled_indices])
    sampled_coords = sampled_coords[:, [2, 1, 0]] + target_origin  # convert to [x, y, z] and then shift

    target_no_negative, target_dim = numpy2tensor(target_no_negative, device)
    # target as [1, 1, z, y, x]
    # target_dim as [z, y, x]

    target_size = np.array(list(map(operator.mul, target_dim, target_steps)))  # in [z, y, x]

    target_no_negative = linear_norm_tensor(target_no_negative)
    #  ### np.save(f"{os.path.dirname(target_vol_path)}/target_filtered_normalized.npy",
    #  ###         target_no_negative.squeeze().detach().cpu().numpy())

    # negative space in target volume
    eligible_volume_tensor = torch.tensor(eligible_volume, device=device).unsqueeze_(0).unsqueeze_(0)
    target = target_no_negative.clone()
    target[~eligible_volume_tensor] = negative_space_value
    # target[~eligible_volume_tensor] = -target[eligible_volume_tensor].mean()

    # ======= create convoluted target volumes
    conv_loops = len(volume_list) - 1
    conv_weights = [1.0] * conv_loops

    target_gaussian_conv_list = [numpy2tensor(vol_np, device)[0] for vol_np in volume_list]

    # ======= get atom coords
    atom_coords_list = mol_coords  # atom coords as [x, y, z]
    mol_centers = [np.mean(coords, axis=0) for coords in atom_coords_list]
    num_molecules = len(atom_coords_list)

    # read simulated map
    sim_map_list = mol_sim_maps
    elements_sim_density_list = sample_sim_map(atom_coords_list, sim_map_list, num_molecules, device)

    # ======= optimization

    # Init params

    e_quaternions = generate_random_quaternions(N_quaternions * N_shifts)

    rotated_centers_array = np.array(rotate_centers(mol_centers, e_quaternions))

    e_shifts = sampled_coords - rotated_centers_array.reshape([num_molecules, N_quaternions, N_shifts, 3])

    e_quaternions = e_quaternions.reshape([N_quaternions, N_shifts, 4])
    e_quaternions = np.repeat(e_quaternions[np.newaxis, :, :, :], num_molecules, axis=0)

    e_shifts = torch.tensor(e_shifts, device=device).float().detach().requires_grad_(True)
    e_quaternions = torch.tensor(e_quaternions, device=device).float().detach().requires_grad_(True)

    # coordinates is in [x, y, z]
    # target_size is in [z, y, x]
    target_size_x_y_z = [target_size[2], target_size[1], target_size[0]]
    target_size_x_y_z_tensor = torch.tensor(target_size_x_y_z, device=device).float()
    target_origin_tensor = torch.tensor(target_origin, device=device).float()

    # Training loop
    log_every = 10

    e_sqd_log = torch.zeros([num_molecules, N_quaternions, N_shifts, int(n_iters / 10) + 2, 12], device=device)
    # [x, y, z, w, -x, -y, -z, occupied_density_sum]

    with torch.no_grad():
        e_sqd_log[:, :, :, 0, 0:3] = e_shifts
        e_sqd_log[:, :, :, 0, 3:7] = e_quaternions

    log_idx = 0


    # Create the optimizer with different learning rates
    optimizer = torch.optim.Adam([
        {'params': [e_shifts], 'lr': target_size.mean() * 0.01},
        {'params': [e_quaternions], 'lr': learning_rate}
    ])

    for epoch in range(n_iters):
        # Forward pass

        first_layer_positive_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device)
        occupied_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device)
        metrics_table = torch.zeros([num_molecules, N_quaternions, N_shifts, 4], device=device)

        for mol_idx in range(num_molecules):
            grid = transform_coords(atom_coords_list[mol_idx],
                                    e_quaternions[mol_idx],
                                    e_shifts[mol_idx],
                                    target_size_x_y_z_tensor, target_origin_tensor, device)
            render = torch.nn.functional.grid_sample(target, grid, 'bilinear', 'border', align_corners=True)

            metrics_table[mol_idx] = calculate_metrics(render, elements_sim_density_list[mol_idx])

            occupied_density_sum[mol_idx] = torch.sum(render, dim=-1).squeeze()
            first_layer_positive_density_sum[mol_idx] = torch.sum(render * (render > 0), dim=-1).squeeze()

            add_conv_density(conv_loops, target_gaussian_conv_list, conv_weights, grid, occupied_density_sum[mol_idx])

            occupied_density_sum[mol_idx] /= len(atom_coords_list[mol_idx])

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
                e_sqd_log[:, :, :, log_idx, 0:3] = e_shifts
                e_sqd_log[:, :, :, log_idx, 3:7] = e_quaternions
                e_sqd_log[:, :, :, log_idx, 7] = first_layer_positive_density_sum
                e_sqd_log[:, :, :, log_idx, 8:12] = metrics_table

                if save_results:
                    with open(f"{out_dir}/log.log", "a") as log_file:
                        log_file.write(f"Epoch: {epoch + 1:05d}, "
                                       f"loss = {loss:.4f}\n")

    timer_stop = datetime.now()

    if save_results:
        with open(f"{out_dir}/log.log", "a") as log_file:
            log_file.write(f"Time elapsed: {timer_stop - timer_start}\n\n")

    # convert quaternion to ChimeraX, Houdini, scipy system and normalize it

    e_sqd_ChimeraX_q = torch.cat([-e_sqd_log[..., 4:7], e_sqd_log[..., 3].unsqueeze(-1)], dim=-1)
    e_sqd_log[:, :, :, :, 3:7] = e_sqd_ChimeraX_q

    q_norms = torch.linalg.vector_norm(e_sqd_log[:, :, :, :, 3:7], dim=-1, keepdim=True)
    e_sqd_log[:, :, :, :, 3:7] /= q_norms

    e_sqd_log_np = e_sqd_log.detach().cpu().numpy()

    if save_results:
        np.savez_compressed(f"{out_dir}/fit_res.npz", mol_centers=mol_centers, opt_res=e_sqd_log_np)
        # np.save(f"{out_dir}/sampled_coords.npy", sampled_coords)

    # e_sqd_log_np = e_sqd_log.detach().cpu().numpy()
    # N_mol, N_quat, N_shift, N_iter, N_metric = e_sqd_log_np
    # e_sqd_log_np = e_sqd_log_np.reshape([N_mol, N_quat * N_shift, N_iter, N_metric])
    # e_sqd_clusters_ordered = cluster_and_sort_sqd_fast(e_sqd_log_np, shift_tolerance=3.0, angle_tolerance=6.0)

    # Each record is in length of 11 as [shift 3, quat 4, quality metric 4]
    # quality metric: occupied_density_avg (idx: 7), overlap (idx: 8), correlation (idx: 9), cam (idx: 10)
    return mol_centers, e_sqd_log_np


def diff_atom_comp(target_vol_path: str,
                   target_surface_threshold: float,
                   min_cluster_size: float,
                   structures_dir: str,
                   structures_sim_map_dir: str,
                   N_shifts: int = 10,
                   N_quaternions: int = 100,
                   negative_space_value: float = -0.5,
                   learning_rate: float = 0.01,
                   n_iters: int = 201,
                   out_dir: str = "out",
                   out_dir_exist_ok: bool = False,
                   conv_loops: int = 10,
                   conv_kernel_sizes: list = (5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                   conv_weights: list = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
                   device: str = "cpu"
                   ):
    timer_start = datetime.now()
    # ======= load target volume to fit into

    target_no_negative, target_steps, target_origin = mrc_to_npy(target_vol_path)
    # target dim is in [z, y, x]
    # target_origin is in [x, y, z]

    target_no_negative, eligible_volume, _ = filter_volume(target_no_negative,
                                                                                target_surface_threshold,
                                                                                min_cluster_size)

    sampled_indices = random_sample_indices(eligible_volume, N_shifts)
    # convert sampled_indices to angstrom space coords
    sampled_coords = np.array([np.array(idx) * np.array(target_steps) for idx in sampled_indices])
    sampled_coords = sampled_coords[:, [2, 1, 0]] + target_origin  # convert to [x, y, z] and then shift

    target_no_negative, target_dim = numpy2tensor(target_no_negative, device)
    # target as [1, 1, z, y, x]
    # target_dim as [z, y, x]

    target_size = np.array(list(map(operator.mul, target_dim, target_steps)))  # in [z, y, x]

    target_no_negative = linear_norm_tensor(target_no_negative)
    # np.save(f"{os.path.dirname(target_vol_path)}/target_filtered_normalized.npy",
    #         target_no_negative.squeeze().detach().cpu().numpy())

    # negative space in target volume
    eligible_volume_tensor = torch.tensor(eligible_volume, device=device).unsqueeze_(0).unsqueeze_(0)
    target = target_no_negative.clone()
    target[~eligible_volume_tensor] = negative_space_value
    # target[~eligible_volume_tensor] = -target[eligible_volume_tensor].mean()

    # ======= create convoluted target volumes

    if len(conv_weights) != conv_loops:
        raise ValueError("Length of conv_weights does not match conv_loops! ")

    target_gaussian_conv_list = conv_volume(target_no_negative, device, conv_loops, conv_kernel_sizes,
                                            negative_space_value, kernel_type="Gaussian")

    atom_coords_list = read_all_files_to_atom_coords_list(structures_dir)  # atom coords as [x, y, z]
    mol_centers = [np.mean(coords, axis=0) for coords in atom_coords_list]
    num_molecules = len(atom_coords_list)

    # read simulated map
    sim_map_list = mrc_folder_to_npy_list(structures_sim_map_dir)
    elements_sim_density_list = sample_sim_map(atom_coords_list, sim_map_list, num_molecules, device)

    # ======= optimization

    # Init params

    e_quaternions = generate_random_quaternions(N_quaternions * N_shifts)

    rotated_centers_array = np.array(rotate_centers(mol_centers, e_quaternions))

    e_shifts = sampled_coords - rotated_centers_array.reshape([num_molecules, N_quaternions, N_shifts, 3])

    e_quaternions = e_quaternions.reshape([N_quaternions, N_shifts, 4])
    e_quaternions = np.repeat(e_quaternions[np.newaxis, :, :, :], num_molecules, axis=0)

    e_shifts = torch.tensor(e_shifts, device=device).float().detach().requires_grad_(True)
    e_quaternions = torch.tensor(e_quaternions, device=device).float().detach().requires_grad_(True)

    # coordinates is in [x, y, z]
    # target_size is in [z, y, x]
    # make the conversion here
    target_size_x_y_z = [target_size[2], target_size[1], target_size[0]]
    target_size_x_y_z_tensor = torch.tensor(target_size_x_y_z, device=device).float()
    target_origin_tensor = torch.tensor(target_origin, device=device).float()

    # Training loop
    log_every = 10

    e_sqd_log = torch.zeros([num_molecules, N_quaternions, N_shifts, int(n_iters / 10) + 2, 12], device=device)
    # [x, y, z, w, -x, -y, -z, occupied_density_sum]

    with torch.no_grad():
        e_sqd_log[:, :, :, 0, 0:3] = e_shifts
        e_sqd_log[:, :, :, 0, 3:7] = e_quaternions

    log_idx = 0
    os.makedirs(out_dir, exist_ok=out_dir_exist_ok)

    with open(f"{out_dir}/log.log", "a") as log_file:
        log_file.write(f"Wall clock time: {datetime.now()}\n")

    # Create the optimizer with different learning rates
    optimizer = torch.optim.Adam([
        {'params': [e_shifts], 'lr': target_size.mean() * 0.01},
        {'params': [e_quaternions], 'lr': learning_rate}
    ])

    for epoch in range(n_iters):
        # Forward pass

        first_layer_positive_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device)
        occupied_density_sum = torch.zeros([num_molecules, N_quaternions, N_shifts], device=device)
        metrics_table = torch.zeros([num_molecules, N_quaternions, N_shifts, 4], device=device)

        for mol_idx in range(num_molecules):
            grid = transform_coords(atom_coords_list[mol_idx],
                                    e_quaternions[mol_idx],
                                    e_shifts[mol_idx],
                                    target_size_x_y_z_tensor, target_origin_tensor, device)
            render = torch.nn.functional.grid_sample(target, grid, 'bilinear', 'border', align_corners=True)

            metrics_table[mol_idx] = calculate_metrics(render, elements_sim_density_list[mol_idx])

            occupied_density_sum[mol_idx] = torch.sum(render, dim=-1).squeeze()
            first_layer_positive_density_sum[mol_idx] = torch.sum(render * (render > 0), dim=-1).squeeze()

            add_conv_density(conv_loops, target_gaussian_conv_list, conv_weights, grid, occupied_density_sum[mol_idx])

            occupied_density_sum[mol_idx] /= len(atom_coords_list[mol_idx])

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
                e_sqd_log[:, :, :, log_idx, 0:3] = e_shifts
                e_sqd_log[:, :, :, log_idx, 3:7] = e_quaternions
                e_sqd_log[:, :, :, log_idx, 7] = first_layer_positive_density_sum
                e_sqd_log[:, :, :, log_idx, 8:12] = metrics_table

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

    np.savez_compressed(f"{out_dir}/fit_res.npz", mol_centers=mol_centers, opt_res=e_sqd_log.detach().cpu().numpy())
    # np.save(f"{out_dir}/sampled_coords.npy", sampled_coords)

    # e_sqd_log_np = e_sqd_log.detach().cpu().numpy()
    # N_mol, N_quat, N_shift, N_iter, N_metric = e_sqd_log_np
    # e_sqd_log_np = e_sqd_log_np.reshape([N_mol, N_quat * N_shift, N_iter, N_metric])
    # e_sqd_clusters_ordered = cluster_and_sort_sqd_fast(e_sqd_log_np, shift_tolerance=3.0, angle_tolerance=6.0)

    # Each record is in length of 11 as [shift 3, quat 4, quality metric 4]
    # quality metric: occupied_density_avg (idx: 7), overlap (idx: 8), correlation (idx: 9), cam (idx: 10)
    return mol_centers, e_sqd_log


def parse_floats(arg):
    try:
        return [float(x) for x in arg.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("The argument should be a comma-separated list of floats")


def parse_ints(arg):
    try:
        return [int(x) for x in arg.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("The argument should be a comma-separated list of ints")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--target_vol', type=str,
                        help="target volume in .map or .mrc format")
    parser.add_argument('--target_surface_threshold', type=float, default=0.8)
    parser.add_argument('--min_cluster_size', type=int, default=100,
                        help="The minimum number of connected voxels that can be considered as a valid cluster")

    parser.add_argument('--structures_dir', type=str,
                        help="directory containing the structures to be fit")
    parser.add_argument('--structures_sim_map_dir', type=str,
                        help="directory containing the simulated map from the structures to be fit")

    parser.add_argument('--out_dir_exist_ok', type=bool,
                        help="if True, output directory will be overwritten when existing")

    parser.add_argument('--N_shifts', type=int, default=10,
                        help="The number of random shift initializations")

    parser.add_argument('--N_quaternions', type=int, default=100,
                        help="The number of random rotation initializations")

    parser.add_argument('--negative_space_value', type=float, default=-0.5,
                        help="The value to set the negative space voxels to")

    args = parser.parse_args()

    # ======= fitting and time it
    timer_start = datetime.now()

    diff_atom_comp(args.target_vol,
                   args.target_surface_threshold,
                   args.min_cluster_size,
                   args.structures_dir,
                   args.structures_sim_map_dir,
                   out_dir_exist_ok=args.out_dir_exist_ok,
                   N_shifts=args.N_shifts,
                   N_quaternions=args.N_quaternions,
                   negative_space_value=args.negative_space_value)

    timer_stop = datetime.now()

    print("Time elapsed: ", timer_stop - timer_start)

