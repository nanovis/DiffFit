import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import shift_difference, quaternion_angle_distance
import numpy as np
from chimerax.core.commands import run


e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains\e_sqd_log.npy")

e_sqd_log.shape

shift_tolerance = 3.0
angle_tolerance = 6.0

N_mol, N_quat, N_shift, _, _ = e_sqd_log.shape

# e_sqd_cluster contents
# e_sqd_cluster[0] = [[mol_idx, quat_idx, shift_idx, shift 3, quat 4, corr],
#                     ...
#                     [mol_idx, quat_idx, shift_idx, shift 3, quat 4, corr]]
# with corr, e_sqd_cluster[0][:, -1] is sorted in descending order


e_sqd_cluster = []

for mol_idx in range(N_mol):
    for quat_idx in range(N_quat):
        for shift_idx in range(N_shift):

            # Choose -1 iter
            # TODO:RL: find the iter with the largest corr
            shift = e_sqd_log[mol_idx, quat_idx, shift_idx, -1, 0:3]
            quat = e_sqd_log[mol_idx, quat_idx, shift_idx, -1, 3:7]

            hit_flag = False
            for cluster_idx in range(len(e_sqd_cluster)):
                for placement in e_sqd_cluster[cluster_idx]:
                    if mol_idx == int(placement[0]):
                        placement_shift = placement[3:6]
                        placement_quat = placement[6:10]

                        shift_diff = shift_difference(shift, placement_shift)
                        angle_diff = quaternion_angle_distance(quat, placement_quat)

                        if shift_diff <= shift_tolerance and angle_diff <= angle_tolerance:
                            hit_flag = True
                            np.vstack((e_sqd_cluster[cluster_idx],
                                       np.array([np.hstack(([mol_idx, quat_idx, shift_idx],
                                                            e_sqd_log[mol_idx, quat_idx, shift_idx, -1]))])))

            if not hit_flag:
                e_sqd_cluster.append(np.array([np.hstack(([mol_idx, quat_idx, shift_idx],
                                                          e_sqd_log[mol_idx, quat_idx, shift_idx, -1]))]))



# Extract shifts and quaternions from a small subset
shifts = e_sqd_log[:, :, :, -1, 0:3]
quaternions = e_sqd_log[:, :, :, -1, 3:7]

# Calculate shift differences for the sample
shift_diffs = np.array([[shift_difference(shifts[i], shifts[j]) for j in range(shifts.shape[0])] for i in range(shifts.shape[0])])

# Calculate quaternion angle distances for the sample
quaternion_angle_diffs = np.array([[quaternion_angle_distance(quaternions_sample[i], quaternions_sample[j]) for j in range(quaternions_sample.shape[0])] for i in range(quaternions_sample.shape[0])])

shift_diffs, quaternion_angle_diffs



e_corr = e_sqd_log[:, :, :, 21, 9]

e_corr.max()

indices_higher_than_0_1 = np.argwhere(e_corr > 0.1)

unique_values, counts = np.unique(indices_higher_than_0_1[:, 0], return_counts=True)

unique_value_counts = dict(zip(unique_values, counts))

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]


# Calculating the indices of the elements in sorted order
sorted_indices_flat = np.argsort(-e_corr, axis=None)

# Reshaping the sorted indices to match the original shape, then converting them to multi-dimensional indices
sorted_multi_indices = np.unravel_index(sorted_indices_flat, e_corr.shape)

sorted_indices = np.vstack(sorted_multi_indices).T


look_at_idx = 0

look_corr = e_corr[*sorted_indices[look_at_idx]]

shift = e_sqd_log[*sorted_indices[look_at_idx]][21, 0:3]
quat = e_sqd_log[*sorted_indices[look_at_idx]][21, 3:7][[1, 2, 3, 0]]  # convert to x,y,z,w

from scipy.spatial.transform import Rotation as R
R_matrix = R.from_quat(quat).as_matrix()

T_matrix = np.zeros([3, 4])
T_matrix[:, :3] = R_matrix
T_matrix[:, 3] = shift

from chimerax.geometry import Place
transformation = Place(matrix=T_matrix)


from chimerax.atomic import AtomicStructure
structures = session.models.list(type=AtomicStructure)
structures[0].scene_position = transformation

