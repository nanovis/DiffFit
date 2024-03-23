import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import shift_difference, quaternion_angle_distance
import numpy as np

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains\e_sqd_log.npy")

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
                            e_sqd_cluster[cluster_idx] = np.vstack((e_sqd_cluster[cluster_idx],
                                                                    np.array([np.hstack(([mol_idx, quat_idx, shift_idx],
                                                                                         e_sqd_log[mol_idx, quat_idx, shift_idx, -1]))])))

            if not hit_flag:
                e_sqd_cluster.append(np.array([np.hstack(([mol_idx, quat_idx, shift_idx],
                                                          e_sqd_log[mol_idx, quat_idx, shift_idx, -1]))]))


