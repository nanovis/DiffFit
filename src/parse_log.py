import numpy as np
from scipy.spatial.transform import Rotation as R
from chimerax.geometry import Place
from chimerax.core.commands import run
import os
from chimerax.atomic import AtomicStructure
import time
from sklearn.cluster import Birch
import math
import asyncio

def shift_difference(shift1, shift2):
    """
    Calculate the Euclidean distance between two shifts.
    """
    return np.linalg.norm(shift1 - shift2)


def quaternion_angle_distance(q1, q2):
    """
    Calculate the angular distance in degrees between two quaternions.
    """
    q1 = q1 / np.linalg.norm(q1)
    q2 = q2 / np.linalg.norm(q2)
    dot = np.dot(q1, q2)
    dot = np.clip(dot, -1.0, 1.0)
    angle_rad = 2 * np.arccos(abs(dot))
    angle_deg = np.degrees(angle_rad)
    return angle_deg


def test_rl():
    print("Test RL")

#async def async_test_loop(session, count):    
#    for x in range(count):
#        session.logger.info("async test {0}".format(x))
#        await asyncio.sleep(0.5)

#async def async_test(session):  
#    task = asyncio.create_task(async_test_loop(session, 10))
#    
#    await task

def animate_MQS(e_sqd_log, mol_folder, MQS, session, clean_scene=True):

    if clean_scene:
        # delete all other structures
        structures = session.models.list(type=AtomicStructure)
        for structure in structures:
            structure.delete()

    mol_files = os.listdir(mol_folder)
    mol_path = os.path.join(mol_folder, mol_files[MQS[0]])
    mol = run(session, f"open {mol_path}")[0]

    N_iter = len(e_sqd_log[0, 0, 0])
    for iter_idx in range(N_iter):
        _, transformation = get_transformation_at_MQS(e_sqd_log, MQS, iter_idx)
        mol.scene_position = transformation
        time.sleep(0.1)

    session.logger.info(f"MQS: {MQS}")

def animate_MQS_2(e_sqd_log, mol_folder, MQS, session, clean_scene=True):

    if clean_scene:
        # delete all other structures
        structures = session.models.list(type=AtomicStructure)
        for structure in structures:
            structure.delete()

    mol_files = os.listdir(mol_folder)
    mol_path = os.path.join(mol_folder, mol_files[MQS[0]])
    mol = run(session, f"open {mol_path}")[0]

    N_iter = len(e_sqd_log[0, 0, 0])
    _, transformation = get_transformation_at_MQS(e_sqd_log, MQS, N_iter - 1)
    
    m = transformation.matrix
    
    print(m)
    
    run(session, "view matrix models #1,{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}".format(
    m[0][0],
    m[1][0],
    m[2][0],
    m[0][1],
    m[1][1],
    m[2][1],
    m[0][2],
    m[1][2],
    m[2][2],
    m[0][3],
    m[1][3],
    m[2][3]))
    #mol.scene_position = transformation
    
    session.logger.info(f"MQS: {MQS}")
    
def animate_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session, clean_scene=True):

    return


def look_at_MQS_idx(e_sqd_log, mol_folder, MQS, session, clean_scene=True):

    if clean_scene:
        # delete all other structures
        structures = session.models.list(type=AtomicStructure)
        for structure in structures:
            structure.delete()

    mol_files = os.listdir(mol_folder)

    look_at_mol_idx, transformation = get_transformation_at_MQS(e_sqd_log, MQS)

    mol_path = os.path.join(mol_folder, mol_files[look_at_mol_idx])
    mol = run(session, f"open {mol_path}")[0]

    mol.scene_position = transformation

    session.logger.info(f"MQS: {MQS}")



def look_at_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session, clean_scene=True):

    if clean_scene:
        # delete all other structures
        structures = session.models.list(type=AtomicStructure)
        for structure in structures:
            structure.delete()

    mol_files = os.listdir(mol_folder)
    # mol_files[idx] pairs with e_sqd_clusters_ordered[:][:, idx]

    look_at_mol_idx, transformation = get_transformation_at_idx(e_sqd_clusters_ordered, cluster_idx)

    mol_path = os.path.join(mol_folder, mol_files[look_at_mol_idx])
    mol = run(session, f"open {mol_path}")[0]

    mol.scene_position = transformation

    session.logger.info(f"Cluster size: {len(e_sqd_clusters_ordered[cluster_idx])}")
    session.logger.info(f"Representative MQS: {e_sqd_clusters_ordered[cluster_idx][0, 0:3].astype(int)}")
    
    return mol


def simulate_volume(session, vol, e_sqd_clusters_ordered, mol_folder, cluster_idx, res=4.0):

    mol_files = os.listdir(mol_folder)
    # mol_files[idx] pairs with e_sqd_clusters_ordered[:][:, idx]

    look_at_mol_idx, transformation = get_transformation_at_idx(e_sqd_clusters_ordered, cluster_idx)

    mol_path = os.path.join(mol_folder, mol_files[look_at_mol_idx])
    mol = run(session, f"open {mol_path}")[0]

    mol.atoms.transform(transformation)

    from chimerax.map.molmap import molecule_map
    mol_vol = molecule_map(session, mol.atoms, res, grid_spacing=vol.data_origin_and_step()[1][0] / 3)

    return mol_vol
    # TODO: Manually change the surface threshold

def zero_cluster_density(session, mol_vol, mol, vol, MQS, zero_iter=0):
    work_vol = run(session, f"volume subtract #{vol.id[0]} #{mol_vol.id[0]} scaleFactors  1.0,1000.0")
    matrix = work_vol.data.matrix()
    matrix[matrix < 0] = 0
    work_vol.update_drawings()
    work_vol.name = "working volume"

    mol_vol.delete()
    mol.delete()

    session.logger.info(f"Zeroing density for MQS: {MQS}")

    return work_vol



def get_transformation_at_MQS(e_sqd_log, MQS, iter_idx=-1):
    shift = e_sqd_log[*MQS, iter_idx][0:3]
    quat = e_sqd_log[*MQS, iter_idx][3:7][[1, 2, 3, 0]]  # convert to x,y,z,w

    R_matrix = R.from_quat(quat).as_matrix()

    T_matrix = np.zeros([3, 4])
    T_matrix[:, :3] = R_matrix
    T_matrix[:, 3] = shift

    transformation = Place(matrix=T_matrix)
    mol_idx = MQS[0]

    return mol_idx, transformation


def get_transformation_at_idx(e_sqd_clusters_ordered, look_at_idx=0):
    shift = e_sqd_clusters_ordered[look_at_idx][0, 3:6]
    quat = e_sqd_clusters_ordered[look_at_idx][0, 6:10][[1, 2, 3, 0]]  # convert to x,y,z,w

    R_matrix = R.from_quat(quat).as_matrix()

    T_matrix = np.zeros([3, 4])
    T_matrix[:, :3] = R_matrix
    T_matrix[:, 3] = shift

    transformation = Place(matrix=T_matrix)
    mol_idx = int(e_sqd_clusters_ordered[look_at_idx][0, 0])

    return mol_idx, transformation


def cluster_and_sort_sqd(e_sqd_log, shift_tolerance: float = 3.0, angle_tolerance: float = 6.0):
    """
    cluster the records in sqd table by thresholding on shift and quaternion
    then sort each cluster by their correlation in descending order
    then choose the one with the highest correlation as the representative of each cluster
    then order the clusters by their representative correlation in descending order
    return the ordered cluster

    e_sqd_clusters contents
    e_sqd_clusters[0] = [[mol_idx, quat_idx, shift_idx, shift 3, quat 4, corr 4],
                         ...
                         [mol_idx, quat_idx, shift_idx, shift 3, quat 4, corr 4]]
    with corr, e_sqd_cluster[0][:, -1] is sorted in descending order

    then e_sqd_clusters_ordered[:][0, 12] is sorted in descending order
    """
    N_mol, N_quat, N_shift, _, _ = e_sqd_log.shape

    e_sqd_clusters = []

    for mol_idx in range(N_mol):
        for quat_idx in range(N_quat):
            for shift_idx in range(N_shift):

                # Choose -1 iter
                # TODO:RL: find the iter with the largest corr
                shift = e_sqd_log[mol_idx, quat_idx, shift_idx, -1, 0:3]
                quat = e_sqd_log[mol_idx, quat_idx, shift_idx, -1, 3:7]

                hit_flag = False
                for cluster_idx in range(len(e_sqd_clusters)):
                    for placement in e_sqd_clusters[cluster_idx]:
                        if mol_idx == int(placement[0]):
                            placement_shift = placement[3:6]
                            placement_quat = placement[6:10]

                            shift_diff = shift_difference(shift, placement_shift)
                            angle_diff = quaternion_angle_distance(quat, placement_quat)

                            if shift_diff <= shift_tolerance and angle_diff <= angle_tolerance:
                                hit_flag = True
                                e_sqd_clusters[cluster_idx] = np.vstack((e_sqd_clusters[cluster_idx],
                                                                         np.array(
                                                                             [np.hstack(([mol_idx, quat_idx, shift_idx],
                                                                                         e_sqd_log[
                                                                                             mol_idx, quat_idx, shift_idx, -1]))])))
                                break

                    if hit_flag:
                        break

                if not hit_flag:
                    e_sqd_clusters.append(np.array([np.hstack(([mol_idx, quat_idx, shift_idx],
                                                               e_sqd_log[mol_idx, quat_idx, shift_idx, -1]))]))

    # e_sqd_clusters_len = [len(cluster) for cluster in e_sqd_clusters]

    # sort within each cluster by descending correlation
    e_sqd_clusters_sorted = [cluster[np.argsort(-cluster[:, 12])] for cluster in e_sqd_clusters]

    # choose the highest correlation
    e_sqd_clusters_representative = np.array([cluster[0, :] for cluster in e_sqd_clusters_sorted])

    # order all the clusters by their representatives' correlation
    clusters_order = np.argsort(-e_sqd_clusters_representative[:, 12])
    e_sqd_clusters_ordered = [e_sqd_clusters_sorted[i] for i in clusters_order]

    # e_sqd_clusters_ordered_len = [len(cluster) for cluster in e_sqd_clusters_ordered]

    return e_sqd_clusters_ordered


def q2_unit_coord(Q):

    rotations = [R.from_quat(q) for q in Q]

    up = np.array([0, 1, 0])
    right = np.array([1, 0, 0])

    rotated_up = np.array([rot.apply(up) for rot in rotations])
    rotated_right = np.array([rot.apply(right) for rot in rotations])

    return np.concatenate((rotated_up, rotated_right), axis=-1)


def cluster_and_sort_sqd_fast(e_sqd_log, shift_tolerance: float = 3.0, angle_tolerance: float = 6.0, sort_column_idx: int = 9):
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


    # Convert angle tolerance to radians and then compute the sin of half the angle
    q_coord_radius_tolerance = math.sin(math.radians(angle_tolerance / 2.0))

    N_mol, N_quat, N_shift, N_iter, N_record = e_sqd_log.shape

    e_sqd_log = e_sqd_log.reshape([N_mol, N_quat * N_shift, N_iter, N_record])

    correlations = e_sqd_log[:, :, 1:22, sort_column_idx]  # remove the 0 iteration, which is before optimization
    max_correlations_idx = np.argmax(correlations, axis=-1)

    # Generate meshgrid for the dimensions you're not indexing through
    dims_0, dims_1 = np.meshgrid(
        np.arange(e_sqd_log.shape[0]),
        np.arange(e_sqd_log.shape[1]),
        indexing='ij'
    )

    # Use the generated meshgrid and max_correlations_idx to index into e_sqd_log
    sqd_highest_corr_np = e_sqd_log[dims_0, dims_1, max_correlations_idx + 1]  # add back 0 iteration

    sqd_clusters = []
    for mol_idx in range(N_mol):
        mol_shift = sqd_highest_corr_np[mol_idx, :, :3]
        mol_q = sqd_highest_corr_np[mol_idx, :, 3:7]

        cluster_shift = Birch(threshold=shift_tolerance/2.0, n_clusters=None).fit(mol_shift)
        # mol_shift_cluster = np.concatenate([mol_shift, cluster_shift.labels_.reshape(-1, 1)], axis=-1)
        # np.save(f"mol{mol_idx}_shift_cluster.npy", mol_shift_cluster)

        mol_q_coord = q2_unit_coord(mol_q)
        cluster_q = Birch(threshold=q_coord_radius_tolerance, n_clusters=None).fit(mol_q_coord)
        # mol_q_coord_cluster = np.concatenate([mol_q_coord, cluster_q.labels_.reshape(-1, 1)], axis=-1)
        # np.save(f"mol{mol_idx}_q_coord_cluster.npy", mol_q_coord_cluster)

        mol_transforma_label = np.concatenate((cluster_shift.labels_.reshape([-1, 1]), cluster_q.labels_.reshape([-1, 1])),
                                              axis=-1)
        unique_labels, indices, counts = np.unique(mol_transforma_label, axis=0, return_inverse=True, return_counts=True)

        for cluster_idx in range(len(unique_labels)):
            sqd_idx = np.argwhere(indices == cluster_idx).reshape([-1])
            max_idx = sqd_idx[np.argsort(-sqd_highest_corr_np[mol_idx, sqd_idx, sort_column_idx])[0]]

            # [mol_idx, max_idx (in e_sqd_log), iter_idx (giving the largest correlation),
            #  cluster size, correlation]
            sqd_clusters.append([mol_idx, max_idx, max_correlations_idx[mol_idx, max_idx],
                                 counts[cluster_idx], sqd_highest_corr_np[mol_idx, max_idx, sort_column_idx]])

    sqd_clusters = np.array(sqd_clusters)
    e_sqd_clusters_ordered = sqd_clusters[np.argsort(-sqd_clusters[:, -1])]

    return e_sqd_clusters_ordered