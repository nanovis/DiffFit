import numpy as np
from scipy.spatial.transform import Rotation as R
from chimerax.geometry import Place
from chimerax.core.commands import run
import os
from chimerax.atomic import AtomicStructure
import time
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
    mol_vol_matrix = mol_vol.data.matrix()
    vol_matrix = vol.data.matrix().copy()
    eligible_indices = np.where(mol_vol_matrix > mol_vol.maximum_surface_level)
    eligible_indices_list = list(zip(*eligible_indices))

    xyz_in_mol = [mol_vol.data.ijk_to_xyz(idx_in_mol_vol) for idx_in_mol_vol in eligible_indices_list]
    idx_in_vol = [vol.data.xyz_to_ijk(xyz).astype(int).tolist() for xyz in xyz_in_mol]

    for idx in idx_in_vol:
        try:
            vol_matrix[*idx] = 0.0
        except IndexError:
            # if idx_in_vol is out of bounds
            continue

    r = vol.subregion(step = [1, 1, 1])
    g = vol.region_grid(r)
    g.array[:, :, :] = vol_matrix
    g.name = vol.name + f"zero_{zero_iter}"

    from chimerax.map.volume import volume_from_grid_data
    v_zeroed = volume_from_grid_data(g, session)
    v_zeroed.copy_settings_from(vol, copy_region=False, copy_colors=True)

    mol_vol.delete()
    mol.delete()

    session.logger.info(f"Zeroing density for MQS: {MQS}")



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
    e_sqd_clusters[0] = [[mol_idx, quat_idx, shift_idx, shift 3, quat 4, corr],
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
