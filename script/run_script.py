import numpy as np
import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import cluster_and_sort_sqd, look_at_cluster, look_at_MQS_idx, animate_MQS, animate_cluster, get_transformation_at_MQS
from chimerax.core.commands import run


from DiffAtomComp import diff_atom_comp


import os

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains_10s20q\e_sqd_log.npy")
e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)

mol_folder = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif"

cluster_idx = 0
look_at_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)

look_at_MQS_idx(e_sqd_log, mol_folder, [0,0,0], session)


animate_MQS(e_sqd_log, mol_folder, [1, 87, 2], session)
animate_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)


e_sqd_log = diff_atom_comp("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc", 0.7, 100, "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif", "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_mrc", out_dir="D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output", negative_space_value=-0.5, N_shifts=10, N_quaternions=20, out_dir_exist_ok=True)

# ======= Code for zeroing out the placed structures =======

from chimerax.map_fit.fitcmd import fitmap as F
from chimerax.atomic import AtomicStructure

structures = session.models.list(type=AtomicStructure)
s = structures[0]
s.display = False
s.display = True




import os
mol_path = os.path.join(mol_folder, os.listdir(mol_folder)[1])
mol = run(session, f"open {mol_path}")[0]

_, transformation = get_transformation_at_MQS(e_sqd_log, [1, 1, 1])

mol.atoms.transform(transformation)

from chimerax.map.molmap import molecule_map
v = molecule_map(session, mol.atoms, 4.0, grid_spacing=vol.data_origin_and_step()[1][0])
v.display = False
v.display = True
v.delete()

v.data.array.shape

v.data_origin_and_step()
vol.data_origin_and_step()

v.data.matrix()
vol.data.matrix()


v.data.size
v.data.origin
v.data.step

v.maximum_surface_level


v_matrix = v.data.matrix()
vol_matrix = vol.data.matrix()
eligible_indices = np.where(v_matrix > v.maximum_surface_level)
eligible_indices_list = list(zip(*eligible_indices))

for idx_in_v in eligible_indices_list:
    xyz = v.data.ijk_to_xyz(idx_in_v)
    idx_in_vol = vol.data.xyz_to_ijk(xyz).astype(int)
    try:
        vol_matrix[idx_in_vol] = 0.0
    except IndexError:
        # if idx_in_vol is out of bounds
        continue

v_matrix.shape
type(v_matrix)

_model_set_position

# v.replace_data()



# playground draft code ========

