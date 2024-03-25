import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\src')
from parse_log import cluster_and_sort_sqd, look_at_cluster, look_at_MQS_idx, animate_MQS, animate_cluster
import numpy as np
from chimerax.core.commands import run
import os

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains_10s20q\e_sqd_log.npy")
e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)

mol_folder = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif"

cluster_idx = 0
look_at_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)  # read from log to get MQS

look_at_MQS_idx(e_sqd_log, mol_folder, [0,0,0], session)


animate_MQS(e_sqd_log, mol_folder, [1, 87, 2], session)
animate_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)



import os
mol_path = os.path.join(mol_folder, os.listdir(mol_folder)[1])  # MQS [0]
mol = run(session, f"open {mol_path}")[0]

_, transformation = get_transformation_at_MQS(e_sqd_log, [1, 1, 1])  # MQS from the above log

mol.atoms.transform(transformation)

from chimerax.map.molmap import molecule_map
mol_vol = molecule_map(session, mol.atoms, 4.0, grid_spacing=vol.data_origin_and_step()[1][0]/3)

# Manually change the surface threshold

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


zero_iter = 0
r = vol.subregion()
g = vol.region_grid(r)
g.array[:, :, :] = vol_matrix
g.name = vol.name + f"zero_{zero_iter}"

from chimerax.map.volume import volume_from_grid_data
v_zeroed = volume_from_grid_data(g, session)
v_zeroed.copy_settings_from(vol, copy_region=False, copy_colors=True)



# playground draft code ========

