import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import cluster_and_sort_sqd, get_transformation_at_idx
import numpy as np
from chimerax.core.commands import run
import os

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains\e_sqd_log.npy")
e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)

mol_folder = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif"
mol_files = os.listdir(mol_folder)
# mol_files[idx] pairs with e_sqd_clusters_ordered[:][:, idx]

look_at_idx = 0
look_at_mol_idx, transformation = get_transformation_at_idx(e_sqd_clusters_ordered, look_at_idx)

mol_path = os.path.join(mol_folder, mol_files[look_at_mol_idx])
mol = run(session, f"open {mol_path}")[0]
mol.scene_position = transformation


# playground draft code ========
from chimerax.atomic import AtomicStructure
structures = session.models.list(type=AtomicStructure)
structures[0].scene_position = transformation

