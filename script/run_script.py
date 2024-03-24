import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import cluster_and_sort_sqd, look_at_cluster, look_at_MQS_idx
import numpy as np
from chimerax.core.commands import run
import os

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains\e_sqd_log.npy")
e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)

mol_folder = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif"

cluster_idx = 0
look_at_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)




# playground draft code ========

