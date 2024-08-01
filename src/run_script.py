import sys
sys.path.append('D:\\GIT\\DiffFit\\src')
from parse_log import cluster_and_sort_sqd, cluster_and_sort_sqd_fast, look_at_cluster, look_at_MQS_idx, animate_MQS, animate_cluster
from parse_log import zero_cluster_density
import numpy as np
from chimerax.core.commands import run
import os

vol_path = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\density2.mrc"
vol = run(session, f"open {vol_path}")[0]

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains_10s20q\e_sqd_log.npy")
# e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)
e_sqd_clusters_ordered = cluster_and_sort_sqd_fast(e_sqd_log)

mol_folder = "D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\input\domain_fit_demo_3domains\subunits_cif"

cluster_idx = 0
look_at_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)  # read from log to get MQS

look_at_MQS_idx(e_sqd_log, mol_folder, [0,0,0], session)


animate_MQS(e_sqd_log, mol_folder, [1, 87, 2], session)
animate_cluster(e_sqd_clusters_ordered, mol_folder, cluster_idx, session)


zero_cluster_density(vol, e_sqd_clusters_ordered, mol_folder, cluster_idx, session, res=4.0, zero_iter=0)

# TODO: res should be specified by the user for once
# TODO: zero_iter should go up in later iterations, this is only for better log/info, does not matter the computation



# playground draft code ========

