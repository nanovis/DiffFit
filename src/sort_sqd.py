import sys
sys.path.append('D:\\Research\\IPM\\PoseEstimation\\DiffFitViewer\\script')
from parse_log import cluster_and_sort_sqd
import numpy as np

e_sqd_log = np.load("D:\Research\IPM\PoseEstimation\DiffFitViewer\dev_data\output\dev_comp_domain_fit_3_domains\e_sqd_log.npy")

e_sqd_clusters_ordered = cluster_and_sort_sqd(e_sqd_log)

e_sqd_clusters_ordered_len = [len(cluster) for cluster in e_sqd_clusters_ordered]

