from scipy.spatial.transform import Rotation as R
from DiffFit_bins import DiffFit_Binned_Transforms
from chimerax.geometry import Place
import numpy as np
from datetime import datetime
from math import pi
import math
from sklearn.cluster import Birch


fit_res = np.load("fit_res_filtered.npz")
fit_res_filtered = fit_res['fit_res_filtered']
# mol_center = fit_res['mol_center']
# mol_center = np.array([0.0, 0.0, 0.0])

mol_shift = fit_res_filtered[:, :3]
mol_q = fit_res_filtered[:, 3:7]

Total_fits = len(mol_shift)
# Total_fits = 1000 # to mimic other records in the reported table

print(f"Clustering {Total_fits} fits")
timer_start = datetime.now()

T = []
for i in range(Total_fits):
    shift = mol_shift[i]
    quat = mol_q[i]
    R_matrix = R.from_quat(quat).as_matrix()

    T_matrix = np.zeros([3, 4])
    T_matrix[:, :3] = R_matrix
    T_matrix[:, 3] = shift

    transformation = Place(matrix=T_matrix)
    T.append(transformation)

print(f"Convert to matrix time: {datetime.now() - timer_start}")
timer_start = datetime.now()

angle_tolerance = 6.0
shift_tolerance = 3.0

bfactor=2
print(f"bfactor: {bfactor}")

ChimeraX_clustering = True
if ChimeraX_clustering:
    b = DiffFit_Binned_Transforms(angle_tolerance * pi / 180, shift_tolerance, bfactor=bfactor)
    mol_transform_label = []
    unique_id = 0
    T_ID_dict = {}
    for i in range(Total_fits):
        ptf = T[i]
        in_cluster = b.one_in_cluster_transform(ptf)
        if in_cluster is None:
            b.add_transform(ptf)
            mol_transform_label.append(unique_id)
            T_ID_dict[id(ptf)] = unique_id
            unique_id = unique_id + 1
        else:
            mol_transform_label.append(T_ID_dict[id(in_cluster)])
            T_ID_dict[id(ptf)] = T_ID_dict[id(in_cluster)]

    print(f"Total unique: {unique_id}")

    print(f"ChimeraX bin clustering: {datetime.now() - timer_start}")
