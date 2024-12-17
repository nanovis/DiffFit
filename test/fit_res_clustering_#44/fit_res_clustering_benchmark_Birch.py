from scipy.spatial.transform import Rotation as R
from chimerax.geometry.bins import Binned_Transforms
from chimerax.geometry import Place
import numpy as np
from datetime import datetime
from math import pi
import math
from sklearn.cluster import Birch

class DiffFit_Binned_Transforms(Binned_Transforms):
    def __init__(self, angle, translation, center=(0, 0, 0), bfactor=2):
        super().__init__(angle, translation, center, bfactor)

    def one_in_cluster_transform(self, tf):

        a, x, y, z = c = self.bin_point(tf)
        clist = self.bins.close_objects(c, self.spacing)
        if len(clist) == 0:
            return None

        itf = tf.inverse()
        d2max = self.translation * self.translation
        for ctf in clist:
            cx, cy, cz = ctf * self.center
            dx, dy, dz = x - cx, y - cy, z - cz
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= d2max:
                dtf = ctf * itf
                a = dtf.rotation_angle()
                if a < self.angle:
                    return ctf

        return None


fit_res = np.load("fit_res_filtered.npz")
fit_res_filtered = fit_res['fit_res_filtered']
# mol_center = fit_res['mol_center']
# mol_center = np.array([0.0, 0.0, 0.0])

mol_shift = fit_res_filtered[:, :3]
mol_q = fit_res_filtered[:, 3:7]

Total_fits = len(mol_shift)
# Total_fits = 10000 # to mimic other records in the reported table

print(f"Clustering {Total_fits} fits")
timer_start = datetime.now()

def q2_unit_coord(Q):
    rotations = [R.from_quat(q) for q in Q]

    up = np.array([0, 1, 0])
    right = np.array([1, 0, 0])

    rotated_up = np.array([rot.apply(up) for rot in rotations])
    rotated_right = np.array([rot.apply(right) for rot in rotations])

    # return np.concatenate((rotated_up, rotated_right), axis=-1)
    return rotated_up, rotated_right

angle_tolerance = 6.0
shift_tolerance = 3.0

# cluster by Birch
timer_start = datetime.now()

q_coord_radius_tolerance = math.sin(math.radians(angle_tolerance / 2.0))

timer_start_tmp = datetime.now()
r_up, r_right  = q2_unit_coord(mol_q)
print(f"q2_unit_coord: {datetime.now() - timer_start_tmp}")
timer_start_tmp = datetime.now()

cluster_shift = Birch(threshold=shift_tolerance / 2.0, n_clusters=None).fit(mol_shift)
print(f"mol_shift cluster: {datetime.now() - timer_start_tmp}")
timer_start_tmp = datetime.now()

cluster_r_up = Birch(threshold=q_coord_radius_tolerance, n_clusters=None).fit(r_up)
print(f"r_up cluster: {datetime.now() - timer_start_tmp}")
timer_start_tmp = datetime.now()

cluster_r_right = Birch(threshold=q_coord_radius_tolerance, n_clusters=None).fit(r_right)
print(f"r_right cluster: {datetime.now() - timer_start_tmp}")
timer_start_tmp = datetime.now()

mol_cluster_label = np.concatenate((cluster_shift.labels_.reshape([-1, 1]),
                                     cluster_r_up.labels_.reshape([-1, 1]),
                                     cluster_r_right.labels_.reshape([-1, 1])), axis=-1)
unique_labels, indices, counts = np.unique(mol_cluster_label, axis=0, return_inverse=True,
                                                   return_counts=True)

print(f"Total unique: {len(unique_labels)}")

print(f"Birch clustering: {datetime.now() - timer_start}")