import cProfile
import pstats
from datetime import datetime

from scipy.spatial.transform import Rotation as R
from chimerax.geometry.bins import Binned_Transforms
from chimerax.geometry import Place
import numpy as np
from datetime import datetime
from math import pi
import math
import io


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
    

def main():
    fit_res = np.load("fit_res_filtered.npz")
    fit_res_filtered = fit_res['fit_res_filtered']
    
    mol_shift = fit_res_filtered[:, :3]
    mol_q = fit_res_filtered[:, 3:7]

    Total_fits = len(mol_shift)
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
    bfactor = 2
    print(f"bfactor: {bfactor}")

    ChimeraX_clustering = True
    if ChimeraX_clustering:
        b = DiffFit_Binned_Transforms(angle_tolerance * pi / 180, shift_tolerance)
        mol_transform_label = []
        unique_id = 0
        T_ID_dict = {}
        for i in range(len(mol_shift)):
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

if __name__ == "__main__":
    repetitions = 2  # Number of times to run the profiling
    profiler_stats = []

    for i in range(repetitions):
        profiler = cProfile.Profile()
        profiler.enable()
        main()
        profiler.disable()

        stats = pstats.Stats(profiler)
        profiler_stats.append(stats)

        print(f"--- Profiling run {i+1} complete ---")

    # Aggregate statistics
    aggregate_stats = profiler_stats[0]
    for i in range(1, repetitions):
        aggregate_stats.add(profiler_stats[i])

    # Print averaged stats
    stream = io.StringIO()
    aggregate_stats.stream = stream
    aggregate_stats.sort_stats('tottime')
    aggregate_stats.print_stats(20)

    print("--- Aggregated Profile Results ---")
    print(stream.getvalue())
