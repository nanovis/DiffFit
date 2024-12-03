import torch
import numpy as np
from .DiffAtomComp import quaternion_to_matrix_batch, normalize_coordinates_to_map_origin_torch

def unit_sphere_vertices(num_vertices):
    from chimerax.surface.shapes import sphere_geometry2
    sphere_vertices = sphere_geometry2(2*num_vertices-4)[0] # 128 points evenly distributed around a unit sphere centred on (0,0,0)
    return sphere_vertices


RANDOM_SEED=1985


def generate_q_shells(mol,
                      points_per_shell=8, max_rad=2.0, step=0.1, num_test_points=128,
                      clustering_iterations=5, include_h=False, randomize_shell_points=True, random_seed=RANDOM_SEED):
    '''
    Implementation of the map-model Q-score as described in Pintille et al. (2020): https://www.nature.com/articles/s41592-020-0731-1.

    If the model is well-fitted to the map, the Q-score is essentially an atom-by-atom estimate of "resolvability"
    (i.e. how much the map tells us about the true position of the atom). If the model is *not* well-fitted, then
    low Q-scores are good pointers to possible problem regions. In practice, of course, usually we have a mixture of
    both situations.

    This version of the algorithm has a few minor modifications compared to the original implementation aimed at
    improving overall speed. As a result the scores it returns are not identical to the original, typically differing by
    +/- 0.04 in individual atom scores and +/- 0.02 in residue averages. This difference can be explained by the different
    choice of test points, and reflects the underlying sampling uncertainty in the method.

    This implementation works as follows:

    - For each atom, define a set of shells in steps of `step` angstroms out to `max_rad` angstroms.
    - For each shell, try to find at least `points_per_shell` probe points closer to the test atom than
      any other atom:

      - For radii smaller than about half a bond length, first try a set of `points_per_shell`
        points evenly spread around the spherical surface.
      - For larger radii (or if this quick approach fails to find enough points on smaller radii),
        start with `num_test_points` evenly spread on the sphere surface, remove points closer to other atoms
        than the test atom. If more than `points_per_shell` points remain, perform up to `clustering_iterations`
        rounds of k-means clustering to find `points_per_shell` clusters, and choose the point nearest to the
        centroid of each cluster. If <= `points_per_shell` points remain, just do the best with what we have.
        By default, the "seed" centroids for each cluster are chosen pseudo-randomly from the input points.
        Using the same `random_seed` will give the same result each time; varying `random_seed` over multiple
        runs can be used to give an idea of the underlying uncertainty in the algorithm. If `randomize_shell_points`
        is False, the seed centroids will instead be the closest point (in spherical coordinates) to each of
        `points_per_shell` evenly-spaced points on a unit sphere. While this may intuitively seem preferable,
        in practice for tightly-packed atoms it leads to oversampling of the "junctions" with other atoms, and
        undersampling of the unhindered space.

    Returns:

    - a numpy array with q shells coordinates
    - radii
    '''
    from datetime import datetime
    global_timer_start = datetime.now()

    from chimerax.geometry import find_close_points, find_closest_points, Places
    from chimerax.atomic import Residues
    import numpy as np

    from chimerax.qscore import _kmeans

    pps_vertices = unit_sphere_vertices(points_per_shell)
    ref_sphere_vertices = unit_sphere_vertices(num_test_points)
    ref_sphere_vertices_large = unit_sphere_vertices(num_test_points * 4)
    ref_sphere_vertices_huge = unit_sphere_vertices(num_test_points * 20)

    radii = np.arange(0, max_rad + step / 2, step)

    query_atoms = mol.atoms

    if not include_h:
        query_atoms = query_atoms[query_atoms.element_names != 'H']


    query_coords = query_atoms.scene_coords


    query_atoms_center = []
    query_atoms_points = []
    not_full_shells = 0

    q_scores = []
    for i, a in enumerate(query_atoms):

        not_full_flag = False

        a_coord = a.scene_coord
        _, nearby_i = find_close_points([a_coord], query_coords, max_rad * 3)
        nearby_a = query_atoms[nearby_i]
        ai = nearby_a.index(a)
        nearby_coords = nearby_a.scene_coords
        shell_rad = step
        local_d_vals = {}

        shell_points = []

        j = 1
        while shell_rad < max_rad + step / 2:
            local_pps = (pps_vertices * shell_rad) + a_coord
            if shell_rad < 0.7:  # about half a C-C bond length
                # Try the quick way first (should succeed for almost all cases unless geometry is seriously wonky)
                i1, i2, near1 = find_closest_points(local_pps, nearby_coords, shell_rad * 1.5)
                closest = near1
                candidates = i1[closest == ai]
                if len(candidates) == points_per_shell:
                    shell_rad += step
                    j += 1

                    shell_points.append(local_pps)

                    continue

            local_sphere = (ref_sphere_vertices * shell_rad) + a_coord
            i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
            closest = near1
            candidates = i1[closest == ai]

            if len(candidates) < points_per_shell:

                local_sphere = (ref_sphere_vertices_large * shell_rad) + a_coord
                i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
                closest = near1
                candidates = i1[closest == ai]

                if len(candidates) < points_per_shell:
                    local_sphere = (ref_sphere_vertices_huge * shell_rad) + a_coord
                    i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
                    closest = near1
                    candidates = i1[closest == ai]

                    if len(candidates) < points_per_shell:
                        not_full_shells += 1
                        not_full_flag = True

                    else:
                        points = local_sphere[candidates]
                        if not randomize_shell_points:
                            labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                                local_pps, clustering_iterations)
                        else:
                            labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                               clustering_iterations, random_seed + j)

                        points = points[closest]

                else:
                    points = local_sphere[candidates]
                    if not randomize_shell_points:
                        labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                            local_pps, clustering_iterations)
                    else:
                        labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                           clustering_iterations, random_seed + j)

                    points = points[closest]

            else:
                points = local_sphere[candidates]
                if not randomize_shell_points:
                    labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                        local_pps, clustering_iterations)
                else:
                    labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                       clustering_iterations, random_seed + j)

                points = points[closest]


            shell_rad += step
            j += 1

            shell_points.append(points)

        if not_full_flag:
            continue

        query_atoms_center.append(a_coord)
        query_atoms_points.append(shell_points)

    if not_full_shells:
        print(f"Not_full_shells: {not_full_shells}")

    query_atoms_points_array = np.stack(
        [np.concatenate(atom_shell_points, axis=0) for atom_shell_points in query_atoms_points])
    query_atoms_center_np = np.stack(query_atoms_center)
    query_atoms_center_repeated = np.repeat(query_atoms_center_np[:, np.newaxis, :], points_per_shell, axis=1)
    q_shell_coords = np.concatenate([query_atoms_center_repeated, query_atoms_points_array], axis=1)

    print(f"Generate q shells timer: {datetime.now() - global_timer_start}")

    return q_shell_coords, radii


def min_max_d(v):
    m = v.data.full_matrix()
    mean, sd, _ = v.mean_sd_rms()
    max_d = min(mean + sd * 10, m.max())
    min_d = max(mean - sd, m.min())
    return min_d, max_d


def q_scores_for_clusters(centered_mol, volume, fit_res_clusters, fit_res_all,
                     ref_sigma=0.6, points_per_shell=8, device="cuda"):

    shifts = []
    quaternions = []

    for cluster_idx in range(len(fit_res_clusters)):
        mol_idx = int(fit_res_clusters[cluster_idx, 0])
        record_idx = int(fit_res_clusters[cluster_idx, 1])
        iter_idx = int(fit_res_clusters[cluster_idx, 2])

        shift = fit_res_all[mol_idx, record_idx, iter_idx, :3]
        quat = np.concatenate(
            ([fit_res_all[mol_idx, record_idx, iter_idx, 6]], -fit_res_all[mol_idx, record_idx, iter_idx, 3:6]),
            axis=0)

        shifts.append(shift)
        quaternions.append(quat)

    shifts = np.stack(shifts)
    quaternions = np.stack(quaternions)
    shifts = torch.tensor(shifts, device=device).float()
    quaternions = torch.tensor(quaternions, device=device).float()
    quaternions_matrices = quaternion_to_matrix_batch(quaternions.unsqueeze(0))


    q_shell_coords, radii = generate_q_shells(centered_mol)
    q_shell_coords = torch.tensor(q_shell_coords, device=device).float()
    q_shell_coords = q_shell_coords.reshape([-1, 3])

    vol_matrix = volume.full_matrix()
    vol_origin_and_step = volume.data_origin_and_step()
    target_origin = vol_origin_and_step[0]
    target_steps = vol_origin_and_step[1]
    target_no_negative = vol_matrix

    target = torch.tensor(target_no_negative, device=device).float()
    target_dim = target.shape

    import operator
    target_size = np.array(list(map(operator.mul, target_dim, target_steps)))

    target_size_x_y_z = [target_size[2], target_size[1], target_size[0]]
    target_size_x_y_z_tensor = torch.tensor(target_size_x_y_z, device=device).float()
    target_origin_tensor = torch.tensor(target_origin, device=device).float()

    target = target.unsqueeze(0).unsqueeze(0)


    min_d, max_d = min_max_d(volume)
    a = max_d - min_d
    b = min_d

    num_shells = len(radii)

    q_reference_gaussian = a * np.exp(-0.5 * (radii / ref_sigma) ** 2) + b
    q_ref = np.concatenate([[q_reference_gaussian[0]] * points_per_shell,
                            *[[q_reference_gaussian[j]] * points_per_shell for j in range(num_shells - 1)]])
    q_ref = torch.tensor(q_ref, device=device, dtype=torch.float32)
    q_ref -= q_ref.mean()


    q_scores = []
    for row in range(len(fit_res_clusters)):
        transformed_coords = torch.matmul(q_shell_coords, quaternions_matrices[:, row:row + 1, :, :])

        transformed_coords += shifts[row, :]

        q_shell_coords_normalized_to_target = normalize_coordinates_to_map_origin_torch(transformed_coords,
                                                                                        target_size_x_y_z_tensor,
                                                                                        target_origin_tensor)

        q_shell_coords_normalized_to_target = q_shell_coords_normalized_to_target.reshape([1, 1, -1, 168, 3])

        q_measure = torch.nn.functional.grid_sample(target, q_shell_coords_normalized_to_target, 'bilinear', 'border',
                                                    align_corners=True)

        q_measure.squeeze_()

        q_measure -= q_measure.mean(dim=-1, keepdim=True)

        inner_product = torch.matmul(q_measure, q_ref)
        q_measure_l2 = torch.norm(q_measure, p=2, dim=-1)
        q_ref_l2 = torch.norm(q_ref, p=2, dim=-1)
        q_score_torch = inner_product / (q_measure_l2 * q_ref_l2)

        q_scores.append(q_score_torch[~torch.isnan(q_score_torch)].mean())

    q_scores_tensor = torch.stack(q_scores)
    # top_10_values, top_10_indices = torch.topk(q_scores_tensor, k=10)
    #
    # # Display the results
    # print("Top 10 values:", top_10_values)
    # print("Indices of top 10 values + 1:", top_10_indices + 1)
    #
    # return top_10_values, top_10_indices

    return q_scores_tensor.cpu().numpy()


from concurrent.futures import ThreadPoolExecutor
def process_atom(atom_idx, query_atoms, query_coords,
                 pps_vertices,
                 ref_sphere_vertices, ref_sphere_vertices_large, ref_sphere_vertices_huge,
                 step, max_rad, points_per_shell,
                 clustering_iterations,
                 randomize_shell_points, random_seed):
    from chimerax.geometry import find_close_points, find_closest_points
    from chimerax.qscore import _kmeans

    not_full_flag = False

    a = query_atoms[atom_idx]
    a_coord = query_coords[atom_idx]
    _, nearby_i = find_close_points([a_coord], query_coords, max_rad * 3)
    nearby_a = query_atoms[nearby_i]
    ai = nearby_a.index(a)
    nearby_coords = nearby_a.scene_coords
    shell_rad = step
    local_d_vals = {}

    shell_points = []

    j = 1
    while shell_rad < max_rad + step / 2:
        local_pps = (pps_vertices * shell_rad) + a_coord
        if shell_rad < 0.7:  # about half a C-C bond length
            # Try the quick way first (should succeed for almost all cases unless geometry is seriously wonky)
            i1, i2, near1 = find_closest_points(local_pps, nearby_coords, shell_rad * 1.5)
            closest = near1
            candidates = i1[closest == ai]
            if len(candidates) == points_per_shell:
                shell_rad += step
                j += 1

                shell_points.append(local_pps)

                continue

        local_sphere = (ref_sphere_vertices * shell_rad) + a_coord
        i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
        closest = near1
        candidates = i1[closest == ai]

        if len(candidates) < points_per_shell:

            local_sphere = (ref_sphere_vertices_large * shell_rad) + a_coord
            i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
            closest = near1
            candidates = i1[closest == ai]

            if len(candidates) < points_per_shell:
                local_sphere = (ref_sphere_vertices_huge * shell_rad) + a_coord
                i1, i2, near1 = find_closest_points(local_sphere, nearby_coords, shell_rad * 1.5)
                closest = near1
                candidates = i1[closest == ai]

                if len(candidates) < points_per_shell:
                    not_full_flag = True

                else:
                    points = local_sphere[candidates]
                    if not randomize_shell_points:
                        labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                            local_pps, clustering_iterations)
                    else:
                        labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                           clustering_iterations, random_seed + j)

                    points = points[closest]

            else:
                points = local_sphere[candidates]
                if not randomize_shell_points:
                    labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                        local_pps, clustering_iterations)
                else:
                    labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                       clustering_iterations, random_seed + j)

                points = points[closest]

        else:
            points = local_sphere[candidates]
            if not randomize_shell_points:
                labels, closest = _kmeans.spherical_k_means_defined(points, a_coord, points_per_shell,
                                                                    local_pps, clustering_iterations)
            else:
                labels, closest = _kmeans.spherical_k_means_random(points, a_coord, points_per_shell,
                                                                   clustering_iterations, random_seed + j)

            points = points[closest]

        shell_rad += step
        j += 1

        shell_points.append(points)

    return None if not_full_flag else (a_coord, shell_points)


def generate_q_shells_parallel(mol,
                               points_per_shell=8, max_rad=2.0, step=0.1,
                               num_test_points=128, clustering_iterations=5,
                               include_h=False, randomize_shell_points=True, random_seed=RANDOM_SEED):
    from datetime import datetime

    from datetime import datetime
    global_timer_start = datetime.now()

    query_atoms = mol.atoms
    if not include_h:
        query_atoms = query_atoms[query_atoms.element_names != 'H']
    query_coords = query_atoms.scene_coords

    pps_vertices = unit_sphere_vertices(points_per_shell)
    ref_sphere_vertices = unit_sphere_vertices(num_test_points)
    ref_sphere_vertices_large = unit_sphere_vertices(num_test_points * 4)
    ref_sphere_vertices_huge = unit_sphere_vertices(num_test_points * 20)

    # Multithreading
    results = []
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                process_atom, atom_idx, query_atoms, query_coords, pps_vertices, ref_sphere_vertices,
                ref_sphere_vertices_large, ref_sphere_vertices_huge, step, max_rad, points_per_shell,
                clustering_iterations, randomize_shell_points, random_seed
            )
            for atom_idx in range(len(query_atoms))
        ]
        for future in futures:
            result = future.result()
            if result is not None:
                results.append(result)

    query_atoms_center, query_atoms_points = zip(*results)
    query_atoms_points_array = np.stack(
        [np.concatenate(atom_shell_points, axis=0) for atom_shell_points in query_atoms_points])
    query_atoms_center_np = np.stack(query_atoms_center)
    query_atoms_center_repeated = np.repeat(query_atoms_center_np[:, np.newaxis, :], points_per_shell, axis=1)
    q_shell_coords = np.concatenate([query_atoms_center_repeated, query_atoms_points_array], axis=1)

    print(f"Generate q shells timer: {datetime.now() - global_timer_start}")

    return q_shell_coords, np.arange(0, max_rad + step / 2, step)
