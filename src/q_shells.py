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
    from math import floor
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

    num_shells = int(floor(max_rad / step))

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

