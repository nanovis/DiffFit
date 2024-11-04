# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import StructureArg, AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg, ColorArg, IntArg, FloatArg, StringArg, SaveFolderNameArg
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers
from chimerax.map.mapargs import MapArg

import os
from datetime import datetime

# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================


def dfit(session, mol, in_map,
         level=None,
         sim_res=5.0,
         num_positions=10,
         num_rotations=100,
         smooth_by="PyTorch iterative Gaussian",
         smooth_loops=3,
         kernel_sizes="[5, 5, 5]",
         Gaussian_mode="Gaussian with negative (shrink)",
         fit_atom_mode="Backbone",
         out_dir="DiffFit_out/interactive",
         device=None):
    """Fit a single structure into a volume map"""

    _save_results = True
    _out_dir_exist_ok = True
    _out_dir = out_dir

    use_level = in_map.maximum_surface_level
    if level is not None:
        use_level = level

    if _save_results:

        os.makedirs(_out_dir, exist_ok=_out_dir_exist_ok)
        with open(f"{_out_dir}/log.log", "a") as log_file:
            log_file.write(f"=======\n"
                           f"Wall clock time: {datetime.now()}\n"
                           f"-------\n"
                           f"Interactive mode\n"
                           f"Target Volume: {in_map.path}\n"
                           f"Structure: {mol.filename}\n"
                           f"Target Surface Threshold: {use_level}\n"
                           f"-------\n"
                           f"Sim-map resolution: {sim_res}\n"
                           f"# positions: {num_positions}\n"
                           f"# rotations: {num_rotations}\n"
                           f"Smooth by: \"{smooth_by}\"\n"
                           f"Smooth loops: {smooth_loops}\n"
                           f"Kernel sizes: \"{kernel_sizes}\"\n"
                           f"Gaussian mode: \"{Gaussian_mode}\"\n"
                           f"Fit atom mode: \"{fit_atom_mode}\"\n"
                           f"-------\n")


dfit_desc = CmdDesc(required=[("mol", StructureArg)],
                    keyword=[("in_map", MapArg),
                             ("level", FloatArg),
                             ("sim_res", FloatArg),
                             ("num_positions", IntArg),
                             ("num_rotations", IntArg),
                             ("smooth_by", StringArg),
                             ("smooth_loops", StringArg),
                             ("kernel_sizes", StringArg),
                             ("Gaussian_mode", StringArg),
                             ("fit_atom_mode", StringArg),
                             ("out_dir", SaveFolderNameArg),
                             ("device", StringArg)],
                    required_arguments=['in_map'])

# Example commands
# dfit #1 in #2
# dfit #1 in #2
#      level 0.7 sim_res 5.0
#      num_p 10 num_r 100
#      smooth_by smooth_loops kernel_sizes gaussian_mode
#      fit_atom_mode out_dir device
#
# dfit multi str_dir sim_dir in map
#      level 0.7 num_p 10 num_r 100
#      smooth_by smooth_loops smooth_weights kernel_sizes gaussian_mode
#      fit_atom_mode out_dir device
#      negative_space
#      learning_rate
#      n_iters



def dfit_disk(session, atoms, color, weighted=False, transformed=True, count=1):
    """Highlight the atoms nearest the center of mass of given atoms."""

    # ``session``     - ``chimerax.core.session.Session`` instance
    # ``atoms``       - ``chimerax.atomic.Atoms`` instance or None
    # ``color``       - ``chimerax.core.colors.Color` instance
    # ``weighted``    - boolean, whether to include atomic mass in calculation
    # ``transformed`` - boolean, use scene rather than original coordinates

    # Compute the center of mass first
    atoms, coords, cofm = _get_cofm(session, atoms, transformed, weighted)

    # Compute the distance of each atom from the cofm
    # using the NumPy vector norm function
    from numpy.linalg import norm
    distances = norm(coords - cofm, axis=1)

    # Sort the array and get the "count" indices to the closest atoms
    if count > len(atoms):
        count = len(atoms)
    from numpy import argsort
    atom_indices = argsort(distances)[:count]

    # Create a collection of atoms from the indices
    chosen = atoms[atom_indices]

    # Update their "colors".  Assigning a single value to an
    # array means assign the same value for all elements.
    chosen.colors = color.uint8x4()


dfit_disk_desc = CmdDesc(required=[("atoms", Or(AtomsArg, EmptyArg)),
                                   ("color", ColorArg)],
                         optional=[("weighted", BoolArg),
                                   ("transformed", BoolArg),
                                   ("count", Bounded(IntArg, 1, 50))])


# ==========================================================================
# Functions intended only for internal use by bundle
# ==========================================================================


def _get_cofm(session, atoms, transformed, weighted):
    # ``session``     - ``chimerax.core.session.Session`` instance
    # ``atoms``       - ``chimerax.atomic.Atoms`` instance
    # ``transformed`` - boolean, use scene rather than original coordinates
    # ``weighted``    - boolean, whether to include atomic mass in calculation

    # If user did not specify the list of atoms, use all atoms
    if atoms is None:
        from chimerax.core.commands import all_objects
        atoms = all_objects(session).atoms

    # We can use either transformed or untransformed coordinates.
    # Transformed coordinates are "scene coordinates", which
    # takes into account translation and rotation of individual
    # models.  Untransformed coordinates are the coordinates
    # read from the data files.
    if transformed:
        coords = atoms.scene_coords
    else:
        coords = atoms.coords

    # ``coords`` is a ``numpy`` float array of shape (N, 3)
    # If we want weighted center, we have to multiply coordinates
    # by the atomic mass
    if not weighted:
        cofm = coords.mean(axis=0)
    else:
        m = atoms.elements.masses
        c = coords * m[:,None]
        cofm = c.sum(axis=0) / m.sum()

    # To get the average coordinates, we use ``numpy.mean``
    # print("DEBUG: center of mass:", cofm)
    return atoms, coords, cofm