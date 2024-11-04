# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import StructureArg, AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg, ColorArg, IntArg, FloatArg, StringArg, SaveFolderNameArg
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers
from chimerax.map.mapargs import MapArg

import os
import numpy as np
import torch
from datetime import datetime

from .DiffAtomComp import diff_atom_comp, cluster_and_sort_sqd_fast, diff_fit, conv_volume, numpy2tensor, \
    linear_norm_tensor
import ast
from chimerax.core.commands import run

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

    _use_level = in_map.maximum_surface_level
    if level is not None:
        _use_level = level

    _use_device = "cpu"
    if torch.cuda.is_available():
        _use_device = "cuda:0"
    if device is not None:
        _use_device = device

    if _save_results:

        os.makedirs(_out_dir, exist_ok=_out_dir_exist_ok)
        with open(f"{_out_dir}/log.log", "a") as log_file:
            log_file.write(f"=======\n"
                           f"Wall clock time: {datetime.now()}\n"
                           f"-------\n"
                           f"Interactive mode\n"
                           f"Target Volume: {in_map.path}\n"
                           f"Structure: {mol.filename}\n"
                           f"Target Surface Threshold: {_use_level}\n"
                           f"-------\n"
                           f"Sim-map resolution: {sim_res}\n"
                           f"# positions: {num_positions}\n"
                           f"# rotations: {num_rotations}\n"
                           f"Smooth by: \"{smooth_by}\"\n"
                           f"Smooth loops: {smooth_loops}\n"
                           f"Kernel sizes: \"{kernel_sizes}\"\n"
                           f"Gaussian mode: \"{Gaussian_mode}\"\n"
                           f"Fit atom mode: \"{fit_atom_mode}\"\n"
                           f"Device: \"{_use_device}\"\n"
                           f"-------\n")

    single_fit_timer_start = datetime.now()

    # Prepare mol and vol
    vol_matrix = in_map.full_matrix()

    # Copy vol and make it clean after thresholding
    vol_copy = in_map.writable_copy()
    vol_copy_matrix = vol_copy.data.matrix()
    vol_copy_matrix[vol_copy_matrix < _use_level] = 0
    vol_copy.data.values_changed()

    # Smooth the volume
    volume_conv_list = _create_volume_conv_list(session,
                                                vol_copy,
                                                smooth_by, smooth_loops, kernel_sizes, Gaussian_mode,
                                                _use_device)
    vol_copy.delete()

    # Apply the user's transformation and center mol
    from chimerax.geometry import Place
    mol.atoms.transform(mol.position)
    mol_center = mol.atoms.coords.mean(axis=0)
    transform = Place(origin=-mol_center)
    mol.atoms.transform(transform)
    mol.position = Place()

    # Simulate a map for the mol
    from chimerax.map.molmap import molecule_map
    mol_vol = molecule_map(session, mol.atoms, sim_res,
                           grid_spacing=in_map.data.step[0])

    input_coords = None
    if fit_atom_mode == "Backbone":
        backbone_atoms = ['N', 'CA', 'C', 'O']
        is_backbone = np.isin(mol.atoms.names, backbone_atoms)

        input_coords = mol.atoms.scene_coords[is_backbone]
    elif fit_atom_mode == "All":
        input_coords = mol.atoms.scene_coords

    # Fit
    timer_start = datetime.now()

    if _save_results:
        with open(f"{_out_dir}/log.log", "a") as log_file:
            log_file.write(f"DiffFit optimization starts: {timer_start}\n")

    (_,
     _,
     mol_paths,
     mol_centers,
     fit_result) = diff_fit(
        volume_conv_list,
        in_map.path,
        _use_level,
        in_map.data.step,
        in_map.data.origin,
        10,
        [input_coords],
        mol.filename,
        [(mol_vol.full_matrix(), mol_vol.data.step, mol_vol.data.origin)],
        N_shifts=num_positions,
        N_quaternions=num_rotations,
        save_results=_save_results,
        out_dir=_out_dir,
        out_dir_exist_ok=_out_dir_exist_ok,
        device=_use_device
    )
    timer_stop = datetime.now()
    print(f"\nDiffFit optimization time elapsed: {timer_stop - timer_start}\n\n")

    if _save_results:
        with open(f"{_out_dir}/log.log", "a") as log_file:
            log_file.write(f"-------\n"
                           f"DiffFit optimization time elapsed: {timer_stop - timer_start}\n")

    mol_vol.delete()


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


def _create_volume_conv_list(session, vol, smooth_by, smooth_loops, smooth_kernel_sizes, Gaussian_mode, device, negative_space_value=-0.5):
    # From here on, there are three strategies for utilizing gaussian smooth
    # 1. with increasing sDev on the same input volume
    # 2. with the same sDev iteratively
    # Combine 1 & 2
    # Need to do experiment to see which one is better

    volume_conv_list = [None] * (smooth_loops + 1)
    volume_conv_list[0] = vol.full_matrix()

    if smooth_by == "PyTorch iterative Gaussian":
        volume_conv_list[0], _ = numpy2tensor(volume_conv_list[0], device)
        volume_conv_list[0] = linear_norm_tensor(volume_conv_list[0])
        volume_conv_list = conv_volume(volume_conv_list[0],
                                       device,
                                       smooth_loops,
                                       ast.literal_eval(smooth_kernel_sizes),
                                       negative_space_value=negative_space_value,
                                       kernel_type="Gaussian",
                                       mode=Gaussian_mode)
        volume_conv_list = [v.squeeze().detach().cpu().numpy() for v in volume_conv_list]
    elif smooth_by == "ChimeraX incremental Gaussian":
        for conv_idx in range(1, smooth_loops + 1):
            vol_gaussian = run(session, f"volume gaussian #{vol.id[0]} sDev {conv_idx}")

            vol_device, _ = numpy2tensor(vol_gaussian.full_matrix(), device)
            vol_device = linear_norm_tensor(vol_device)

            eligible_volume_tensor = vol_device > 0.0
            vol_device[~eligible_volume_tensor] = negative_space_value

            volume_conv_list[conv_idx] = vol_device.squeeze().detach().cpu().numpy()

            vol_gaussian.delete()
    elif smooth_by == "ChimeraX iterative Gaussian":
        kernel_sizes = ast.literal_eval(smooth_kernel_sizes)
        vol_current = vol
        for conv_idx in range(1, smooth_loops + 1):
            vol_gaussian = run(session, f"volume gaussian #{vol_current.id[0]} sDev {kernel_sizes[conv_idx - 1]}")

            if conv_idx > 1:
                vol_current.delete()

            vol_device, _ = numpy2tensor(vol_gaussian.full_matrix(), device)
            vol_device = linear_norm_tensor(vol_device)

            eligible_volume_tensor = vol_device > 0.0
            vol_device[~eligible_volume_tensor] = negative_space_value

            volume_conv_list[conv_idx] = vol_device.squeeze().detach().cpu().numpy()

            vol_current = vol_gaussian

        vol_current.delete()


    return volume_conv_list