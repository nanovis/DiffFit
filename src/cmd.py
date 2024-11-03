# vim: set expandtab shiftwidth=4 softtabstop=4:

from chimerax.core.commands import CmdDesc      # Command description
from chimerax.atomic import AtomsArg            # Collection of atoms argument
from chimerax.core.commands import BoolArg      # Boolean argument
from chimerax.core.commands import ColorArg     # Color argument
from chimerax.core.commands import IntArg       # Integer argument
from chimerax.core.commands import EmptyArg     # (see below)
from chimerax.core.commands import Or, Bounded  # Argument modifiers


# ==========================================================================
# Functions and descriptions for registering using ChimeraX bundle API
# ==========================================================================


def dfit(session, atoms, weighted=False, transformed=True):
    """Report center of mass of given atoms."""

    # ``session``     - ``chimerax.core.session.Session`` instance
    # ``atoms``       - ``chimerax.atomic.Atoms`` instance or None
    # ``weighted``    - boolean, whether to include atomic mass in calculation
    # ``transformed`` - boolean, use scene rather than original coordinates

    atoms, coords, cofm = _get_cofm(session, atoms, transformed, weighted)
    session.logger.info("%s center of mass: %s" %
                        ("weighted" if weighted else "unweighted", cofm))


dfit_desc = CmdDesc(required=[("atoms", Or(AtomsArg, EmptyArg))],
                    optional=[("weighted", BoolArg),
                              ("transformed", BoolArg)])



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