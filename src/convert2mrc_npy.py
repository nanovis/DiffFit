#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script for ChimeraX to extract each chain from a structure,
save them as separate PDB files, and generate MRC files for each chain.
"""

# Import necessary modules
import os, sys
from chimerax.core.commands import run
from datetime import datetime
import numpy as np

timer_start = datetime.now()

# Input parameters
structures_dir = sys.argv[1]  # Path to the input structure file
out_mrc_dir = sys.argv[2]  # Directory to save output mrc files
out_npy_dir = sys.argv[3]  # Directory to save output npy files
resolution = float(sys.argv[4])  # Resolution for simulated MRC files
gridSpacing = float(sys.argv[5])  # gridSpacing for simulated MRC files

for file_name in os.listdir(structures_dir):
    print(f"\n======= Processing {file_name} =======")
    full_path = os.path.join(structures_dir, file_name)
    # Check if the current path is a file and not a directory
    if os.path.isfile(full_path):
        # Open the input structure
        structure = run(session, f'open {full_path}')[0]

        # Get the base name of the input structure for naming output files
        structure_basename = os.path.basename(full_path).split('.')[0]

        # Save the structure's coordinates as a npy file
        npy_filename = f"{structure_basename}.npy"
        npy_filepath = os.path.join(out_npy_dir, npy_filename)
        atom_coordinates = []
        for chain in structure.chains:
            # Iterate through each residue in the chain to access its atoms
            for residue in chain.residues:
                # Check if the residue is not None
                if residue is not None:
                    for atom in residue.atoms:
                        # Append the atom's coordinates to our list
                        atom_coordinates.append(atom.coord)
        # Convert the list of coordinates to a NumPy array
        coordinates_array = np.array(atom_coordinates)
        np.save(npy_filepath, coordinates_array)

        # Generate and save the MRC file for the structure
        mrc_filename = f"{structure_basename}.mrc"
        mrc_filepath = os.path.join(out_mrc_dir, mrc_filename)

        vol = run(session, f'molmap #{structure.id[0]} {resolution} gridSpacing {gridSpacing}')
        run(session, f"save {mrc_filepath} #{vol.id[0]}")

        run(session, f"close #{structure.id[0]}")
        run(session, f"close #{vol.id[0]}")


print("Process completed.")

timer_stop = datetime.now()
print(f"Time elapsed: {timer_stop - timer_start}\n")
