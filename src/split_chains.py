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
input_model = sys.argv[1]  # Path to the input structure file
output_dir = sys.argv[2]  # Directory to save output files
resolution = float(sys.argv[3])  # Resolution for simulated MRC files
gridSpacing = float(sys.argv[4])  # gridSpacing for simulated MRC files


# Open the input structure
structure = run(session, f'open {input_model}')[0]

# Get the base name of the input structure for naming output files
structure_basename = os.path.basename(input_model).split('.')[0]

# Iterate over each chain in the structure
chain_id_name_list = []
for chain in structure.chains:
    chain_id = chain.chain_id
    print(f"\n======= Processing {chain_id} =======")

    chain_id_name = f"{chain_id}"

    if chain_id_name.upper() in chain_id_name_list:
        chain_id_name = f"{chain_id}_"

    chain_id_name_list.append(chain_id_name.upper())

    if chain.num_existing_residues <= 500:
        continue

    # Save the chain's coordinates as a npy file
    npy_filename = f"{structure_basename}_chain_{chain_id_name}.npy"
    npy_filepath = os.path.join(output_dir, npy_filename)

    atom_coordinates = []
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

    # Save the chain as a cif file
    chain_filename = f"{structure_basename}_chain_{chain_id_name}.cif"
    chain_filepath = os.path.join(output_dir, chain_filename)

    run(session, f"select #{structure.id[0]}/{chain_id}")
    run(session, f"save {chain_filepath} selectedOnly true")
    run(session, "select clear")

    # Generate and save the MRC file for the chain
    mrc_filename = f"{structure_basename}_chain_{chain_id_name}.mrc"
    mrc_filepath = os.path.join(output_dir, mrc_filename)

    vol = run(session, f'molmap #{structure.id[0]}/{chain_id} {resolution} gridSpacing {gridSpacing}')
    run(session, f"save {mrc_filepath} #{vol.id[0]}")
    run(session, f"close #{vol.id[0]}")

print("Process completed.")

timer_stop = datetime.now()
print(f"Time elapsed: {timer_stop - timer_start}\n")
