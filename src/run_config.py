import os
import json
from pathlib import Path
import numpy as np
from difffit import (process_volume,
                     parse_precision,
                     prepare_atoms,
                     optimize_fitting)


def get_filtered_files(folder_path, include_patterns, exclude_patterns):
    """
    Get files from a folder, applying inclusion and exclusion rules.
    :param folder_path: Path to the folder.
    :param include_patterns: List of patterns to include (e.g., *.pdb).
    :param exclude_patterns: List of filenames to exclude.
    :return: List of filtered file names.
    """
    included_files = set()

    # Add files matching include patterns
    for pattern in include_patterns:
        included_files.update(folder_path.glob(pattern))

    # Convert to set of strings for easier exclusion
    included_files = {str(file) for file in included_files}

    # Remove files explicitly excluded
    excluded_files = {str(folder_path / name) for name in exclude_patterns}
    final_files = included_files - excluded_files

    return [Path(file) for file in sorted(final_files)]

def process_files(config_path):
    """
    Main function to process files based on the configuration.
    :param config_path: Path to the configuration JSON file.
    """
    with open(config_path, 'r') as config_file:
        config = json.load(config_file)

    print(f"Processing volume: {Path(config['target_volume'])}")
    (target_gaussian_conv_list,
     target,
     target_no_negative,
     target_size_x_y_z_tensor,
     target_origin_tensor,
     sampled_coords) = process_volume(target_vol_path=Path(config["target_volume"]),
                                      target_surface_threshold=config.get("target_surface_threshold"),
                                      N_shifts=config.get("num_positions", 30),
                                      negative_space_value=config.get("negative_space", -0.5),
                                      conv_loops=config.get("conv_loops", 3),
                                      conv_kernel_sizes=config.get("conv_kernel_sizes", [5, 5, 5]),
                                      conv_weights=config.get("conv_weights", [1.0, 1.0, 1.0]),
                                      Gaussian_mode=config.get("Gaussian_mode", "Gaussian with negative (shrink)"),
                                      device=config.get("gpu_device", "cuda:0"),
                                      precision=parse_precision(config.get("precision", "float32"))
                                      )

    output_directory = Path(config["output_directory"])

    for folder in config["structure_folders"]:
        folder_path = Path(folder["path"])
        include_patterns = folder.get("include", [])
        exclude_patterns = folder.get("exclude", [])

        filtered_files = get_filtered_files(folder_path, include_patterns, exclude_patterns)

        for structure_path in sorted(filtered_files):
            try:
                out_sub_folder = Path(f"{output_directory}/{Path(structure_path).stem}")
                out_sub_folder.mkdir(parents=True, exist_ok=True)

                (atom_coords_list,
                 mol_centers,
                 mol_num_atoms) = prepare_atoms(structure_path, config.get("fit_atom_mode", "Backbone"))

                e_sqd_log = optimize_fitting(
                    target,
                    target_gaussian_conv_list,
                    atom_coords_list,
                    sampled_coords,
                    target_size_x_y_z_tensor,
                    target_origin_tensor,
                    config.get("num_rotations", 100),
                    config.get("num_positions", 30),
                    conv_loops=config.get("conv_loops", 3),
                    conv_weights=config.get("conv_weights", [1.0, 1.0, 1.0]),
                    num_molecules=1,
                    n_iters=config.get("n_iters", 101),
                    learning_rate=config.get("learning_rate", 0.01),
                    device=config.get("gpu_device", "cuda:0"),
                    precision=parse_precision(config.get("precision", "float32")),
                    out_dir=out_sub_folder,
                )

                print(f"Completed processing {structure_path}")

                np.savez_compressed(f"{out_sub_folder}/fit_res.npz",
                                    target_vol_path=config['target_volume'],
                                    target_surface_threshold=config.get("target_surface_threshold"),
                                    mol_paths=[str(structure_path)],
                                    mol_num_atoms=mol_num_atoms,
                                    opt_res=e_sqd_log.detach().cpu().numpy())
            except Exception as e:
                print(f"Error processing {structure_path}: {e}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run DiffFit with configuration file.")
    parser.add_argument("--config", required=True, help="Path to the configuration file.")

    args = parser.parse_args()

    from datetime import datetime
    timer_start = datetime.now()
    process_files(args.config)
    print(f"Time elapsed: {datetime.now() - timer_start}\n\n")

