import os
import json
from pathlib import Path


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

    # Load the target volume once
    volume_path = Path(config["target_volume"])

    output_directory = Path(config["output_directory"])

    for folder in config["folders"]:
        folder_path = Path(folder["path"])
        include_patterns = folder.get("include", [])
        exclude_patterns = folder.get("exclude", [])

        filtered_files = get_filtered_files(folder_path, include_patterns, exclude_patterns)

        for file_path in sorted(filtered_files):
            try:
                print(f"Processing {file_path}")
                out_sub_folder = Path(f"{output_directory}/{Path(file_path).stem}")
                out_sub_folder.mkdir(parents=True, exist_ok=True)
                # diff_atom_comp(
                #     target_vol_path=str(volume_path),
                #     target_surface_threshold=config.get("target_surface_threshold", 0.02),
                #     min_cluster_size=config.get("min_cluster_size", 100),
                #     structures_dir=str(file_path),
                #     fit_atom_mode=config.get("fit_atom_mode", "Backbone"),
                #     Gaussian_mode=config.get("Gaussian_mode", "Gaussian with negative (shrink)"),
                #     N_shifts=config.get("num_positions", 10),
                #     N_quaternions=config.get("num_rotations", 100),
                #     negative_space_value=config.get("negative_space", -0.5),
                #     learning_rate=config.get("learning_rate", 0.01),
                #     n_iters=config.get("n_iters", 201),
                #     out_dir=str(output_directory),
                #     out_dir_exist_ok=True,
                #     conv_loops=config.get("conv_loops", 3),
                #     conv_kernel_sizes=config.get("conv_kernel_sizes", [5, 5, 5]),
                #     conv_weights=config.get("conv_weights", [1.0, 1.0, 1.0]),
                #     device=config.get("gpu_device", "cuda:0")
                # )
                print(f"Completed processing {file_path}")
            except Exception as e:
                print(f"Error processing {file_path}: {e}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run DiffFit with configuration file.")
    parser.add_argument("--config", required=True, help="Path to the configuration file.")

    args = parser.parse_args()

    process_files(args.config)
