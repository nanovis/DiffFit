import os
import subprocess

# Define the paths
segmented_maps_dir = r"D:/GIT/DiffFit/case_study_data/7QIZ/100-500aa/segmented_maps"
structures_dir = r"D:/GIT/DiffFit/case_study_data/7QIZ/100-500aa/subunits_cif"
structures_sim_map_dir = r"D:/GIT/DiffFit/case_study_data/7QIZ/100-500aa/subunits_mrc"
out_base_dir = r"D:/GIT/DiffFit/case_study_data/7QIZ/100-500aa/out_cmd"
device = "cuda"
N_shifts = 100
target_surface_threshold = 0.02

# Make sure the output base directory exists
os.makedirs(out_base_dir, exist_ok=True)

# Loop through all files in the segmented_maps directory
for map_file in os.listdir(segmented_maps_dir):
    if map_file.endswith(".mrc"):
        # Construct the full path to the map file
        target_vol = os.path.join(segmented_maps_dir, map_file)
        
        # Define the output directory based on the map file name
        out_dir = os.path.join(out_base_dir, os.path.splitext(map_file)[0])
        os.makedirs(out_dir, exist_ok=True)
        
        # Construct the command
        cmd = [
            "python", "D:/GIT/DiffFit/src/DiffAtomComp.py",
            "--target_vol", target_vol,
            "--target_surface_threshold", str(target_surface_threshold),
            "--structures_dir", structures_dir,
            "--structures_sim_map_dir", structures_sim_map_dir,
            "--out_dir", out_dir,
            "--out_dir_exist_ok", "True",
            "--N_shifts", str(N_shifts),
            "--device", device
        ]
        
        # Execute the command
        print(f"Running command for {map_file}:")
        subprocess.run(cmd)