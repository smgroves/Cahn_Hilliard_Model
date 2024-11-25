import os
import re

# Define the folder path
folder_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"


# Define a function to extract variables and rename files
def rename_files(folder_path):
    # Navigate to the folder
    for filename in os.listdir(folder_path):
        # Match filenames with the required prefixes
        if filename.startswith("discrete_norm_e") or filename.startswith("ave_mass"):
            # Extract the variables using regex
            match = re.match(
                r".*_(\d+)_(\d+)_(\d+(?:\.\d+)?(?:e[-+]?\d+)?)_dt_(\d+(?:\.\d+)?(?:e[-+]?\d+)?)",
                filename,
            )
            if match:
                Nx, timesteps, tol, dt = match.groups()
                # Construct the new filename
                prefix = "energy" if "discrete_norm_e" in filename else "mass"
                new_name = (
                    f"MG_{timesteps}_dt_{dt}_Nx_{Nx}_eps_0.015_tol_{tol}_{prefix}.txt"
                )
                # Define full paths
                old_path = os.path.join(folder_path, filename)
                new_path = os.path.join(folder_path, new_name)
                # Rename the file
                os.rename(old_path, new_path)
                print(f"Renamed: {old_path} -> {new_path}")
                # print(old_path)
                # print(new_path)


# Run the function
rename_files(folder_path)
