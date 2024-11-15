import os
import shutil
import filecmp

def combine_folders(source1, source2, destination):
    # Ensure the destination directory exists
    if not os.path.exists(destination):
        os.makedirs(destination)
    
    # Open the output file to log duplicate files
    with open('/nonlinear_multigrid/utils/output.txt', 'w') as output_file:
        
        # Get the set of all unique project names in both source folders, filtering by "phi"
        projects = {p for p in os.listdir(source1) if p.startswith('phi')} | \
                   {p for p in os.listdir(source2) if p.startswith('phi')}

        for project in projects:
            project_source1 = os.path.join(source1, project)
            project_source2 = os.path.join(source2, project)
            project_dest = os.path.join(destination, project)
            
            # Ensure the project directory exists in the destination
            if not os.path.exists(project_dest):
                os.makedirs(project_dest)
            
            # Copy files from the first source directory
            if os.path.exists(project_source1):
                for file_name in os.listdir(project_source1):
                    src_file = os.path.join(project_source1, file_name)
                    dest_file = os.path.join(project_dest, file_name)
                    if os.path.isfile(src_file):
                        if not os.path.exists(dest_file):
                            shutil.copy(src_file, dest_file)
                        else:
                            # Check if the files are identical
                            if filecmp.cmp(src_file, dest_file, shallow=False):
                                output_file.write(f"Duplicate file: {src_file} -> {dest_file}\n")
                            else:
                                output_file.write(f"Conflict (different content): {src_file} -> {dest_file}\n")

            # Copy files from the second source directory
            if os.path.exists(project_source2):
                for file_name in os.listdir(project_source2):
                    src_file = os.path.join(project_source2, file_name)
                    dest_file = os.path.join(project_dest, file_name)
                    if os.path.isfile(src_file):
                        if not os.path.exists(dest_file):
                            shutil.copy(src_file, dest_file)
                        else:
                            # Check if the files are identical
                            if filecmp.cmp(src_file, dest_file, shallow=False):
                                output_file.write(f"Duplicate file: {src_file} -> {dest_file}\n")
                            else:
                                output_file.write(f"Conflict (different content): {src_file} -> {dest_file}\n")


# Define the source folders and destination folder
source_folder1 = '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/radii_lineplots_kymographs/domain_0_2_from_rivanna_kymographs'
source_folder2 = '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/radii_lineplots_kymographs/domain_0_2_from_rivanna_kymographs/new'
destination_folder = '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/radii_lineplots_kymographs/domain_0_2_from_rivanna_kymographs/combined_folder'

# Combine the folders
combine_folders(source_folder1, source_folder2, destination_folder)

print("Folder combining is complete. Check output.txt for details.")