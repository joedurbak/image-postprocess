import os
from zipfile import ZipFile
from shutil import copy2

input_dir = r'F:\Data'
os.chdir(input_dir)
output_dir = r'Z:\prime1\data'
compression_prefix = ''
compression_suffix = ''
directories_only = True
skip_names = ['20210713_to_20210719']
compression_folders = \
    [folder for folder in os.listdir(input_dir) if os.path.isdir(folder) and folder.startswith(compression_prefix) and folder.endswith(compression_suffix) and folder not in skip_names]
# compression_folders_abs = [os.path.join(input_dir, folder) for folder in compression_folders_rel]
zip_names = [folder+'.zip' for folder in compression_folders]
# zip_names_abs = [folder+'.zip' for folder in compression_folders_abs]
zip_names_moved = [os.path.join(output_dir, f) for f in zip_names]

for compression_folder, zip_name, zip_name_moved in zip(compression_folders, zip_names, zip_names_moved):
    print(compression_folder, zip_name, zip_name_moved)
    with ZipFile(zip_name, 'w') as zip_object:
        for folder_name, sub_folders, file_names in os.walk(compression_folder):
            for filename in file_names:
                # Create filepath of files in directory
                file_path = os.path.join(folder_name, filename)
                print('compressing {}'.format(file_path))
                # Add files to zip file
                zip_object.write(file_path, file_path)
    if os.path.isfile(zip_name):
        zip_size = os.path.getsize(zip_name)
        print('zip size {}'.format(zip_size))
        if zip_size > 10000:
            print('copying zipfile {} to {}'.format(zip_name, zip_name_moved))
            copy2(zip_name, zip_name_moved)
            copy_size = os.path.getsize(zip_name_moved)
            if os.path.getsize(zip_name_moved) == zip_size:
                print('deleting {}'.format(zip_name))
                os.remove(zip_name)
            else:
                print("copysize {} didn't match zipsize {}, skipping delete...".format(copy_size, zip_size))
        else:
            print('zipsize too small, skipping...')
    else:
        print("zipfile {} doesn't exist".format(zip_name))
