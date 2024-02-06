import os
from zipfile import ZipFile
from shutil import copy2, rmtree

folder_dir = r'F:\Data'
os.chdir(folder_dir)
zip_dir = r'Z:\prime1\data'
zipfiles = [f.replace('.zip', '') for f in os.listdir(zip_dir) if f.endswith('.zip')]
folders = [folder for folder in os.listdir(folder_dir) if os.path.isdir(folder)]

for folder in folders:
    if folder in zipfiles:
        print('removing {}'.format(folder))
        rmtree(folder)
