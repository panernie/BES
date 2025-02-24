# BES
Binding energy screening of small molecules and target materials

This tool is based on rdkit, unidock tool implementation, before using this code, please make sure that these programs have been installed properly.

rdkit==2022.3.3 (https://github.com/rdkit/rdkit)

Uni-Dock v1.1 (https://github.com/dptech-corp/Uni-Dock)

Instructions for use:

1. Structure generation

batch_smiles2pdbqt_build.py is used to generate pdbqt files based on SMILES of small molecules, please prepare a CSV file storing the names of small molecules and smiles before using this code.

python batch_smiles2pdbqt_build.py is ready to run, store up to 5000 small molecule pdbqt files in each folder.


2. Batch docking

Batch docking can be performed in batch_Uni-dock.ipynb. config1.json holds the configuration files required for Uni-dock docking, including box size and box coordinates


