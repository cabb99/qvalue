"""
Location of data files
======================

Use as ::

    from qvalue.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources
 
data_directory = importlib.resources.files("qvalue") / "data"
lammps_files = data_directory / 'lammps'
openmm_files = data_directory / 'openmm'


MDANALYSIS_LOGO = data_directory / "mda.txt"
DCD = openmm_files / 'movie.dcd'
PDB = openmm_files / 'native.pdb'
REF = openmm_files / 'crystal_structure-openmmawsem.pdb'
CIF = openmm_files / '1r69.cif'
INFO = openmm_files / 'info.dat'
DUMP = lammps_files / 'dump.lammpstrj'
QW = lammps_files / 'qw'
QO = lammps_files / 'qo'

