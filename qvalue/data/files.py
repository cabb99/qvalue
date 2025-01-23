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
DCD = openmm_files / '1R69'/ 'movie.dcd'
PDB = openmm_files / '1R69'/ 'native.pdb'
REF = openmm_files / '1R69'/ 'crystal_structure-openmmawsem.pdb'
CIF = openmm_files / '1R69'/ '1r69.cif'
INFO = openmm_files / '1R69' / 'info.dat'
DUMP = lammps_files / '1R69' /'dump.lammpstrj'
QW = lammps_files / '1R69' / 'qw'
QO = lammps_files / '1R69' / 'qo'
QINTERFACE = lammps_files / 'p50p65' / 'q_interface_cb'
QINTERFACE_REF = lammps_files / 'p50p65' / 'Structure_p50p65.pdb'
QINTERFACE_DUMP = lammps_files / 'p50p65' / 'dump.lammpstrj'

