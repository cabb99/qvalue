"""
Unit and regression test for the qvalue package.
"""

# Import package, test suite, and other packages as needed
import qvalue
import pytest
import sys
import MDAnalysis as mda
from qvalue.data.files import DCD, PDB, REF, INFO
import numpy as np

class TestQValue(object):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(PDB, DCD)
    
    @pytest.fixture()
    def reference(self):
        return mda.Universe(REF)
    
    @pytest.fixture()
    def q_value_array(self):
        q_value_array = []
        with open(INFO, 'r') as infile:
            lines = infile.readlines()
        
        # Skip the first line (header)
        data_lines = lines[1:]
        
        for line in data_lines:
            # Split on *any* amount of whitespace
            parts = line.split()
            
            # Make sure there's at least two columns
            if len(parts) >= 2:
                q_value_array.append(float(parts[1]))
        
        return np.array(q_value_array,dtype=float)



    def test_qvalue(self, universe, reference,q_value_array):
        qvalues = qvalue.qValue(universe, reference)
        qvalues.run()
        print(qvalues.timeseries.shape)
        print(qvalues.qvalues.shape)
        print(q_value_array.shape)
        q_value_array

@pytest.fixture
def dcdfile():
    return qvalue.data.files.DCD

@pytest.fixture
def pdbfile():
    return qvalue.data.files.PDB

def test_import_dcd_file(dcdfile, pdbfile):
    mda.Universe(pdbfile, dcdfile)

def test_qvalue_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "qvalue" in sys.modules


def test_mdanalysis_logo_length(mdanalysis_logo_text):
    """Example test using a fixture defined in conftest.py"""
    logo_lines = mdanalysis_logo_text.split("\n")
    assert len(logo_lines) == 46, "Logo file does not have 46 lines!"
