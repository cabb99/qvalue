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
import csv

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
        with open(INFO, mode='r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                q_value_array.append(row)
        return q_value_array


    def test_qvalue(self, universe, reference):
        qvalues = qvalue.qValue(universe, reference)
        qvalues.run()
        print(qvalues.timeseries.shape)
        print(qvalues.qvalues.shape)
        test_qvalues = np.ones(20)

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
