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

    @staticmethod
    def print_sorted_differences(array1, array2, tolerance=0.1):
        """
        Finds, sorts, and prints the indices, differences, and values of two arrays
        based on the magnitude of their differences.
        
        Parameters:
            array1 (np.ndarray): The first array to compare.
            array2 (np.ndarray): The second array to compare.
            tolerance (float): The acceptable difference to highlight.
            
        Prints:
            A sorted list of indices, differences, value1, and value2 from the arrays.
        """
        # Ensure the arrays have the same shape
        if array1.shape != array2.shape:
            raise ValueError("Arrays must have the same shape to compare.")
        
        # Compute the absolute differences
        differences = np.abs(array1 - array2)
        
        # Create a list of tuples (index, difference, value1, value2)
        results = [(index, diff, array1[index], array2[index]) 
                for index, diff in enumerate(differences)]
        
        # Sort the results by the difference in descending order
        sorted_results = sorted(results, key=lambda x: x[1], reverse=True)
        
        # Print the sorted differences
        print(f"{'Index':<8} {'Difference':<12} {'Array1':<8} {'Array2':<8}")
        print("-" * 40)
        for idx, diff, val1, val2 in sorted_results:
            if diff > tolerance:
                print(f"{idx:<8} {diff:<12.4f} {val1:<8.4f} {val2:<8.4f}")
            else:
                print(f"{idx:<8} {diff:<12.4f} {val1:<8.4f} {val2:<8.4f} (Within tolerance)")

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
            parts = line.split()
            if len(parts) >= 2:
                q_value_array.append(float(parts[1]))
        
        return np.array(q_value_array,dtype=float)

    def test_qvalue(self, universe, reference,q_value_array):
        qvalues = qvalue.qValue(universe, reference)
        qvalues.run()
        #self.print_sorted_differences(qvalues.qvalues, q_value_array, 0.1)
        assert np.allclose(qvalues.qvalues, q_value_array, atol=0.005)

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
