"""
Qvalue
"Calculates the fraction of native contacts conserved during the simulation"
"""

# Add imports here
from importlib.metadata import version

#__version__ = version("qvalue")
__version__ = "0.0.1"
from qvalue.qvalue import qValue

__all__ = ['Qvalue']