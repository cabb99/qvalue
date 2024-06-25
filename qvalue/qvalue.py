"""
q-value analysis
================================================================

Analyze the q-value *q*, a measure of the fraction of native contacts 
conserved during a simulation. The native contacts are defined as the
contacts between the selected atoms in the reference structure. This module
can be considered as an extension of the contact analysis module in
MDAnalysis. It implements the following additional varieties of Q-value:

1. *Wolynes q*: Defines the q-value as the fraction of native contacts, where
    the native contacts are defined as the contacts between the selected atoms
    in the reference structure without any cutoff. Additionally, each contact
    is weighted by a gaussian function centered aroung the distance between the 
    atoms in the reference structure and with the standard deviation of the separation
    based on the sequence among the atoms in the reference structure.

The "fraction of native contacts" *q(t)* is a number between 0 and 1, where 1 means that
all native contacts are conserved and 0 means that no native contacts are conserved.

"""
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import distance_array
import numpy as np


class qValue(AnalysisBase):
    def __init__(self, u, select, refgroup, method="Wolynes", kwargs=None, **basekwargs):
        self.u = u
        super(qValue, self).__init__(self.u.trajectory, **basekwargs)

        self.fraction_kwargs = kwargs if kwargs is not None else {}

        if method == 'Wolynes':
            self.qvalue_function = qvalue_wolynes


        self.select = select
        
        # contacts formed in reference
        self.r0 = []
        self.native_contacts = []
        self.sequence_separation=[]


def qvalue_wolynes(rij, rijn, sigmaij):
    return np.exp(-(rij - rijn)**2 / (2.0 * sigmaij**2))