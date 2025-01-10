"""
q-value analysis
================================================================

Analyze the q-value *q*, a measure of the fraction of native contacts 
conserved during a simulation. The native contacts are defined as the
contacts between the selected atoms in the reference structure. This module
is an extension of the contact analysis module in MDAnalysis. 
It implements the following additional varieties of Q-value:

1. *Wolynes q*: Defines the q-value as the fraction of native contacts, where
    the native contacts are defined as the contacts between the selected atoms
    in the reference structure without any cutoff. Additionally, each contact
    is weighted by a gaussian function centered aroung the distance between the 
    atoms in the reference structure and with the standard deviation of the separation
    based on the sequence among the atoms in the reference structure.

    .. math::
    
        q(t) = \frac{1}{N} \sum_{i,j} e^{-\frac{(r_{ij} - r_{ij}^N)^2}{2\sigma_{ij}^2}}

    where :math:`r_{ij}` is the distance between atoms :math:`i` and :math:`j` in the

2. "Onuchic q": 

    

The "fraction of native contacts" *q(t)* is a number between 0 and 1, where 1 means that
all native contacts are conserved and 0 means that no native contacts are conserved.

"""
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import distance_array
import numpy as np


class qValue(AnalysisBase):
    def __init__(self, u, refgroup=None, selection=None, selectionB=None, method="Wolynes", pbc=True, cutoff_radius=4.5, kwargs=None, **basekwargs):
        self.u = u
        super(qValue, self).__init__(self.u.trajectory, **basekwargs)

        # Arguments for the qvalue function
        self.qvalue_kwargs = kwargs if kwargs is not None else {}

        # Select the qvalue function
        if method == 'Wolynes':
            self.qvalue_function = qvalue_wolynes
        else:
            if not callable(method):
                raise ValueError("method has to be callable or 'Wolynes'")
            self.qvalue_function = method


        # If the selection is not provided, use all atoms
        if selection is None:
            self.selection = u.atoms
        else:
            self.selection = selection

        # If the selection group B is not provided, use selection group A
        if selectionB is None:
            self.selectionB = self.selection
        else:
            self.selectionB = selectionB

        # Get dimension of box if pbc set to True
        self.pbc = pbc
        if self.pbc:
            self._get_box = lambda ts: ts.dimensions
        else:
            self._get_box = lambda ts: None

        # If the reference group is not provided, use the first frame of the trajectory
        if refgroup is None:
            self.refgroup = u.atoms
        else:
            self.refgroup = refgroup


        # Measure the distance of the contacts formed in the reference
        self.r0 = []
        self.native_contacts = []
        self.sequence_separation=[]

        # Measure the distance of the contacts in the reference structure using periodic boundary conditions
        self.r0.append(distance_array(self.refgroup.select(self.selectionA).positions, self.refgroup.select(self.selectionB).positions, box=self._get_box(self.refgroup.universe)))
        
        # Select the contacts formed in the reference
        self.native_contacts.append(self.r0[-1] <= cutoff_radius)

        # Measure the sequence separation of the contacts in the reference structure. 
        self.sequence_separation.append(np.abs(self.refgroup.select(self.selectionA).resids[:, np.newaxis] - self.refgroup.select(self.selectionB).resids) - 1)

        #If the chains are not continuous, the sequence separation is np.inf
        self.sequence_separation[-1][self.sequence_separation[-1] < 0] = np.inf
        
        def _prepare(self):
            self.timeseries = np.empty((self.n_frames, len(self.r0)+1))

        def _single_frame(self):
            self.timeseries[self._frame_index][0] = self._ts.frame
            
            # compute distance array for a frame
            d = distance_array(self.grA.positions, self.grB.positions,
                                box=self._get_box(self._ts))
            
            for i, (initial_contacts, r0) in enumerate(zip(self.initial_contacts,
                                                        self.r0), 1):
                # select only the contacts that were formed in the reference state
                r = d[initial_contacts]
                r0 = r0[initial_contacts]
                q = self.fraction_contacts(r, r0, **self.fraction_kwargs)
                self.timeseries[self._frame_index][i] = q

def qvalue_wolynes(rij, rijn, sigmaij):
    return np.exp(-(rij - rijn)**2 / (2.0 * sigmaij**2))