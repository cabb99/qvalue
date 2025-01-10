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
    """
    Calculate Q-values (fraction of native contacts) across a trajectory.

    This class compares distances between pairs of atoms (or residues) in the
    current trajectory against a reference structure. It can compute Q-values
    using one of several methods (e.g. Wolynes or user-defined functions).
    
    Parameters
    ----------
    universe : MDAnalysis.core.universe.Universe
        The universe (trajectory) on which Q-values are calculated.
    reference_universe : MDAnalysis.core.universe.Universe, optional
        A separate reference structure/universe. If not provided, the first
        frame of ``u`` is used as the reference.
    primary_selection : str, MDAnalysis.core.groups.AtomGroup, optional
        A selection string (e.g., ``'name CA'``) or an AtomGroup for the
        primary group of atoms. If None, uses all atoms in ``u``.
    secondary_selection : str, MDAnalysis.core.groups.AtomGroup, optional
        A selection string or an AtomGroup for the secondary group of atoms.
        If None, it defaults to the same selection as ``primary_selection``.
        This allows you to compute either intra-selection contacts or
        inter-selection contacts.
    q_method : {'Wolynes'} or callable, optional
        Method to compute the Q-values. By default, the "Wolynes" method is
        used. If a callable is provided, it must be a function that takes
        distances and returns a Q-value or fraction of native contacts.
    use_pbc : bool, optional
        Whether to apply periodic boundary conditions when computing distances.
        Default is True.
    contact_cutoff : float, optional
        Distance (in Å) below which two atoms are considered in contact when
        building the reference contact map. Default is 4.5 Å.
    method_kwargs : dict, optional
        Additional keyword arguments for the Q-value method (e.g., parameters
        for Wolynes or user-defined methods).
    **basekwargs : dict
        Additional keyword arguments passed to the :class:`MDAnalysis.analysis.base.AnalysisBase`
        constructor (e.g., `start`, `stop`, `step`).

    Attributes
    ----------
    qvalue_function : callable
        The function used internally to compute Q-values each frame.
    pbc : bool
        Whether to apply periodic boundary conditions.
    selectionA : MDAnalysis.core.groups.AtomGroup or str
        The primary AtomGroup (or selection string).
    selectionB : MDAnalysis.core.groups.AtomGroup or str
        The secondary AtomGroup (or selection string). Defaults to ``selectionA``.
    refgroup : MDAnalysis.core.groups.AtomGroup
        The reference AtomGroup. If not provided, uses the first frame of ``u``.
    r0 : list of ndarray
        Reference distances (in the reference structure) used for Q-value
        calculations.
    native_contacts : list of ndarray
        Boolean masks of which contacts are considered "native" in the reference.
    sequence_separation : list of ndarray
        The sequence separation (in residue index space) for the reference
        contacts. Used in certain Q-value calculations.

    Notes
    -----
    - The Q-value calculation compares distances in each frame to the reference
      distances (``r0``). 
    - If you define your own method in ``q_method``, make sure it handles input
      arrays of distances appropriately.
    - By default, sequence separation values < 0 are replaced with ``np.inf``
      to account for discontinuous residue numbering or multiple chains.

    Examples
    --------
    >>> from MDAnalysis import Universe
    >>> from qvalue_module import qValue
    >>> u = Universe("trajectory.psf", "trajectory.dcd")
    >>> # Use the first frame of u as the reference
    >>> q_calc = qValue(u, primary_selection='name CA')
    >>> q_calc.run()
    >>> # Access the final Q-values
    >>> final_q_values = q_calc.timeseries

    See Also
    --------
    MDAnalysis.analysis.base.AnalysisBase : The base class for MDAnalysis analyses.
    """
        
    def __init__(self, 
                 universe, 
                 reference_universe=None, 
                 primary_selection='name CA', 
                 secondary_selection=None, 
                 q_method="Wolynes", 
                 use_pbc=True, 
                 contact_cutoff=9.5,
                 min_seq_sep=3,
                 method_kwargs=None, 
                 **basekwargs):
        """
        Initialize the Q-value analysis.

        This constructor sets up the references, selections, method of
        calculation, and other necessary attributes before performing the
        Q-value computation in :meth:`run`.

        Parameters
        ----------
        universe : MDAnalysis.core.universe.Universe
            The trajectory (Universe) on which Q-values will be computed.
        reference_universe : MDAnalysis.core.universe.Universe, optional
            A reference Universe for contacts. If None, the first frame of ``u``
            is used as the reference.
        primary_selection : str or MDAnalysis.core.groups.AtomGroup, optional
            A selection or an AtomGroup to define the primary set of atoms.
        secondary_selection : str or MDAnalysis.core.groups.AtomGroup, optional
            A selection or an AtomGroup to define the secondary set of atoms.
            If None, uses the same as the primary selection.
        q_method : {'Wolynes'} or callable, optional
            A built-in method name or a user-supplied function to compute Q-values.
        use_pbc : bool, optional
            Whether to apply periodic boundary conditions when computing distances.
        contact_cutoff : float, optional
            Distance cutoff (in Å) for defining native contacts.
        method_kwargs : dict, optional
            Additional arguments for the Q-value method.
        **basekwargs
            Keyword arguments passed to the :class:`MDAnalysis.analysis.base.AnalysisBase`
            constructor, such as `start`, `stop`, `step`.

        Raises
        ------
        ValueError
            If the provided `q_method` is not 'Wolynes' or a callable function.
        """
        self.universe = universe
        super(qValue, self).__init__(self.universe.trajectory, **basekwargs)

        # Arguments for the qvalue function
        self.method_kwargs = method_kwargs if method_kwargs is not None else {}

        # Select the qvalue function
        if q_method == 'Wolynes':
            self.method_function = qvalue_wolynes
        else:
            if not callable(q_method):
                raise ValueError("method has to be callable or 'Wolynes'")
            self.method_function = q_method

        # If the reference group is not provided, use the first frame of the trajectory
        if reference_universe is None:
            self.reference_universe = universe
        else:
            self.reference_universe = reference_universe


        # If the selection is not provided, use all atoms
        if primary_selection is None:
            self.selection_a = universe.atoms
            self.reference_selection_a = self.reference_universe.atoms
        else:
            self.selection_a = universe.select_atoms(primary_selection)
            self.reference_selection_a = self.reference_universe.select_atoms(primary_selection)

        # If the selection group B is not provided, use selection group A
        if secondary_selection is None:
            self.selection_b = self.selection_a
            self.reference_selection_b = self.reference_selection_a
        else:
            self.selection_b = universe.select_atoms(secondary_selection)
            self.reference_selection_b = self.reference_universe.select_atoms(secondary_selection)

        # Get dimension of box if pbc set to True
        self.use_pbc = use_pbc
        if self.use_pbc:
            self._get_box = lambda ts: ts.dimensions
        else:
            self._get_box = lambda ts: None

        r0 = distance_array(self.reference_selection_a.positions, self.reference_selection_b.positions, box=self._get_box(self.reference_universe.universe))
        seq_sep = np.abs(self.reference_selection_a.resids[:, np.newaxis] - self.reference_selection_b.resids[np.newaxis, :], dtype=float)
        seq_sep [self.reference_selection_a.chainIDs[:, np.newaxis] != self.reference_selection_b.chainIDs[np.newaxis, :]] = np.inf
        
        #Select the indices
        self.i,self.j = np.where((r0<contact_cutoff) & (seq_sep >= min_seq_sep) & ~np.isnan(seq_sep))
        
        self.r0 = r0[self.i, self.j]
        self.seq_sep = seq_sep[self.i, self.j]
        # print(len(self.r0))
        # print(len(self.seq_sep))
        
    def _prepare(self):
        self.timeseries = np.empty((self.n_frames, len(self.r0)))
        # print(self.timeseries.shape)

    def _single_frame(self):
        self.timeseries[self._frame_index][0] = self._ts.frame
        # print(self.timeseries)
        
        # compute distance array for a frame
        r = distance_array(self.selection_a.positions, self.selection_b.positions, box=self._get_box(self._ts))
        q=self.method_function(r[self.i, self.j], self.r0, self.seq_sep)
        self.timeseries[self._frame_index] = q

    def _conclude(self):
        self.qvalues = np.mean(self.timeseries, axis=1)
        return super()._conclude()

def qvalue_wolynes(rij, rijn, seq_sep, sigma_exp=0.15):
    sigmaij = np.power(seq_sep,sigma_exp)
    return np.exp(-(rij - rijn)**2 / (2.0 * sigmaij**2))