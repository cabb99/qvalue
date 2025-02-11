"""
Q-Value Analysis
=================

The Q-value is a measure of the fraction of native contacts conserved during a simulation. 
Native contacts are defined from the selected atoms in a reference structure. 
This module implements two primary Q-value calculation methods:

1. **Wolynes Q**: A weighted measure of native contact preservation.
    - Native contacts are defined without applying a specific distance cutoff.
    - Each contact is weighted using a Gaussian function centered on the reference structure's inter-atomic distance and scaled by the sequence separation:

    .. math::
        
        q(t) = \\frac{1}{N} \\sum_{i,j} \\exp\\left(-\\frac{(r_{ij} - r_{ij}^N)^2}{2\\sigma_{ij}^2}\\right)

    where:
    - :math:`r_{ij}` is the distance between atoms :math:`i` and :math:`j` at time :math:`t`.
    - :math:`r_{ij}^N` is the reference distance between atoms :math:`i` and :math:`j`.
    - :math:`\sigma_{ij}` is the separation in the sequence between residues :math:`i` and :math:`j`.

2. **Onuchic Q**: Similar to Wolynes Q but applies additional constraints:
    - Excludes contacts with a sequence separation below a threshold.
    - Incorporates a distance cutoff to restrict the evaluation to native contacts.
    - :math:`\sigma_{ij}` is the separation in the sequence between residues :math:`i` and :math:`j` plus 1.

Usage Example
-------------

.. code-block:: python

    from MDAnalysis import Universe
    from qvalue_module import qValue
    u = Universe("native.pdb", "trajectory.dcd")
    q_calc = qValue(u)
    q_calc.run()
    q_wolynes = q_calc.results['Wolynes']['q_values']
    q_onuchic = q_calc.results['Onuchic']['q_values']
"""


from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import distance_array
import numpy as np

class qValue(AnalysisBase):
    """
    Calculates Q-values (fraction of native contacts) across a trajectory by comparing inter-atomic distances in a trajectory to those in a reference structure.

    
    Parameters
    ----------
    universe :  MDAnalysis.core.universe.Universe 
        The universe (trajectory) on which Q-values are calculated.
    reference_universe :  MDAnalysis.core.universe.Universe , optional
        A separate reference structure/universe. If not provided, the first frame of  universe  is used as the reference.
    primary_selection :  str  or  MDAnalysis.core.groups.AtomGroup , optional
        A selection string (e.g.,  'name CA' ) or an AtomGroup for the primary group of atoms.
    secondary_selection :  str  or  MDAnalysis.core.groups.AtomGroup , optional
        A selection string or an AtomGroup for the secondary group of atoms.
    q_method: {'Wolynes','Onuchic', 'Interface','Interchain', callable}, optional
        Method to compute the Q-values. Default is "Wolynes". If callable, it must accept distances and return Q-values.
    use_pbc : bool, optional
        Whether to apply periodic boundary conditions. Default is  False .
    contact_cutoff : float, optional
        Maximum distance (in Å) below which two atoms are considered in contact. Default is  np.inf .
    method_kwargs : dict, optional
        Additional arguments for the Q-value method.
    **basekwargs : dict
        Additional keyword arguments passed to the base  AnalysisBase  class.

    Attributes
    ----------
    qvalue_function :  callable 
        Function used to compute Q-values for each frame.
    r0 :  np.ndarray 
        Reference inter-atomic distances.
    seq_sep :  np.ndarray 
        Sequence separation for atom pairs.
    qvalues :  np.ndarray 
        Computed Q-values for each frame.

    See Also
    --------
    MDAnalysis.analysis.base.AnalysisBase : The base class for MDAnalysis analyses.
    """
        
    def __init__(self, 
                 universe, 
                 reference_universe=None, 
                 method=["Wolynes","Onuchic"], 
                 use_pbc=False, 
                 custom_CA_selection = None,
                 custom_CB_selection = None,
                 **basekwargs):
        """
        Initialize the Q-value calculation.
        
        Parameters
        ----------
        universe : MDAnalysis.core.universe.Universe
            The universe (trajectory) on which Q-values are calculated.
        reference_universe : MDAnalysis.core.universe.Universe, optional
            A separate reference structure/universe. If not provided, the first frame of universe is used as the reference.
        method : {'Wolynes','Onuchic', 'Interface','Intrachain', callable}, optional
            Method to compute the Q-values. Default is "Wolynes". If callable, it must accept distances and return Q-values.
        use_pbc : bool, optional
            Whether to apply periodic boundary conditions. Default is False.
        custom_CA_selection : str, optional
            Custom selection string for the CA atoms. Default is 'name CA'.
        custom_CB_selection : str, optional
            Custom selection string for the CB atoms. Default is 'name CB or (name CA and resname GLY IGL)'.
        **basekwargs : dict
            Additional keyword arguments passed to the method that calculates the q value.
        
            """
        self.universe = universe
        super(qValue, self).__init__(self.universe.trajectory, **basekwargs)

        # Get dimension of box if pbc set to True
        self.use_pbc = use_pbc
        if self.use_pbc:
            self._get_box = lambda ts: ts.dimensions
        else:
            self._get_box = lambda ts: None

        # If the reference group is not provided, use the first frame of the trajectory
        if reference_universe is None:
            self.reference_universe = universe
        else:
            self.reference_universe = reference_universe

        self.reference_selection_CA = self.select_CA(self.reference_universe, custom_CA_selection)
        self.selection_CA = self.select_CA(self.universe, custom_CA_selection)
        self.reference_selection_CB = self.select_CB(self.reference_universe, custom_CB_selection)
        self.selection_CB = self.select_CB(self.universe, custom_CB_selection)
        
        # Check that the selections are the same length
        self.assert_same_length(self.reference_selection_CA, self.selection_CA, "CA selection length mismatch between reference and trajectory")
        self.assert_same_length(self.reference_selection_CB, self.selection_CB, "CB selection length mismatch between reference and trajectory")

        self.nframes = len(self.universe.trajectory)
        
        self.methods=[]
        #If the method is an iterable add each method, else add the single method
        if hasattr(method, '__iter__'):
            for method in method:
                self.add_method(method)
        else:
            self.add_method(method)

        
    def select_CA(self, universe, custom_CA_selection=None):
        """
        Select the alpha carbon atoms from the universe.
        """
        #CA and CB selection
        if custom_CA_selection is not None:
            selection_CA = universe.select_atoms(custom_CA_selection)
        else:
            try:
                selection_CA = universe.select_atoms('name CA')
            except AttributeError as e:
                print("CA selection failed for trajectory. Assuming lammps input.")
                selection_CA = universe.select_atoms('type 1')
        return selection_CA

    def select_CB(self, universe, custom_CB_selection):
        """
        Select the beta carbon atoms from the universe.
        """
        if custom_CB_selection is not None:
            selection_CB = universe.select_atoms(custom_CB_selection)
        else:
            try:
                selection_CB = universe.select_atoms('name CB or (name CA and resname GLY IGL)')
            except AttributeError as e:
                print("CB selection failed for trajectory. Assuming lammps input.")
                selected_indices = []
                for atom in universe.atoms:
                    if atom.type == '1':
                        CA_atom=atom
                    elif atom.type == '4':
                        selected_indices.append(atom.index)
                    elif atom.type == '5':
                        selected_indices.append(CA_atom.index)
                selection_CB = universe.atoms[selected_indices]
        return selection_CB

    def assert_same_length(self, selection_a, selection_b, message=None):
        """
        Check that two selections have the same length.
        """
        if len(selection_a) != len(selection_b):
            # If the selections are not the same length go through the selections and find the residues that are not in the other selection
            resids_a = selection_a.atoms.resids
            resids_b = selection_b.atoms.resids
            if len(resids_a) > len(resids_b):
                missing_residues = [resid for resid in resids_a if resid not in resids_b]
            else:
                missing_residues = [resid for resid in resids_b if resid not in resids_a]

            print(message)
            print(f"Selections have different lengths.")
            print(f"Length of selection A: {len(selection_a)}")
            print(f"Length of selection B: {len(selection_b)}")
            print(f"Missing residues: {missing_residues}")
            
            raise ValueError("Selections must have the same length.")


    def _updated_method_name(self, method_name):
        """
        Generate a unique method name by adding a suffix
        if the given method_name already exists in self.methods.

        Parameters:
            method_name (str): The base name of the method to check.

        Returns:
            str: A unique method name.

        Raises:
            ValueError: If a unique name cannot be generated.
        """

        if any(m['name'] == method_name for m in self.methods):
            i = 1
            while f"{method_name}_{i}" in self.methods:
                i += 1
                if i >= 1000:
                    raise ValueError(f"Unable to generate a unique name for '{method_name}'.")
            return f"{method_name}_{i}"
        return method_name
    
    def _force_selection(self, reference_atoms, atoms, selection_str):
        ref_indices = reference_atoms.select_atoms(selection_str).indices
        indices = atoms.select_atoms(selection_str).indices

        if len(ref_indices) == 0:
            raise ValueError(f"No atoms selected in the reference for selection '{selection_str}'.")

        #Select the indices that are in the reference indices
        if len(ref_indices) != len(indices):
            indices = np.array([a for a,b in zip(atoms.indices, reference_atoms.indices) if b in ref_indices])
        
        if len(ref_indices) != len(indices):
            raise ValueError(f"Selection length mismatch between reference and trajectory for selection '{selection_str}'.")
        
        return ref_indices, indices


    def add_method(self, method, method_name=None, 
                   selection=None, complementary_selection=None, atoms=None,
                   contact_cutoff=None, min_seq_sep=None, max_seq_sep=None,
                   store_per_residue=False, store_per_contact=False,
                   alternative_reference = None,
                   **kwargs):
        """
        Add a new method to the Q-value calculation.

        Parameters
        ----------
        method : {'Wolynes','Onuchic', 'Interface','Intrachain', callable}
            Method to compute the Q-values. If callable, it must accept distances and return Q-values.
        method_name : str, optional
            Name of the method. If not provided, a unique name is generated.
        selection : str, optional
            Selection string for the primary group of atoms. Default is 'all'.
        complementary_selection : str, optional
            Selection string for the secondary group of atoms. Default is 'all'.
        atoms : str, optional
            Atom selection string. Default is 'name CA'.
        contact_cutoff : float, optional
            Maximum distance (in Å) below which two atoms are considered in contact. Default is np.inf.
        min_seq_sep : int, optional
            Minimum sequence separation between residues. Default is 3.
        max_seq_sep : int, optional
            Maximum sequence separation between residues. Default is np.inf.
        store_per_residue : bool, optional
            Whether to store Q-values per residue. Default is False.
        store_per_contact : bool, optional
            Whether to store Q-values per contact. Default is False.
        alternative_reference : MDAnalysis.core.universe.Universe, optional
            An alternative reference structure/universe. If not provided, the reference universe is used.
        **kwargs : dict
        """
        
        method_description = {'name':method_name, 
                              'function':method, 
                              'selection':selection,
                              'complementary_selection':complementary_selection,
                              'atoms':atoms,
                              'contact_cutoff':contact_cutoff, 
                              'min_seq_sep':min_seq_sep, 
                              'max_seq_sep':max_seq_sep,
                              'reference':alternative_reference,
                              'store_per_contact':store_per_contact,
                              'store_per_residue':store_per_residue,
                              'kwargs': kwargs}
        
        #General defaults
        if selection is None:
            method_description['selection'] = 'all'
        if complementary_selection is None:
            method_description['complementary_selection'] = 'all'
        if contact_cutoff is None:
            method_description['contact_cutoff'] = np.inf
        if atoms is None:
            method_description['atoms'] = 'name CA'
        if min_seq_sep is None:
            method_description['min_seq_sep'] = 3
        if max_seq_sep is None:
            method_description['max_seq_sep'] = np.inf
        if alternative_reference is None:
            method_description['reference'] = self.reference_universe
        
            
        #Generate a unique method name
        if method_name is not None:
            method_description['name'] = self._updated_method_name(method_name)
        elif type(method) is str:
            method_description['name'] = self._updated_method_name(method)
        else:
            method_description['name'] = self._updated_method_name('Custom')
        reference = method_description['reference']

        if len(reference.trajectory) > 1:
            raise ValueError("Alternative reference must be a single frame.")
                  
        # Define the method function based on the method name
        if method == 'Wolynes':
            method_description['function'] = qvalue_pair_wolynes
        elif method == 'Onuchic':
            if contact_cutoff is None:
                method_description['contact_cutoff'] = 9.5
            if min_seq_sep is None:
                method_description['min_seq_sep'] = 4
            method_description['function'] = qvalue_pair_onuchic
        elif method == 'Interface':
            raise NotImplementedError()
            # chainIDs = np.unique(self.reference_universe.atoms.segids)
            # for chain_a in chainIDs:
            #     for chain_b in chainIDs:
            #         if chain_a != chain_b:
            #             selection = f"segid {chain_a}"
            #             complementary_selection = f"segid {chain_b}"
            #             self.add_method('Wolynes', method_name=f'Interface_{chain_a}_{chain_b}', selection=selection, complementary_selection=complementary_selection, atoms='name CA', **kwargs)
            # return
        elif method == 'Interface_CB':
            method_description['min_seq_sep'] = 0 if min_seq_sep is None else min_seq_sep
            method_description['contact_cutoff'] = 9.5 if contact_cutoff is None else contact_cutoff
            if 'N' not in kwargs:
                method_description['kwargs']['N'] = len(self.reference_selection_CB)

            chainIDs = np.unique(self.reference_universe.atoms.segids)
            for chain_a in chainIDs:
                for chain_b in chainIDs:
                    if chain_a != chain_b:
                        selection = f"segid {chain_a}"
                        complementary_selection = f"segid {chain_b}"
                        self.add_method(method = qvalue_pair_interface_CB, method_name=f'Interface_{chain_a}_{chain_b}', selection=selection, complementary_selection=complementary_selection, atoms='name CB', **method_description['kwargs'],
                                        contact_cutoff=method_description['contact_cutoff'], min_seq_sep=method_description['min_seq_sep'])
            return
        elif method == 'Intrachain':
            raise NotImplementedError()
            # chainIDs = np.unique(self.reference_universe.atoms.segids)
            # for chain in chainIDs:
            #     selection = f"segid {chain}"
            #     self.add_method('Wolynes', method_name=f'Intrachain_{chain}', selection=selection, complementary_selection=selection, atoms='name CA', **kwargs)
            # return
        elif method == 'Contact':
            method_description['min_seq_sep'] = 10 if min_seq_sep is None else min_seq_sep
            if 'sigma' not in kwargs:
                method_description['kwargs']['a'] = 2
            method_description['function'] = qvalue_pair_wolynes

        elif callable(method):
            pass
        else:
            raise ValueError(f"Method '{method}' is not recognized.")
        
        # Get the indices of the selected atoms
        if method_description['atoms']=='name CA' and alternative_reference is None:
            ref_indices_a, indices_a = self._force_selection(self.reference_selection_CA, self.selection_CA, method_description['selection'])
            ref_indices_b, indices_b = self._force_selection(self.reference_selection_CA, self.selection_CA, method_description['complementary_selection'])
        elif method_description['atoms']=='name CB' and alternative_reference is None:
            ref_indices_a, indices_a = self._force_selection(self.reference_selection_CB, self.selection_CB, method_description['selection'])
            ref_indices_b, indices_b = self._force_selection(self.reference_selection_CB, self.selection_CB, method_description['complementary_selection'])
        else:
            ref_indices_a, indices_a = self._force_selection(reference.select_atoms(method_description['atoms']), self.universe.select_atoms(method_description['atoms']), method_description['selection'])
            ref_indices_b, indices_b = self._force_selection(reference.select_atoms(method_description['atoms']), self.universe.select_atoms(method_description['atoms']), method_description['complementary_selection'])


        ref_indices_a_query = f'index {" ".join(map(str, ref_indices_a))}'
        ref_indices_b_query = f'index {" ".join(map(str, ref_indices_b))}'

        method_description.update({'ref_indices_a':ref_indices_a, 'ref_indices_b':ref_indices_b, 'indices_a':indices_a, 'indices_b':indices_b})


        if alternative_reference is None:
            alternative_reference = self.reference_universe

        #Calculate the reference distances and sequence separations
        r0 = distance_array(alternative_reference.select_atoms(ref_indices_a_query).positions, alternative_reference.select_atoms(ref_indices_b_query).positions, box=self._get_box(alternative_reference.universe))
        seq_sep = np.abs(alternative_reference.select_atoms(ref_indices_a_query).resids[:, np.newaxis] - alternative_reference.select_atoms(ref_indices_b_query).resids[np.newaxis, :], dtype=float)
        seq_sep[alternative_reference.select_atoms(ref_indices_a_query).atoms.segindices[:, np.newaxis] != alternative_reference.select_atoms(ref_indices_b_query).atoms.segindices[np.newaxis, :]] = len(seq_sep)

        assert r0.shape == seq_sep.shape, "Reference distances and sequence separations must have the same shape."
        assert len(ref_indices_a) == len(r0), "Reference distances must have the same length as the reference indices."

        # Filter the pairs based on the cutoffs
        max_r0=method_description['contact_cutoff']
        min_seq_sep=method_description['min_seq_sep']
        max_seq_sep=method_description['max_seq_sep']
        i,j = np.where((r0<max_r0) & (seq_sep >= min_seq_sep) & (seq_sep <= max_seq_sep))

        #Select the indices from the trajectory
        ti = indices_a[i]
        tj = indices_b[j]
        
        if method_description['store_per_residue']:
            ti_elements, ti_counts = np.unique(ti,return_counts=True)
            tj_elements, tj_counts = np.unique(tj,return_counts=True)
            tij_elements = np.unique(np.concatenate((ti_elements, tj_elements)))
            n_residues = len(tij_elements)
           
            ti_sort_idx = np.argsort(ti_elements)
            ti_elements = ti_elements[ti_sort_idx]
            ti_counts = ti_counts[ti_sort_idx]

            tj_sort_idx = np.argsort(tj_elements)
            tj_elements = tj_elements[tj_sort_idx]
            tj_counts = tj_counts[tj_sort_idx]

            tij_elements.sort()

            ti_prime = np.searchsorted(ti_elements, ti)
            tj_prime = np.searchsorted(tj_elements, tj)
            ti_elements_prime = np.searchsorted(tij_elements, ti_elements)
            tj_elements_prime = np.searchsorted(tij_elements, tj_elements)

            q_count = np.zeros(n_residues)

            q_count[ti_elements_prime] += ti_counts
            q_count[tj_elements_prime] += tj_counts

        #Select r0 and seq_sep for the pairs
        r0 = r0[i, j]
        seq_sep = seq_sep[i, j]
  
        method_description.update({'ti':ti, 'tj':tj})
        method_description.update({'n_pairs':len(i)})
        method_description.update({'r0':r0, 'seq_sep':seq_sep})
        if method_description['store_per_residue']:
            method_description['chain_res'] = [f'{atom.segid}_{atom.resid}' for atom in self.universe.atoms if atom.index in tij_elements]
            method_description.update({'ti_prime':ti_prime, 'tj_prime':tj_prime})
            #method_description.update({'ti_elements':ti_elements, 'tj_elements':tj_elements})
            method_description.update({'ti_elements_prime':ti_elements_prime, 'tj_elements_prime':tj_elements_prime})
            #method_description.update({'ti_counts':ti_counts, 'tj_counts':tj_counts})
            method_description.update({'q_count':q_count})
            method_description.update({'n_i':len(ti_elements), 'n_j':len(tj_elements)})
            method_description.update({'n_residues':n_residues})

        self.methods.append(method_description)

    def _prepare(self):
        #Define the selections for each method
        self.all_ti = []
        self.all_tj = []
        self.results = {}
        for method in self.methods:
            self.all_ti.append(method['ti'])
            self.all_tj.append(method['tj'])
            result_template = {'q_values':np.empty(self.nframes)}
            if method['store_per_contact']:
                result_template['q_per_contact'] = np.empty((self.n_frames, method['n_pairs']))
            if method['store_per_residue']:
                result_template['chain_res'] = method['chain_res']
                result_template['q_per_residue'] = np.empty((self.n_frames, method['n_residues']))
            self.results.update({method['name']:result_template})

        self.all_ti=np.unique(np.concatenate(self.all_ti))
        self.all_tj=np.unique(np.concatenate(self.all_tj))
        self.all_ti.sort()
        self.all_tj.sort()
        self.all_ti_query = f'index {" ".join(map(str, self.all_ti))}'
        self.all_tj_query = f'index {" ".join(map(str, self.all_tj))}'

        #Reindex ti
        for method in self.methods:
            method['ti_reindex'] = [self.all_ti.tolist().index(x) for x in method['ti']]
            method['tj_reindex'] = [self.all_tj.tolist().index(x) for x in method['tj']]
        
    def _single_frame(self):
        # distances in the current frame
        r_matrix = distance_array(
            self.universe.select_atoms(self.all_ti_query).positions,
            self.universe.select_atoms(self.all_tj_query).positions,
            box=self._get_box(self._ts)
        )
        # Pull out only the pairs we decided were “native” in the reference
        for method in self.methods:
            i=method['ti_reindex']
            j=method['tj_reindex']
            rij_frame = r_matrix[i, j]
            
            # Evaluate the pairwise exponent for each contact
            q_per_contact = method['function'](
                rij_frame,    # current distances (Å)
                method['r0'],      # reference distances (Å)
                method['seq_sep'], # sequence separations
                **method['kwargs']
            )
            
            if method['store_per_contact']:
                self.results[method['name']]['q_per_contact'][self._frame_index] = q_per_contact
            if method['store_per_residue']:
                ti_prime = method['ti_prime']
                tj_prime = method['tj_prime']
                ti_elements_prime = method['ti_elements_prime']
                tj_elements_prime = method['tj_elements_prime']
                n_residues = method['n_residues']
                n_i = method['n_i']
                n_j = method['n_j']
                q_count = method['q_count']
                
                q_matrix = np.zeros((n_i,n_j))
                q_matrix[ti_prime,tj_prime] = q_per_contact

                q_sum = np.zeros(n_residues)

                q_sum[ti_elements_prime] += q_matrix.sum(axis=1)
                q_sum[tj_elements_prime] += q_matrix.sum(axis=0)
                q_per_residue = q_sum / q_count

                self.results[method['name']]['q_per_residue'][self._frame_index] = q_per_residue
            self.results[method['name']]['q_values'][self._frame_index] = np.mean(q_per_contact)

    def _conclude(self):
        pass

    def plot(self):
        """
        Plot the Q-values for each method.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        for method in self.results:
            ax.plot(self.results[method]['q_values'], label=method)
        ax.legend()
        plt.show()
        return ax

def qvalue_pair_wolynes(rij, rijn, seq_sep, a = 1.0, sigma_exp = 0.15):
    """
    Wolynes Q value contribution for each pair of atoms.
    
    Parameters
    ----------
    rij : array
        Distances between the i-j pairs in the frame (Å).
    rijn : array
        Distances between the i-j pairs in the reference structure (Å).
    seq_sep : array
        Sequence separations between residues i and j.
    a : float
        Base “sigma” in Å for one residue separation, e.g. a=1.0 means 1 Å
        for a 1-residue separation (equivalent to 0.1 nm if you were in nm).
    sigma_exp : float
        Exponent for the sequence separation: sigma_ij = a * (|i-j|^sigma_exp).
    
    Returns
    -------
    q_per_pair : array
        An array of shape (n_contacts,) giving exp(-((rij-rijn)^2)/(2*sigma^2)).
    """
    sigma_ij = a*np.power(seq_sep, sigma_exp)
    return np.exp(-(rij - rijn)**2 / (2.0 * sigma_ij**2))

def qvalue_pair_onuchic(rij, rijn, seq_sep, a = 1.0, sigma_exp = 0.15):
    """
    Onuchic Q value contribution for each pair of atoms.
    
    Parameters
    ----------
    rij : array
        Distances between the i-j pairs in the frame (Å).
    rijn : array
        Distances between the i-j pairs in the reference structure (Å).
    seq_sep : array
        Sequence separations between residues i and j.
    a : float
        Base “sigma” in Å for one residue separation, e.g. a=1.0 means 1 Å
        for a 1-residue separation (equivalent to 0.1 nm if you were in nm).
    sigma_exp : float
        Exponent for the sequence separation: sigma_ij = a * (1+|i-j|)^sigma_exp.
        """
    sigma_ij = a*np.power(1+seq_sep, sigma_exp)
    return np.exp(-(rij - rijn)**2 / (2.0 * sigma_ij**2))

def qvalue_pair_interface_CB(rij, rijn, seq_sep, N=None, sigma_exp=0.15):
    # N must be the total number of CB atoms (the old script’s “len(cb_atoms_pdb)”)
    # Use the old script’s “(1 + (N//2))**sigma_exp” for *all* pairs
    fixed_sigma = (1 + (N//2))**sigma_exp
    return np.exp(-(rij - rijn)**2/(2.0 * fixed_sigma**2))