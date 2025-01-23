Developer Notes for Q-Value Analysis
====================================

Overview
--------

**Q-value** is a common measure in protein folding and biomolecular simulations. It quantifies how many native contacts are preserved (relative to a reference/native structure) at each simulation frame. In this implementation, **native contacts** are defined based on residue pair distances in a *reference structure*. Various “flavors” of Q-value calculations are provided, including:

1. **Wolynes Q**
2. **Onuchic Q**
3. **Interface Q**
4. **Intrachain Q**
5. **Contact** or “Qc”
6. **Custom** (user-defined) Q functions

The code utilizes `MDAnalysis <https://www.mdanalysis.org/>`_ to handle trajectory data (``Universe`` objects), inter-atomic distances, and selections. The main analysis is performed by the ``qValue`` class, which extends ``AnalysisBase`` from MDAnalysis.

Code Organization
-----------------

``qvalue.py``
~~~~~~~~~~~~~

Below is a high-level look at the file structure:

1. **``class qValue(AnalysisBase)``**  
   - Orchestrates Q-value computations across a trajectory.  
   - Handles multiple “methods” or flavors of Q calculations in one pass.  
   - Stores results in ``self.results[method_name]``.
   - Introduces a flexible “add_method” that allows both developers or end-users to supply custom selection logic or Q-value functional forms.

2. **Helper Functions**  
   - ``qvalue_pair_wolynes``  
   - ``qvalue_pair_onuchic``  
   - ``qvalue_pair_interface_CB``

   These are lower-level functions that compute the pairwise contribution to Q for each contact. Some of the “flavor” differences lie in how the Gaussian exponent is computed, which can depend on sequence separation, distance cutoffs, but also on the equation to calculate the q per pair.

qvalue.add_method()
-------------------

This function defines **which contacts** to track and **how** to score them.

**Key parameters**:

- **``method``**: ``'Wolynes'``, ``'Onuchic'``, ``'Interface_CB'``, ``'Intrachain'``, ``'Qc'``, or any user-defined callable function.
- **``method_name``**: Defines a method-specific key in the results dictionary.
- **``selection`` / ``complementary_selection``**: Which atoms (e.g., which chain or residue subset) to consider. By default, ``'all'``. Using two groups will consider only the pairs between them.
- **``contact_cutoff``**: Distance (Å) below which pairs are considered “native.” For ``Onuchic``, defaults to ~9.5 Å; for ``Wolynes``, defaults to ``np.inf``.
- **``min_seq_sep``** and **``max_seq_sep``**: Restrict pairs by sequence separation, if needed (e.g., ignoring local contacts).
- **``store_per_contact``**: Store a Q-value for every contact in every frame (useful for visualizing contact-specific Q-values in heatmaps).
- **``store_per_residue``**: Store a Q-value for each residue in every frame (by averaging over that residue’s contacts). Useful for coloring residues in a visualization as local Q-values.
- **``alternative_reference``**: An optional different reference structure (single frame) against which to calculate Q-values.
- **``extra_kwargs``**: Additional parameters for the Q-value function, e.g., sigma exponent.

After each call to ``add_method``, an internal record of the method’s parameters is stored (re-indexed for faster distance lookups).  
*(Note: multiple reference frames are not yet supported, but multiple universes can be used to define a reference.)*

Execution
~~~~~~~~~

.. code-block:: python

   q_calc = qValue(universe, reference_universe=ref)
   # Optionally add more specialized methods:
   # Each method's parameters can be customized
   q_calc.add_method('Onuchic', method_name="Custom_Onuchic",
                     min_seq_sep=4, contact_cutoff=9.5)
   q_calc.run()

Internally, MDAnalysis will iterate over the trajectory and call the following methods:

1. **``_prepare()``**: Collates all selected pairs across all methods, pre-computes reference distances and sequence separations, then merges them if needed.
2. **``_single_frame()``**: For each frame:
   - Computes distances for *all* pairs in a single pass (``distance_array``).
   - For each method, slices out the relevant i–j pairs and evaluates the chosen Q formula.
3. **``_conclude()``**: Final post-processing steps, if any.

Results
~~~~~~~

- ``q_calc.results[method_name]['q_values']``: NumPy array of shape ``(n_frames,)`` with the average Q-value per frame.
- Optionally, if ``store_per_contact = True`` or ``store_per_residue = True``, the arrays are stored as ``q_per_contact`` or ``q_per_residue`` in the same dictionary.

Helper q-value Functions
------------------------

1. **``qvalue_pair_wolynes``**  

   .. math::

      \sigma_{ij} = a \times |i-j|^{\text{sigma_exp}}

   .. math::

      q_{ij}(t) = \exp\Bigl(- \frac{\bigl(r_{ij}(t) - r_{ij}^N\bigr)^2}{2\,\sigma_{ij}^2}\Bigr)

2. **``qvalue_pair_onuchic``**  

   Similar to Wolynes, but includes a +1 in sequence separation:

   .. math::

      \sigma_{ij} = a \times \bigl(1 + |i-j|\bigr)^{\text{sigma_exp}}

3. **``qvalue_pair_interface_CB``**  

   Used for computing Q across chain interfaces using CB selections. Typically uses a fixed :math:`\sigma` based on half the total number of CB atoms:

   .. math::

      \sigma = \bigl(1 + (N/2)\bigr)^{\text{sigma_exp}}

   where ``N`` is the total count of :math:`C_B` atoms in the reference.

Developers can create additional variants by providing their own callable with the signature:

.. code-block:: python

   def my_qfunction(rij, rijn, seq_sep, **kwargs):
       # rij    - array of shape (n_contacts,) for current distances
       # rijn   - array of shape (n_contacts,) for reference distances
       # seq_sep- array of shape (n_contacts,) for sequence separations
       # kwargs - any user-defined parameters

       return q_per_contact  # shape (n_contacts,)

Example Usage
-------------

.. code-block:: python

   import MDAnalysis as mda
   from qvalue_module import qValue

   # Load trajectory + reference
   u = mda.Universe("native.pdb", "trajectory.dcd")
   ref = mda.Universe("native.pdb")  # single-frame

   # Instantiate qValue analysis
   q_calc = qValue(u, reference_universe=ref)

   # Add some methods
   q_calc.add_method(qvalue_pair_wolynes, method_name="Wolynes_CB",
                     selection='segid A',
                     complementary_selection='segid A',
                     atoms='name CB',
                     min_seq_sep=3, contact_cutoff=12.0)

   q_calc.add_method(qvalue_pair_onuchic, method_name="Onuchic_CA",
                     selection='name CA and segid A',
                     complementary_selection='name CA and segid A',
                     min_seq_sep=10, contact_cutoff=9.5)

   q_calc.add_method('Interface_CB', method_name="Interface_AB",
                     selection='segid A', complementary_selection='segid B',
                     # Extra kwargs for the custom interface function
                     N=200, sigma_exp=0.2, store_per_contact=True)

   # Run
   q_calc.run()

   # Extract results
   q_wolynes = q_calc.results["Wolynes_CA"]["q_values"]
   q_onuchic = q_calc.results["Onuchic_CA"]["q_values"]
   q_iface   = q_calc.results["Interface_AB"]["q_values"]

   # Optionally, if store_per_contact was True:
   # q_per_contact_iface = q_calc.results["Interface_AB"]["q_per_contact"]

Testing & Validation
--------------------

- **``test_qvalue.py``**:  
  We validate that Q-values match known reference data (e.g., from openAWSEM). If you add a new Q method, add corresponding tests to ensure consistent behavior.

- **Continuous Integration**:  
  The repository includes GitHub Actions for automated testing. Any push or pull request triggers the test suite, verifying correctness and preventing regression.

Future Directions
-----------------

- **mmCIF support**:  
  Pending updates from MDAnalysis will enable cif compatibility. In the meantime, you can use the custom `MDAnalysis_CIFReading plugin <https://github.com/cabb99/MDAnalysis_CIFReading>`_.

- **Multiple frames**:  
  Currently, only one reference frame is supported. We plan to extend this to multiple reference frames, which can be useful for ensemble-based Q-value calculations.

- **Performance optimizations**:
  - **Parallelization**: Distribute the Q-value computation across multiple cores.
  - **GPU acceleration**: GPU-accelerated distance calculations for ensemble systems using numba or similar libraries.

- **Visualization**:
  - **Interactive plots**: Use Bokeh or Plotly to create interactive plots for Q-value trajectories.
  - **Heatmaps**: Visualize contact-specific Q-values as heatmaps.

- **User Interface**:
  - **Command-line interface**: Create a CLI for running Q-value analysis on trajectories.

- **Documentation**:
  - **Tutorials**: Write tutorials for common Q-value analyses.
  - **API Reference**: Document the ``qValue`` class and helper functions.

- **Integration with other tools**:
  - **OpenAWSEM / AWSEMtools**: Integrate Q-value analysis with OpenAWSEM for easy calculations.

----

*Last updated: 2025-01-22*