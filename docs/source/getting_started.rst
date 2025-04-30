Getting Started
===============

Q-value is a metric used in protein folding and biomolecular simulations that quantifies the fraction of native contacts preserved relative to a reference structure at each simulation frame.
The qvalue package provides a flexible and efficient implementation for computing q-values, leveraging MDAnalysis for trajectory handling and atom selections. 
Multiple Q-value "flavors" are supported, including Wolynes, Onuchic, Interface, Intrachain, and Contact (Qc) methods.

Features
--------
- **Multiple Q-value methods**: Wolynes, Onuchic, Interface, Intrachain, Contact (Qc), and user-defined custom Q functions.
- **Flexible atom/group selections**: Analyze Q-values for specific chains, selections, or interfaces.
- **Per-contact and per-residue Q-values**: Store and visualize Q-values for individual contacts or residues.
- **Customizable parameters**: Control distance cutoffs, sequence separation, reference structures, and more.
- **Integration with MDAnalysis**: Use familiar Universe objects and selection syntax.
- **Extensible**: Add new Q-value methods by supplying your own Python function.

Installation
------------

To install, use pip:

.. code-block:: bash

    pip install md-qvalue

Or for development:

.. code-block:: bash

    git clone https://github.com/your-repo/md-qvalue.git
    cd md-qvalue
    pip install -e .

Once available on conda-forge:

.. code-block:: bash

    conda install -c conda-forge md-qvalue

For developers, see the `Developer Notes <dev.html>`_.

Basic Usage
-----------

The Q-value for a contact is computed as:

.. math::
    q_{ij} = \exp\left(-\frac{(r_{ij} - r_{ij}^N)^2}{2\sigma_{ij}^2}\right)

where:
    - :math:`r_{ij}` is the distance between atoms :math:`i` and :math:`j` in the current frame
    - :math:`r_{ij}^N` is the reference (native) distance
    - :math:`\sigma_{ij}` is a the expected fluctuation of the contact that generally depends on the sequence separation

The total Q-value is the average over all native contacts:

.. math::

    Q = \frac{1}{N_{\text{contacts}}} \sum_{(i,j) \in \text{native}} q_{ij}

where :math:`N` is the number of CB atoms in the system.

### Example: Q-value Calculation for a Trajectory
---------------------------------------------

.. code-block:: python

    import MDAnalysis as mda
    from md_qvalue import qValue, qvalue_pair_wolynes, qvalue_pair_onuchic

    u = mda.Universe("native.pdb", "trajectory.dcd")
    ref = mda.Universe("native.pdb")
    q_calc = qValue(u, reference_universe=ref)

    # Add Wolynes Q for CB atoms in chain A
    q_calc.add_method(qvalue_pair_wolynes, method_name="Wolynes_CB",
                     selection='segid A', complementary_selection='segid A',
                     atoms='name CB', min_seq_sep=3, contact_cutoff=12.0)

    # Add Onuchic Q for CA atoms in chain A
    q_calc.add_method(qvalue_pair_onuchic, method_name="Onuchic_CA",
                     selection='name CA and segid A', complementary_selection='name CA and segid A',
                     min_seq_sep=10, contact_cutoff=9.5)

    # Add interface Q between chains A and B (CB atoms)
    q_calc.add_method('Interface_CB', method_name="Interface_AB",
                     selection='segid A', complementary_selection='segid B',
                     N=200, sigma_exp=0.2, store_per_contact=True)

    q_calc.run()
    q_wolynes = q_calc.results["Wolynes_CB"]["q_values"]
    q_onuchic = q_calc.results["Onuchic_CA"]["q_values"]
    q_iface   = q_calc.results["Interface_AB"]["q_values"]

    # If store_per_contact was True:
    # q_per_contact_iface = q_calc.results["Interface_AB"]["q_per_contact"]

Custom Q-value Methods
----------------------

You can define your own Q-value function and use it with add_method:

.. code-block:: python

    def my_qfunction(rij, rijn, seq_sep, **kwargs):
        # Custom logic here
        return np.exp(-((rij - rijn)**2) / (2.0 * (0.5 + seq_sep)**2))

    q_calc.add_method(my_qfunction, method_name="CustomQ", ...)

Advanced Features
-----------------
- **Flexible selections**: Use any valid MDAnalysis selection string for atoms or groups.
- **Multiple methods**: Add as many Q-value methods as you want in a single analysis.
- **Per-contact and per-residue storage**: Enable store_per_contact or store_per_residue for detailed output.
- **Alternative references**: Use a different structure as the reference for any method.

Testing
-------

To ensure everything is working correctly, you can run the test suite. First, install the test dependencies:

.. code-block:: bash

    pip install -r devtools/conda-envs/test_env.yaml

Then, run the tests using pytest:

.. code-block:: bash

    pytest -v .

References
----------

The q-value methodologies have been described in the following papers:

- Wei Lu, Carlos Bueno, Nicholas P. Schafer, Joshua Moller, Shikai Jin, Xun Chen, Mingchen Chen, Xinyu Gu, Aram Davtyan, Juan J. de Pablo, Peter G. Wolynes. "OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations." *PLoS Computational Biology*, 17(2): e1008308, 2021. https://doi.org/10.1371/journal.pcbi.1008308
- Samuel S. Cho, Yaakov Levy, Peter G. Wolynes. "P versus Q: Structural reaction coordinates capture protein folding on smooth landscapes." *PNAS*, 103(3): 586-591, 2006. https://doi.org/10.1073/pnas.0509768103
- Bobby L. Kim, Nicholas P. Schafer, Peter G. Wolynes. "Predictive energy landscapes for folding α-helical transmembrane proteins." *PNAS*, 111(30): 11031-11036, 2014. https://doi.org/10.1073/pnas.1410529111
- Corey Hardin, Michael P. Eastwood, Michael C. Prentiss, Peter G. Wolynes. "Associative memory Hamiltonians for structure prediction without homology: α/β proteins." *PNAS*, 100(4): 1679-1684, 2003. https://doi.org/10.1073/pnas.252753899
