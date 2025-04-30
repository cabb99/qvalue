Contact Q-value Analysis for Biomolecular Simulations
=======================================================================

The q-value is the fraction of native residue–residue contacts retained when compared to a reference structure. The package includes some of the most common Q-value definitions:

- **Wolynes Q**  
  A smooth, distance-weighted score: every reference contact contributes via a Gaussian centered at its native distance, with width growing as a power of sequence separation.  

- **Onuchic Q**  
  Similar Gaussian form as Wolynes, but only considers contacts beyond a minimum sequence gap and within a hard distance cutoff. 

- **Interface Q**  
  Computes Wolynes-style Q only for pairs that lie in different chains or domains—ideal for monitoring binding or assembly.  

- **Intrachain Q**  
  Restricts Wolynes Q to contacts within each chain separately—handy for multi-domain proteins or polymers.  

- **Contact Q (Qc)**  
  A version of Wolynes Q where only interactions between residues in contact in the native structure are considered based on a distance cutoffs

- **Custom Q**  
  Other definitions can be implemented by providing a callable function.

Each method can be applied to Cα, Cβ, or any atom selection, with optional per-residue or per-contact outputs, and integrates with MDAnalysis.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   api
   dev

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
