:orphan:

=============
API Reference
=============


.. contents:: Table of Contents
   :depth: 1
   :local:
   :backlinks: none

:mod:`cocoatree`: Overall functions
===================================

.. automodule:: cocoatree
    :no-members:
    :no-inherited-members:

Functions
---------

.. currentmodule:: cocoatree
 

.. autosummary::

  :toctree: generated/
  :template: function.rst

  perform_sca


:mod:`cocoatree.io`: Functions to import and export files
===================================================================

.. automodule:: cocoatree.io
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree
 

.. autosummary::

  :toctree: generated/
  :template: function.rst

  io.load_MSA
  io.load_tree_ete3
  io.export_fasta
  io.load_pdb
  io.export_xcor_for_pymol

.. _io_ref:



:mod:`cocoatree.datasets`: Utilities to load popular datasets.
==============================================================

.. automodule:: cocoatree.datasets
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

   :toctree: generated/
   :template: function.rst

   datasets.load_S1A_serine_proteases
   datasets.load_rhomboid_proteases
   datasets.load_DHFR

.. _datasets_ref:


:mod:`cocoatree.msa`: Manipulating multiple sequence alignment and sequences
============================================================================

.. automodule:: cocoatree.msa
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

   :toctree: generated/
   :template: function.rst

   msa.filter_sequences
   msa.filter_seq_id
   msa.filter_ref_seq
   msa.map_to_pdb
   msa.compute_seq_identity
   msa.compute_seq_weights
   msa.map_msa_positions
   msa.compute_seq_similarity
   msa.compute_normalized_seq_similarity


.. _msa_ref:


:mod:`cocoatree.statistics`: Computation of position-specific and pairwise statistics
=====================================================================================

.. automodule:: cocoatree.statistics
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

   :toctree: generated/
   :template: function.rst

   statistics.compute_all_frequencies
   statistics.pairwise.compute_sca_matrix
   statistics.pairwise.compute_mutual_information_matrix
   statistics.pairwise.compute_apc
   statistics.pairwise.compute_entropy_correction
   statistics.position.compute_entropy
   statistics.position.compute_conservation


.. _statistics_ref:


:mod:`cocoatree.deconvolution`: Co-evolution matrix deconvolution
=================================================================

.. automodule:: cocoatree.deconvolution
  :no-members:
  :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::
    
    :toctree: generated/
    :template: function.rst

    deconvolution.extract_independent_components
    deconvolution.extract_principal_components
    deconvolution.extract_xcors_from_ICs
    deconvolution.extract_xcors
    deconvolution.remove_global_correlations


.. _deconvolution_ref:


:mod:`cocoatree.visualization`: Visualization with ete3
==========================================================

.. automodule:: cocoatree.visualization
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

    :toctree: generated/
    :template: function.rst

    visualization.update_tree_ete3_and_return_style
    visualization.add_xcor_sequences_to_tree
    visualization.add_heatmap_to_tree


.. _visualization_ref:

