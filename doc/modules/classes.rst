:orphan:

=============
API Reference
=============


.. contents:: Table of Contents
   :depth: 1
   :local:
   :backlinks: none


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
  io.export_sector_for_pymol

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
   msa.map_to_pdb
   msa.compute_seq_identity
   msa.compute_seq_weights


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

    deconvolution.eigen_decomp
    deconvolution.compute_ica
    deconvolution.choose_num_components
    deconvolution.extract_positions_from_IC
    deconvolution.eigen_decomp
    deconvolution.compute_ica

.. _deconvolution_ref:


:mod:`cocoatree.randomize`: Alignment randomization
====================================================

.. automodule:: cocoatree.randomize
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

  :toctree: generated/
  :template: function.rst


  randomize.randomization

.. _randomize_ref:


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


.. _visualization_ref:

