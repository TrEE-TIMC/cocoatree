:orphan:

=========
Reference
=========


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

.. _io_ref:



:mod:`cocoatree.datasets`: Datasets
===================================

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

.. _datasets_ref:


:mod:`cocoatree.msa`: Multiple sequence alignment
=================================================

.. automodule:: cocoatree.msa
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

   :toctree: generated/
   :template: function.rst

   msa.filter_gap_pos
   msa.filter_gap_seq
   msa.filter_ref_seq
   msa.filter_seq_id
   msa.choose_ref_seq
   msa.compute_seq_identity


.. _msa_ref:


:mod:`cocoatree.statistics`: Computation of statistics
======================================================

.. automodule:: cocoatree.statistics
    :no-members:
    :no-inherited-members:


Functions
---------

.. currentmodule:: cocoatree

.. autosummary::

   :toctree: generated/
   :template: function.rst

   statistics.pairwise.aa_joint_freq
   statistics.pairwise.compute_seq_identity
   statistics.pairwise.compute_sca_matrix
   statistics.pairwise.compute_mutual_information_matrix
   statistics.pairwise.compute_apc
   statistics.pairwise.compute_entropy_correction
   statistics.position.aa_freq_at_pos
   statistics.position.compute_background_frequencies
   statistics.position.compute_entropy
   statistics.position.compute_rel_entropy
   statistics.position.aa_freq_at_pos
   statistics.position.compute_background_frequencies
   statistics.position.compute_rel_entropy
   statistics.sequence.compute_seq_weights


.. _statistics_ref:


:mod:`cocoatree.deconvolution`: Matrix deconvolution
=====================================================

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

  deconvolution.eigen_decomp
  deconvolution.compute_ica

    randomize.randomization

.. _randomize_ref:


:mod:`cocoatree.visualization`: Visualization
====================================================

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

