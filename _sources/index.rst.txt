.. cocoatree documentation master file, created by
   sphinx-quickstart on Thu Apr 29 13:51:21 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
  :hidden:
  :glob:
  :caption: Contents

  interface_cocoatree_pymol.rst

#########################################################################################
COCOA-Tree: COllaborative COevolution Analysis Toolbox including tree-based vizualisation
#########################################################################################

COCOA-Tree is a `Python <https://www.python.org>`_ library designed to perform coevolution analysis of protein sequence data and to compare the results with phylogenetic trees and metadata. By doing so, it facilitates the investigation of the origins of coevolution patterns, their relationships to phylogeny and functional protein annotation, and ultimately aims to democratize the use and development of so-called *protein sectors* :cite:p:`rivoire_evolution-based_2016` and *specificity-determining positions (SDPs)* :cite:p:`de_juan_emerging_2013`.

The library includes several coevolution metrics (and allows users to define their own), along with various correction methods, enabling easy comparison between approaches. It also integrates with the molecular visualization software `PyMOL <https://pymol.org>`_, allowing users to map results onto 3D structures, more specifically, extremal co-evolving residues that form the basis of protein sectors :cite:p:`jullien_cocoatree_2025`.

The software is organized into different modules, detailed below:

.. image:: schemas/draft_fig1_v5.png
	:width: 100%

List of abbreviations
---------------------
- MI: Mutual Information
- MSA: Multiple Sequence Alignment
- SCA: Statistical Coupling Analysis, the method underlying the identification of protein sectors :cite:p:`halabi_protein_2009,rivoire_evolution-based_2016`
- XCoR: eXtremal Co-evolving Residues

Example gallery
---------------

Various examples of COCOA-Tree's use can be found in the `Gallery <auto_examples/index.html>`_.

If you wish to perform a minimal SCA analysis, go to
`Minimal SCA example <../examples/a_quick_start/plot_minimal_sca.py>`_.

For a detailed SCA analysis, go to
`Full SCA analysis <../examples/b_advanced/plot_full_SCA_analysis.py>`_.

Citing COCOA-Tree
-----------------

If you use COCOA-Tree in a scientific publication, we would appreciate citations
to the following paper:

   BEST PAPER EVER

Bibtex entry::

   @article {}

License information
-------------------

Conditions on the use and redistribution of this package.

.. literalinclude:: ../LICENSE

.. bibliography::

