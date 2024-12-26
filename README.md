# COCOA-Tree: COllaborative COevolution Analysis Toolbox

**COCOA-Tree** is a Python library to perform coevolution analyses of proteins and integrate phylogenetic information
to better understand the coevolution signals.

It regroups various coevolution metrics and corrections, such as statistical coevolution analysis (SCA) or mutual
information (MI).

Website: [tree-bioinfo-intra.timc.fr/projects/cocoa/](http://tree-bioinfo-intra.timc.fr/projects/cocoa/index.html)

## Library organization

COCOA-Tree is organized in several modules, each allowing to perform different tasks of a coevolution analysis
pipeline:

[cocoatree_orga](../cocoatree_orga.pdf)

## Installation

### Dependencies

COCOA-Tree requires:
- Python (>= 3.9)
- NumPy
- scikit-learn
- biopython
- ete3

Matplotlib is also required for running the examples.

### User installation

To install for development purposes, use::

    python setup.py develop


All changes in the python module will be repercuted on the installed version
(in practice the installed version is a symlink towards the package in
development).

