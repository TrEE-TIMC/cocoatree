# COCOA-Tree: COllaborative COevolution Analysis Toolbox

**COCOA-Tree** is a Python library to perform coevolution analyses of proteins and integrate phylogenetic information
to better understand the coevolution signals.

Website: [tree-bioinfo-intra.timc.fr/projects/cocoa/](http://tree-bioinfo-intra.timc.fr/projects/cocoa/index.html)

## Library organization

[cocoatree_orga](../cocoatree_orga.pdf)

## Installation

### Dependencies

COCOA-Tree requires:
- Python (>= 3.9)
- NumPy
- scikit-learn
- biopython
- ete3

For running the examples, Matplotlib is required.

### User installation

To install for development purposes, use::

    python setup.py develop


All changes in the python module will be repercuted on the installed version
(in practice the installed version is a symlink towards the package in
development).

