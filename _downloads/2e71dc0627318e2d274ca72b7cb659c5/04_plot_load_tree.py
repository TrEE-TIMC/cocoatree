"""
================
Load a tree file
================

This is a small example expliciting how to load a tree file in :mod:`cocoatree`
using :func:`cocoatree.io.load_tree_ete3` function.
"""

# %%
from cocoatree.io import load_tree_ete3

# %%
# Let's start by creating a simplified tree:
simple_tree = "(A:1,(B:1,(C:1,D:1):0.5):0.5);"
t = load_tree_ete3(simple_tree)
print(t)

# %%
# Here, we have a tree with distances and leaf names. The distances are
# specified after the ``:`` character and the leaf names are A, B, C, and D.

# %%
# :func:`cocoatree.io.load_tree_ete3` also accepts trees with more information
# such as branch support values encoded as follow:
bootstrap_tree = "(A:1,(B:1,(C:1,D:1)95:0.5)98:0.5);"
t_boot = load_tree_ete3(bootstrap_tree)

# %%
# The branch support values are placed following their corresponding closing
# bracket.

# %%
# .. warning::
#       Only one type of branch support value is accepted (i.e. a file with
#       bootstrap and jack-knife branch support values will not be accepted).
#
#       Branch support values should not be between square brackets as those
#       are not accepted.

# %%
# .. seealso::
#       :ref:`sphx_glr_auto_examples_c_visualizations_plot_simple_tree.py` for
#       a simple tree visualization or :ref:`sphx_glr_auto_examples_\
#       c_visualizations_plot_tree_metadata_xcor_seq_and_coevol.py` for a more
#       complete visualization.
