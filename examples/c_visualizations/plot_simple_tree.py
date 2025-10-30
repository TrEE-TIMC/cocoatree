"""
=======================================
Simple tree visualization with metadata
=======================================

This example details how to plot a simple tree visualization with cocoatree.

"""

# %%
from cocoatree.io import load_tree_ete3
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.visualization import update_tree_ete3_and_return_style

# %%
# Load metadata
# -------------
# We use the S1A serine protease dataset that is included in cocoatree
# but you can use your own metadata by importing a csv file as a pandas
# dataframe.
#
# For more details on the S1A serine proteases dataset, go to
# :ref:`sphx_glr_auto_examples_d_datasets_plot_s1A_serine_proteases.py`.
serine_dataset = load_S1A_serine_proteases('halabi')
df_annot = serine_dataset["metadata"]
print(df_annot)

# %%
# .. note::
#       The dataframe must have a ``Seq_ID`` column with identifiers that are
#       identical to the ones in the MSA (although no MSA is needed in the
#       present example), and in the newick tree.::

# %%
# Load tree file
# --------------
# Use :func:`cocoatree.io.load_tree_ete3` dedicated function
tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)
print(tree_ete3)

# %%
# .. seealso::
#       Check :ref:`sphx_glr_auto_examples_a_quick_start_plot_load_tree.py`
#       for details on the Newick format that is accepted.

# %%
# Plot figure
# -----------
# Create the tree style that will be applied to your tree with
# :func:`cocoatree.visualization.update_tree_ete3_and_return_style`
tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    fig_title='A tree visualization with metadata',
    )

# %%
# Save the figure:
tree_ete3.render("simple_tree.png", tree_style=tree_style)

# %%
# You can also use ete3's ``tree.show()`` method for displaying the figure
# in ete3's interactive GUI.

# %%
# .. warning::
#       If you wish to modify the style that has been applied to the
#       ``tree_ete3`` object, you need to reload the tree in order to clear
#       out the previous tree style. Otherwise, there will be a
#       superimposition of tree styles.
