"""
==========================================================
Plot XCoR together with (phylogenetic) tree and metadata
==========================================================

A small example that shows how to plot the XCoR composition within a tree and
to add metadata.

"""

# %%
# Import necessary packages
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.visualization import update_tree_ete3_and_return_style

# %%
# Import datasets
# ---------------
serine_dataset = load_S1A_serine_proteases('halabi')
df_annot = serine_dataset["metadata"]
print(df_annot)

# %%
# For more details on the S1A serine proteases dataset, go to
# :ref:`sphx_glr_auto_examples_d_datasets_plot_s1A_serine_proteases.py`.

# %%
# To use your own metadata file, import a csv file as a pandas dataframe.
#
# The dataframe must have a 'Seq_ID' column, which corresponds to the sequence
# identifiers used in the fasta and in the phylogenetic tree. The other columns
# can contain qualitative data that will be displayed as categories alongside
# the phylogenetic tree. In this example, we will use the last 3 columns
# 'Protein_type', 'Subphylum', and 'Class'.

# %%
# Import tree file
# ----------------
# The file must be in Newick format and can include bootstrap values (only
# one).
tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)
print(tree_ete3)

# %%
# Import XCoR sequences
# ---------------------
# Load the sequences you wish to visualize with :func:`cocoatree.io.load_MSA`
# as a fasta file. The sequence names must correspond to `Seq_ID` and to the
# leaf names in the tree file.
xcor_file = 'data/halabi_xcor_1_SCA.fasta'
data = load_MSA(xcor_file, 'fasta')
xcor_id = data["sequence_ids"]
xcor_seq = data["alignment"]

# %%
# Plot figure
# -----------
# Generate the tree style that will be applied to your tree.
#
# Here, we will show various elements:
#   - the tree with its leaf names (corresponding to `Seq_ID`)
#   - metadata as colored columns, in order: *Protein_type*, *Subphylum*,
# and *Class*. The default colormap (`jet`), will be used. See
# :ref:`sphx_glr_auto_examples_c_visualizations_plot_tree_metadata_colors.py`
# example on how to modify it.
#   - XCoR sequences colored by amino acid physico-chemical properties
# (`t_xcor_seq=True`)
#   - a heatmap of pairwise sequence identity computed on the XCoR\
# sequences (`t_xcor_heatmap=True`) using the `GnBu` colormap
tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, xcor_id, xcor_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    fig_title='Visualization example',
    t_xcor_seq=True,
    t_xcor_heatmap=True,
    colormap='GnBu'
    )

# %%
# Save the image file
tree_ete3.render("xcor_phylogeny.png", tree_style=tree_style)

# %%
# You can use ete3's `tree.show()` method for displaying the figure in ete3's
# interactive GUI.
