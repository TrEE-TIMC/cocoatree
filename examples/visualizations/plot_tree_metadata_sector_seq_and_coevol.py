"""
==========================================================
Plot sector together with (phylogenetic) tree and metadata
==========================================================

A small example that shows how to plot the sector composition
within a tree and to add metadata

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
# To use your own metadata file, import a csv file as a pandas dataframe.
#
# The dataframe must have a 'Seq_ID' column, which corresponds to the sequence
# identifiers used in the fasta and in the phylogenetic tree. The other columns
# can contain qualitative data that will be displayed as categories alongside
# the phylogenetic tree. In this example, we will use the last 3 columns
# 'Protein_type', 'Class', and 'Family'.

# %%
# Import tree file
# ----------------
# The file must be in Newick format and can include bootstrap values (only
# one).
tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)
print(tree_ete3)

# %%
# Import sector sequences
# -----------------------
sector_file = 'data/halabi_sector_1_SCA.fasta'
data = load_MSA(sector_file, 'fasta')

sector_id = data["sequence_ids"]
sector_seq = data["alignment"]

# %%
# Plot figure
# -----------
tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id, sector_seq,
    meta_data=('Protein_type', 'Class', 'Family'),
    fig_title='Visualization example',
    t_sector_seq=True,
    t_sector_heatmap=True,
    colormap='GnBu'
    )

# %%
# Save the image file
tree_ete3.render("sector_phylogeny.png", tree_style=tree_style)

# %%
# Use ete3's `tree.show()` method for displaying the figure in ete3's
# interactive GUI.
