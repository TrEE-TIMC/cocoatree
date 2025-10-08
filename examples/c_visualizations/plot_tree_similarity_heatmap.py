"""
=================================================================
Plot a similarity heatmap of a sector along the phylogenetic tree
=================================================================

Here we present how you can plot a heatmap of sequence similarity
ordered following a phylogenetic tree.

"""

# %%
# Import necessary packages
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.visualization import update_tree_ete3_and_return_style

# %%
# Import metadata
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
# 'Protein_type', 'Subphylum', and 'Class'.
#
# We will also use a personnalized colormap defined as follows:
halabi_cmap = {
    'vertebrate': '#798e87',
    'invertebrate': '#c27d38',
    'fungi': '#ccc591',
    'bacteria': '#29211f',
    'Chymotrypsin': '#00a08a',
    'Trypsin': '#ff0000',
    'Tryptase': '#f2ad00',
    'Kallikrein': '#f98400',
    'Granzyme': '#5bbcd6',
    'Mammalia': '#C969A1',
    'Actinopterygii': '#CE4441',
    'Amphibia': '#EE8577',
    'Malacostraca': '#EB7926',
    'other': 'lightgrey',
    'Insecta': '#FFBB44',
    'Actinobacteria': '#859B6C',
    'Arachnida': '#62929A',
    'Oligochaeta': '#004F63'
}

# %%
# Import tree file
# ----------------
# The file must be in Newick format and can include confidence scores such as
# bootstrap or jack-knife (only one type of confidence score at a time).
tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)
print(tree_ete3)

# %%
# Import sector sequences
# -----------------------
# Load the sequences you wish to visualize with `cocoatree.io.load_msa()` as
# a fasta file. The sequence names must correspond to `Seq_ID` and to the leaf
# names in the tree file.
sector_file = 'data/halabi_sector_1_SCA.fasta'
data = load_MSA(sector_file, 'fasta')
sector_id = data["sequence_ids"]
sector_seq = data["alignment"]

# %%
# Plot tree with sequence similarity heatmap
# ------------------------------------------
# Generate the tree style that will be applied to your tree.
#
# Here, we will show various elements:
#   - the tree without its leaf names (`show_leaf_name=False`)
#   - metadata as colored columns, in order: *Protein_type*, *Subphylum*,
# and *Class*.
#   - sector sequences colored by amino acid physico-chemical properties
# (`t_sector_seq=True`)
#   - a heatmap of pairwise sequence similarity computed on the sector
# sequences (`t_sector_heatmap=True`, `matrix_type='similarity'`) and
# using the `GnBu` colormap (`colormap='GnBu'`).
tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id, sector_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    show_leaf_name=False,
    fig_title='Heatmap of sequence similarity',
    linewidth=3,
    linecolor="#000000",
    bootstrap_style={},
    tree_scale=200,
    metadata_colors=halabi_cmap,
    t_sector_seq=True,
    t_sector_heatmap=True,
    matrix_type='similarity',
    colormap='GnBu'
    )

# %%
# Save the image file
tree_ete3.render("sector_phylogeny.png", tree_style=tree_style)

# %%
# You can use ete3's `tree.show()` method for displaying the figure in
# ete3's interactive GUI.

# %%
# Plot tree with normalized similarity heatmap
# --------------------------------------------
# You first need to reload the tree in order to clean the tree style
# defined above
tree_ete3 = load_tree_ete3(tree_file)
tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id, sector_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    show_leaf_name=False,
    fig_title='Heatmap of normalized sequence similarity',
    linewidth=3,
    linecolor="#000000",
    bootstrap_style={},
    tree_scale=200,
    metadata_colors=halabi_cmap,
    t_sector_seq=True,
    t_sector_heatmap=True,
    matrix_type='norm_similarity',
    colormap='GnBu'
    )
tree_ete3.render("sector_phylogeny.png", tree_style=tree_style)

# %%
# Plot heatmap of a single sector position
# ----------------------------------------
# You can compute sequence similarity on a subset of positions, here
# we look at position 189 (based on 3TGI pdb sequence), which is the
# third residue in the sector sequence.

tree_ete3 = load_tree_ete3(tree_file)
subsector_seq = [res[2] for res in sector_seq]
tree_style, _ = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id,
    sector_seq=subsector_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    show_leaf_name=False,
    fig_title='Normalized sequence similarity of position 189',
    linewidth=3,
    linecolor="#000000",
    bootstrap_style={},
    tree_scale=200,
    metadata_colors=halabi_cmap,
    t_sector_seq=True,
    t_sector_heatmap=True,
    matrix_type='norm_similarity',
    colormap='GnBu'
    )
tree_ete3.render("sector_phylogeny.png", tree_style=tree_style)
