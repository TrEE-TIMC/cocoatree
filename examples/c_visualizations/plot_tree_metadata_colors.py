"""
============================
Specifying metadata's colors
============================

The user can provide a dictionnary of colors for metadata.

"""
# %%
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.visualization import update_tree_ete3_and_return_style

# %%
# Load metadata
# -------------
# Here we use metadata associated with the serine protease dataset but you
# can use your own metadata by importing a csv file as a pandas dataframe.
# Note that your dataframe must have a `Seq_ID` column with identifiers that
# are identical to the ones in your MSA and in your newick tree.
serine_dataset = load_S1A_serine_proteases('halabi')
df_annot = serine_dataset["metadata"]
print(df_annot)

# %%
# Define a custom colormap as a dictionary. Each element in the dictionary
# corresponds to a category in the metadata dataframe that you wish to
# display. In this example, we use the columns `Protein_type`, `Subphylum`,
# and `Class`.
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
# Load sequence data
# ------------------
# Load the sequences you wish to visualize with `cocoatree.io.load_msa()` as
# a fasta file. The sequence names must correspond to `Seq_ID` and to the leaf
# names in the tree file.
#
# See the Perform full SCA analysis example on how to export a fasta of your
# sectors.
sector_file = 'data/halabi_sector_1_SCA.fasta'
data = load_MSA(sector_file, 'fasta')
sector_id = data["sequence_ids"]
print(sector_id[:5])
sector_seq = data["alignment"]
print(sector_seq[:5])

# %%
# Load tree file
# --------------
tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)
print(tree_ete3)

# %%
# Plot figure
# -----------
# Create the tree style that will be applied to your tree
tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot,
    sector_id=sector_id,
    sector_seq=sector_seq,
    meta_data=('Protein_type', 'Subphylum', 'Class'),
    fig_title='Personnalized colormap example',
    metadata_colors=halabi_cmap,
    linewidth=3,
    show_leaf_name=False,
    t_sector_seq=True,
    t_sector_heatmap=False,
    colormap='GnBu'
    )

# %%
# Save the figure:
tree_ete3.render("sector_phylogeny_colors.png", tree_style=tree_style)
