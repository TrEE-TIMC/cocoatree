"""
============================
Specifying metadata's colors
============================

The user can provide a dictionnary of colors for metadata.

"""

import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.visualization import update_tree_ete3_and_return_style


serine_dataset = load_S1A_serine_proteases('halabi')
df_annot = serine_dataset["metadata"]

tree_file = 'data/halabi_82_seqs.txt'
tree_ete3 = load_tree_ete3(tree_file)

sector_file = 'data/halabi_sector_1_SCA.fasta'
data = load_MSA(sector_file, 'fasta')

sector_id = data["sequence_ids"]
sector_seq = data["alignment"]


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


tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id, sector_seq,
    meta_data=('Protein_type', 'Class', 'Family'),
    fig_title='Personnalized colormap example',
    metadata_colors=halabi_cmap,
    linewidth=3,
    show_leaf_name=False,
    t_sector_seq=True,
    t_sector_heatmap=False,
    colormap='GnBu'
    )


tree_ete3.render("sector_phylogeny_colors.png", tree_style=tree_style)
