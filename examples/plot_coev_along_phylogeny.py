"""
==============================
Plot final figure
==============================

A small example that shows how to use the function
plot_coev_along_phylogeny of the visualization module.

"""

import pandas as pd
from cocoatree.io import load_MSA, load_tree
from cocoatree.visualization import plot_coev_along_phylogeny

annot_file = 'data/random_seq_tree_annotations.csv'
df_annot = pd.read_csv(annot_file)
tree_file = 'data/random_seq_tree_kpicsg_ete3.treefile'
tree, leaves_id = load_tree(tree_file)
sector_file = 'data/SCA_sector_1_core_tagged_seq.fasta'
sector_id, sector_seq = load_MSA(sector_file, 'fasta')

tree, tree_style = plot_coev_along_phylogeny(
    tree, df_annot, sector_id, sector_seq,
    attributes=('HMM_annotation',),
    fig_title='Test HMM annot',
    rectface=True,
    seqmotif=True,
    heatmap=True,
    colormap='inferno')
