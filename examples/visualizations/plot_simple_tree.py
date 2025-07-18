"""
=======================================
Simple tree visualization with metadata
=======================================

A simple tree visualization with cocoatree.

"""

import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style


annot_file = 'data/random_seq_tree_annotations.csv'
df_annot = pd.read_csv(annot_file)

tree_file = 'data/random_seq_tree_kpicsg_ete3.treefile'
tree_ete3 = load_tree_ete3(tree_file)

sector_file = 'data/SCA_sector_1_core_tagged_seq.fasta'
data = load_MSA(sector_file, 'fasta')

sector_id = data["sequence_ids"]
sector_seq = data["alignment"]

tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot,
    meta_data=("superkingdom", 'class', 'HMM_annotation'),
    fig_title='A tree visualization with metadata',
    )


tree_ete3.render("simple_tree.png", tree_style=tree_style)
