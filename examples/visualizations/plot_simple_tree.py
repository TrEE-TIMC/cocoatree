"""
=======================================
Simple tree visualization with metadata
=======================================

A simple tree visualization with cocoatree.

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

tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot,
    meta_data=('Protein_type', 'Class', 'Family'),
    fig_title='A tree visualization with metadata',
    )


tree_ete3.render("simple_tree.png", tree_style=tree_style)
