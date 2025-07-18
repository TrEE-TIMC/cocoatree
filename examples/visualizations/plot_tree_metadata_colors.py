"""
============================
Specifying metadata's colors
============================

The user can provide a dictionnary of colors for metadata.

"""

import pandas as pd
from cocoatree.io import load_MSA, load_tree_ete3
from cocoatree.visualization import update_tree_ete3_and_return_style
from cocoatree.visualization import add_heatmap_to_tree


annot_file = 'data/random_seq_tree_annotations.csv'
df_annot = pd.read_csv(annot_file)

tree_file = 'data/random_seq_tree_kpicsg_ete3.treefile'
tree_ete3 = load_tree_ete3(tree_file)

sector_file = 'data/SCA_sector_1_core_tagged_seq.fasta'
data = load_MSA(sector_file, 'fasta')

sector_id = data["sequence_ids"]
sector_seq = data["alignment"]


colors = {
    "Bacteria": "black",
    "Eukaryota": "#cccccc",
    "FMO_Cyano": "#fff5f0",
    "UbiL": "#fdd4c2",
    "UbiI": "#fca082",
    "UbiH": "#fb694a",
    "UbiF": "#e32f27",
    "UbiN": "#b11218",
    "UbiFHILM_Cyano": "#67000d",
    }


tree_style, column_layout = update_tree_ete3_and_return_style(
    tree_ete3, df_annot, sector_id, sector_seq,
    #meta_data=("superkingdom", 'class', 'HMM_annotation'),
    meta_data=("superkingdom", "HMM_annotation"),
    fig_title='Visualization example',
    metadata_colors=colors,
    t_sector_seq=True,
    t_sector_heatmap=False,
    colormap='inferno'
    )


tree_ete3.render("sector_phylogeny_colors.png", tree_style=tree_style)
