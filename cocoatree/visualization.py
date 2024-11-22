"""Module to visualize phylogenetic trees along with sectors"""

# User provided file:
# - phylogenetic tree in newick format
# - multiple sequence alignment used to generate the tree in fasta format
# - annotation table in csv format

# Import necessary packages
from ete3 import Tree, ProfileFace, TreeStyle, NodeStyle, TextFace, \
    add_face_to_node, SeqMotifFace, RectFace, NCBITaxa
import json
import pandas as pd  # type: ignore
from pandas.api.types import is_numeric_dtype  # type: ignore
import numpy as np
from Bio import AlignIO, Entrez
import math
from PyQt5 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

from .msa import filter_seq_id
from .statistics.pairwise import compute_seq_identity


def import_tree(tree):
    """
    Import tree (Newick format) and get list of sequence IDs

    Arguments
    ---------
    tree : path to the newick file

    Returns
    -------
    t : tree object

    id_lst : list of the tree leaves (sequence IDs)

    Example
    -------
    t, id_lst = import_tree(tree)
    """
    t = Tree(tree, format=0)
    id_lst = t.get_leaf_names()
    return t, id_lst


def _annot_to_color(attribute, tree, annot_file):
    """
    Reads in the attributes specified by the user in the annotation csv file
    and attributes a color palette for each.

    Arguments
    ---------
    tree : path to the tree (Newick format)

    attributes : list of column names to grab

    annot_file : path to the annotation file

    Returns
    -------
    att_dict:

    color_dict:
    """
    t, id_lst = import_tree(tree)
    df_annot = pd.read_csv(annot_file)
    df_annot = df_annot.fillna('unknown')
    if is_numeric_dtype(df_annot['Seq_ID']):
        df_annot['Seq_ID'] = df_annot['Seq_ID'].astype('str')
    df_annot = df_annot[df_annot['Seq_ID'].isin(id_lst)]

    att_dict = {}
    df_annot = df_annot[['Seq_ID', attribute]]
    color_dict = _get_color_palette(df_annot[attribute].unique())
    df_annot[str(attribute + '_color')] = df_annot.apply(
        lambda row: color_dict[row[attribute]], axis=1)
    for i in range(0, len(df_annot['Seq_ID'])):
        row = df_annot.iloc[i].tolist()
        att_dict[row[0]] = row[2]

    return att_dict, color_dict


# TO DO: add check if value == 'unknown': color = 'white'
def _get_color_palette(values):
    # color palettes (modified from colorbrewer set1, expanded to 50)
    colors_50 = ["#E41A1C", "#C72A35", "#AB3A4E", "#8F4A68", "#735B81",
                 "#566B9B", "#3A7BB4", "#3A85A8", "#3D8D96", "#419584",
                 "#449D72", "#48A460", "#4CAD4E", "#56A354", "#629363",
                 "#6E8371", "#7A7380", "#87638F", "#93539D", "#A25392",
                 "#B35A77", "#C4625D", "#D46A42", "#E57227", "#F67A0D",
                 "#FF8904", "#FF9E0C", "#FFB314", "#FFC81D", "#FFDD25",
                 "#FFF12D", "#F9F432", "#EBD930", "#DCBD2E", "#CDA12C",
                 "#BF862B", "#B06A29", "#A9572E", "#B65E46", "#C3655F",
                 "#D06C78", "#DE7390", "#EB7AA9", "#F581BE", "#E585B8",
                 "#D689B1", "#C78DAB", "#B791A5", "#A8959F", "#999999"]
    simple_palette = ["HotPink", "LimeGreen", "DodgerBlue", "Turquoise",
                      "Indigo", "MediumTurquoise", "Sienna", "LightCoral",
                      "LightSkyBlue", "Indigo", "Tan", "Coral",
                      "OliveDrab", "Teal"]

    nvals = len(values)
    if nvals > 14:
        if nvals > 50:
            color_list = colors_50 + int(nvals/50) * colors_50
        else:
            # every nth colour
            color_list = colors_50[::int(math.floor(50/nvals))]
    else:
        color_list = simple_palette
    color_dict = {}  # key = value, value = colour id
    for i in range(0, nvals):
        if values[i] == 'unknown':
            color_dict[values[i]] = 'white'
        else:
            color_dict[values[i]] = color_list[i]

    return color_dict
