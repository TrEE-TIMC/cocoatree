"""Module to visualize phylogenetic trees along with sectors"""

# User provided file:
# - phylogenetic tree in newick format
# - multiple sequence alignment used to generate the tree in fasta format
# - annotation table in csv format

# Import necessary packages
from ete3 import ProfileFace, TreeStyle, NodeStyle, TextFace, \
    add_face_to_node, SeqMotifFace, RectFace
import pandas as pd  # type: ignore
from pandas.api.types import is_numeric_dtype  # type: ignore
import numpy as np
from Bio import AlignIO
import math
from PyQt5 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

from .io import load_tree
from .msa import filter_seq_id
from .statistics.pairwise import compute_seq_identity


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
    t, id_lst = load_tree(tree)
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


def _get_color_gradient(self):
    """
    Function which allows to use matplotlib colormaps in ete3 heatmap
    Adapted from:
    https://github.com/lthiberiol/virfac/blob/master/get_color_gradient.py
    """
    cNorm = colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm,
                                   cmap=plt.get_cmap(self.colorscheme))
    color_scale = []
    for scale in np.linspace(0, 1, 201):
        [r, g, b, a] = scalarMap.to_rgba(scale, bytes=True)
        color_scale.append(QtGui.QColor(r, g, b, a))

    return color_scale


def _get_sector_seq(sector_fasta):
    """
    Get the amino acid sequences of the sector (fasta format)

    Arguments
    ---------
    sector_fasta : path to the fasta of the sector sequences

    Returns
    -------
    sector : list of sequences

    sector_length : number of residues in the sector
    """
    sector = {}
    fasta = AlignIO.read(sector_fasta, "fasta")
    sector_length = len(fasta[0].seq)

    for i in range(0, len(fasta)):
        sector[fasta[i].id] = str(fasta[i].seq)

    return sector, sector_length


# DON'T REVIEW YET, I STILL HAVE THINGS TO CHECK
def plot_coev_along_phylogeny(tree_file, annot_file, sector_fasta, attributes,
                              fig_title, rectface=True, seqmotif=True,
                              heatmap=True, colormap='inferno'):
    """
    Wrapping function that draws the phylogenetic tree along with specified
    sector characteristics.

    Arguments
    ---------
    tree_file : path to the tree (Newick format)

    annot_file : path to the annotation file (csv format)

    sector_fasta : path to the fasta containing the sequences to display

    attributes : list of annotations to display, should be the same as in the
        annotation file (should be a list, even if there is only one attribute)

    fig_title : figure title (str)

    rectface : boolean,
        whether to add a RectFace of each attribute

    seqmotif : boolean,
        whether to add a SeqMotifFace of the sequences

    heatmap : boolean,
        whether to add a heatmap of the identity matrix between sector
        sequences
    """

    t, id_lst = load_tree(tree_file)
    nb_seq = len(id_lst)

    ts = TreeStyle()
    ts.layout_fn = []

    # Add bootstrap support NodeStyle
    boot_style = NodeStyle()
    boot_style["fgcolor"] = "darkred"
    boot_style["size"] = 10
    empty_style = NodeStyle()
    empty_style["size"] = 0
    for node in t.traverse():
        if node.support >= 95:
            node.set_style(boot_style)
        else:
            node.set_style(empty_style)

    col_rectface = 0
    col_legend_rectface = 0
    if rectface is True:
        for att in attributes:
            attribute_colors, col_dict = _annot_to_color(att, tree_file,
                                                         annot_file)

            def layout_RectFace(node):
                if node.is_leaf():
                    name = node.name
                    square = RectFace(50, 20, fgcolor=attribute_colors[name],
                                      bgcolor=attribute_colors[name])
                    square.margin_left = 10
                    # TO DO: For now, works only when there is only one
                    # attribute to represent
                    add_face_to_node(square, node, column=col_rectface,
                                     position='aligned')
            ts.layout_fn.append(layout_RectFace)
            col_rectface += 1
            # Add legend
            ts.legend.add_face(TextFace(att, fsize=10, bold=True),
                               column=col_legend_rectface)
            # otherwise text is not in front of RectFace
            ts.legend.add_face(TextFace(""), column=col_legend_rectface + 1)
            for gene in col_dict.keys():
                legend_face = RectFace(50, 20, fgcolor=col_dict[gene],
                                       bgcolor=col_dict[gene])
                legend_face.margin_right = 5
                ts.legend.add_face(legend_face, column=col_legend_rectface)
                ts.legend.add_face(TextFace(gene, fsize=10),
                                   column=col_legend_rectface + 1)
            col_legend_rectface += 2

    col_seqmotif = 0
    if seqmotif is True:
        sector_seq, sector_length = _get_sector_seq(sector_fasta)
        if rectface is True:
            col_seqmotif = col_rectface

        def layout_SeqMotifFace(node):
            if node.is_leaf():
                name = node.name
                if name in sector_seq:
                    seq = sector_seq[name]
                else:
                    seq = '-' * sector_length
                seqFace = SeqMotifFace(seq,
                                       motifs=[[0, sector_length, "seq", 20,
                                                20, None, None, None]],
                                       scale_factor=1)
                seqFace.margin_left = 10
                seqFace.margin_right = 10
                add_face_to_node(seqFace, node, column=col_seqmotif,
                                 position='aligned')
        ts.layout_fn.append(layout_SeqMotifFace)
        col_seqmotif += 1

    col_heatmap = 0
    if heatmap is True:
        # allow to chose among Matplotlib's colormaps
        ProfileFace.get_color_gradient = _get_color_gradient
        # Check that sequences in the similarity matrix are ordered as in the
        # tree leaves
        msa = AlignIO.read(sector_fasta, 'fasta')
        reorder_msa = filter_seq_id(msa, id_lst)
        id_mat = compute_seq_identity(reorder_msa[2])
        # Define the column in which the heatmap will be
        if (rectface is True) & (seqmotif is False):
            col_heatmap = col_rectface + 1
        elif (rectface is True) & (seqmotif is True):
            col_heatmap = col_seqmotif + 1
        elif (rectface is False) & (seqmotif is False):
            col_heatmap = 0
        elif (rectface is False) & (seqmotif is True):
            col_heatmap = col_seqmotif + 1

        count = 0
        # Add heatmap profile to each leaf
        for lf in t.iter_leaves():
            lf.add_features(profile=id_mat[count])
            count += 1
            lf.add_features(deviation=[0 for x in range(id_mat.shape[0])])
            lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5,
                                    width=(nb_seq*20), height=20,
                                    style='heatmap',
                                    colorscheme=colormap),
                        column=col_heatmap, position="aligned")

    # Add title
    ts.title.add_face(TextFace(fig_title, fsize=20), column=0)

    return t.show(tree_style=ts)
