"""Module to visualize phylogenetic trees along with sectors"""

# User provided file:
# - phylogenetic tree in newick format
# - multiple sequence alignment used to generate the tree in fasta format
# - annotation table in csv format

# Import necessary packages
from ete3 import ProfileFace, TreeStyle, NodeStyle, TextFace, \
    add_face_to_node, SeqMotifFace, RectFace
from pandas.api.types import is_numeric_dtype  # type: ignore
import pandas as pd
import numpy as np
from PyQt5 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

from .msa import compute_seq_identity, compute_seq_similarity, \
    compute_normalized_seq_similarity


def _annot_to_color(attribute, tree, df_annot, cmap='jet'):
    """
    Reads in the attributes specified by the user in the annotation csv file
    and attributes a color palette for each.

    Parameters
    ----------
    tree : ete3's tree object,
            as imported by io.load_tree_ete3()

    attributes : list of column names to grab

    df_annot : pandas dataframe of the annotation file

    Returns
    -------
    att_dict : dictionnary in which keys are the sequence IDs and the values
            are the colors associated with it

    color_dict : dictionnary in which keys are the attribute's categories and
            the values are the colors associated to each category
    """
    id_lst = tree.get_leaf_names()
    df_annot = df_annot.fillna('unknown')
    if is_numeric_dtype(df_annot['Seq_ID']):
        df_annot['Seq_ID'] = df_annot['Seq_ID'].astype('str')
    df_annot = df_annot[df_annot['Seq_ID'].isin(id_lst)]

    att_dict = {}
    df_annot = df_annot[['Seq_ID', attribute]]

    if isinstance(cmap, str):
        color_dict = _get_color_palette(
            list(df_annot[attribute].unique()), cmap)
    else:
        color_dict = {n: cmap[n] for n in list(df_annot[attribute].unique())}

    df_annot[str(attribute + '_color')] = df_annot.apply(
        lambda row: color_dict[row[attribute]], axis=1)
    for i in range(0, len(df_annot['Seq_ID'])):
        row = df_annot.iloc[i].tolist()
        att_dict[row[0]] = row[2]

    return att_dict, color_dict


def _generate_colors_from_colormaps(n_colors, cmap="jet", as_hex=True):
    """
    Generate a list of n colors from colormap
    """

    colormap = plt.get_cmap(str(cmap))
    indx = np.linspace(0, 1, n_colors)
    indexed_colors = [colormap(i) for i in indx]
    if as_hex:
        indexed_colors = [colors.to_hex(i) for i in indexed_colors]
    return indexed_colors


# TO DO: allow user to choose which color to use for 'unknown'
# (currently: white by default)
def _get_color_palette(values, cmap):

    nvals = len(values)
    colors = _generate_colors_from_colormaps(nvals, cmap=cmap, as_hex=True)

    color_dict = {}  # key = value, value = colour id
    for i in range(0, nvals):
        if values[i] == 'unknown':
            color_dict[values[i]] = '#FFFFFF'
        else:
            color_dict[values[i]] = colors[i]

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


def update_tree_ete3_and_return_style(
        tree_ete3, df_annot,
        sector_id=None,
        sector_seq=None,
        meta_data=None,
        show_leaf_name=True,
        fig_title='',
        linewidth=1,
        linecolor="#000000",
        bootstrap_style={},
        tree_scale=200,
        metadata_colors=None,
        t_sector_seq=False,
        t_sector_heatmap=False,
        matrix_type='identity',
        colormap='inferno'
        ):
    """
    Update ete3 tree with sector info and attributes
    and return tree_style for further visualization.

    Parameters
    ----------
    tree_ete3 : ete3's tree object,
            as imported by io.load_tree_ete3()

    annot_file : pandas dataframe of the annotation file

    sector_id : list of sector identifiers, as imported by io.load_msa()
            the ids must match with the tree's leaves id

    sector_seq : corresponding list of sector sequences to display,
            as imported by io.load_msa()

    meta_data : tuple of annotations to display
                 (from annotation file's header)

    show_leaf_name : boolean, optional, default: True
        whether to show leaf names.

    linewidth : int, optional, default: 1
        width of the lines in the tree

    linecolor : str, optional, default: "#000000"
        color of the lines

    bootstrap_style : dict, optional,
        `fgcolor`: color of the bootstrap node, default: "darkred"
        `size`: size of the bootstrap node, default: 10
        `support`: int between 0 and 100, minimum support level for display

    tree_scale : int, optional, default: 200
        sets the scale of the tree in ETE3: the higher, the larger the tree
        will be (in width)

    metadata_colors : dict, str, or None, optional, default: None
        colors for the metadata:
            - None: generates automatically the colors
            - str: uses a Matplotlib colormap to generate the colors
            - dict: specifies colors for each matadata entry
                {key: color}

    fig_title : figure title (str)

    t_sector_seq : boolean,
        whether to show the sequences of the sector

    t_sector_heatmap : boolean,
        whether to add a heatmap of the identity or similarity matrix between
        sector sequences

    matrix_type : str, default='identity'
        whether to compute pairwise sequence identity ('identity'), similarity
        ('similarity'), or normalized similarity ('norm_similarity').

    colormap : str, default='inferno'
        the matplotlib colormap to use for the heatmap

    Returns
    -------
    tree_style : TreeStyle class from ete3

    column_end : int, the number of columns after the tree. If you want to
        plot anything else alongside the tree, the column number should be
        equal to this value.

    """

    tree_style = TreeStyle()
    tree_style.scale = tree_scale
    tree_style.layout_fn = []
    # tree_style.branch_vertical_margin = 20
    tree_style.show_leaf_name = show_leaf_name

    # Add bootstrap support NodeStyle
    boot_style = NodeStyle()
    boot_style["fgcolor"] = \
        bootstrap_style["fgcolor"] if "fgcolor" in bootstrap_style \
        else "darkred"
    boot_style["size"] = \
        bootstrap_style["size"] if "size" in bootstrap_style else 10
    support = \
        bootstrap_style["support"] if "support" in bootstrap_style else 95

    boot_style["hz_line_width"] = linewidth
    boot_style["vt_line_width"] = linewidth
    boot_style["vt_line_color"] = linecolor
    boot_style["hz_line_color"] = linecolor

    empty_style = NodeStyle()
    empty_style["size"] = 0
    empty_style["vt_line_width"] = linewidth
    empty_style["hz_line_width"] = linewidth
    empty_style["vt_line_color"] = linecolor
    empty_style["hz_line_color"] = linecolor

    for node in tree_ete3.traverse():
        if node.support >= support:
            node.set_style(boot_style)
        else:
            node.set_style(empty_style)

    column_layout = 0
    col_legend_rectface = 0

    if metadata_colors is None:
        metadata_colors = "jet"

    # If no metadata, do nothing
    if meta_data:

        def layout_attribute(node, column=column_layout):
            if node.is_leaf():
                name = node.name
                rect_faces = [None for i in range(len(meta_data))]
                for i, col in enumerate(meta_data):
                    colors, _ = _annot_to_color(col,
                                                tree_ete3,
                                                df_annot,
                                                cmap=metadata_colors)

                    rect_faces[i] = RectFace(50, 20,
                                             fgcolor=colors[name],
                                             bgcolor=colors[name])
                    rect_faces[i].margin_left = 5
                    rect_faces[i].margin_right = 0
                    if i == len(meta_data) - 1:
                        rect_faces[i].margin_right = 30
                    add_face_to_node(rect_faces[i], node, column=column,
                                     position='aligned')
                    column += 1

        tree_style.layout_fn.append(layout_attribute)

        # Add legend
        legend_face = [None for i in range(len(meta_data))]
        for i, col in enumerate(meta_data):
            _, col_dict = _annot_to_color(col, tree_ete3,
                                          df_annot, cmap=metadata_colors)
            tree_style.legend.add_face(TextFace(col,
                                                fsize=10,
                                                bold=True),
                                       column=col_legend_rectface)
            # otherwise text is not in front of RectFace
            tree_style.legend.add_face(TextFace(""),
                                       column=col_legend_rectface + 1)

            legend_face[i] = {key: None for key in col_dict.keys()}
            for key in col_dict.keys():
                legend_face[i][key] = RectFace(50, 20, fgcolor=col_dict[key],
                                               bgcolor=col_dict[key])
                legend_face[i][key].margin_right = 5
                legend_face[i][key].margin_left = 10
                tree_style.legend.add_face(legend_face[i][key],
                                           column=col_legend_rectface)
                tree_style.legend.add_face(TextFace(key, fsize=10),
                                           column=col_legend_rectface + 1)
            col_legend_rectface += 2
    column_layout += len(meta_data) if meta_data else 0

    if t_sector_seq:
        tree_style, column_layout = add_sector_sequences_to_tree(
            tree_style, tree_ete3, sector_id,
            sector_seq, column_start=column_layout)

    if t_sector_heatmap:
        tree_style, column_layout = add_heatmap_to_tree(
            tree_style, tree_ete3, sector_id, sector_seq,
            matrix_type=matrix_type,
            column_start=column_layout, colormap=colormap)

    # Add title
    tree_style.title.add_face(TextFace(fig_title, fsize=20), column=0)

    return tree_style, column_layout


def add_sector_sequences_to_tree(tree_style, tree_ete3, sector_id, sector_seq,
                                 column_start=0):
    """
    Add sector sequence to ETE3's tree style

    Parameters
    ----------
    tree_style : ETE3's tree_style object

    tree_ete3 : ete3's tree object,
            as imported by io.load_tree_ete3()

    sector_id : list of sector identifiers, as imported by io.load_msa()
            the ids must match with the tree's leaves id

    sector_seq : corresponding list of sector sequences to display,
            as imported by io.load_msa()

    column_start : int, optional, default : 0
        the column on which to start plotting

    Returns
    -------
    tree_style : TreeStyle class from ete3

    column_end : int, the number of columns after the tree. If you want to
        plot anything else alongside the tree, the column number should be
        equal to this value.

    """
    sector_dict = {
        sector_id[i]: str(sector_seq[i]) for i in range(len(sector_id))}

    def layout_SeqMotifFace(node, column=column_start):
        if node.is_leaf():
            if node.name in sector_dict:
                seq = sector_dict[node.name]
            else:
                seq = '-' * len(sector_seq[0])
            seqFace = SeqMotifFace(seq,
                                   motifs=[[0, len(sector_seq[0]), "seq",
                                            20, 20, None, None, None]],
                                   scale_factor=1)
            seqFace.margin_right = 30
            add_face_to_node(seqFace, node, column=column,
                             position='aligned')
    tree_style.layout_fn.append(layout_SeqMotifFace)
    column_start += 1
    return tree_style, column_start


def add_heatmap_to_tree(tree_style, tree_ete3, sector_id, sector_seq,
                        matrix_type="identity",
                        column_start=0, width=20, colormap="inferno"):
    """
    Add heatmap to ETE3's tree style

    Parameters
    ----------
    tree_style : ETE3's tree_style object

    tree_ete3 : ete3's tree object,
            as imported by io.load_tree_ete3()

    sector_id : list of sector identifiers, as imported by io.load_msa()
            the ids must match with the tree's leaves id

    sector_seq : corresponding list of sector sequences to display,
            as imported by io.load_msa()

    matrix_type : str, default='identity'
            whether to compute pairwise matrix identity ('identity'),
            similarity ('similarity'), or normalized similarity
            ('norm_similarity')

    column_start : int, optional, default : 0
        the column on which to start plotting

    width : int, optional, default : 20
        the width of each square of the heatmap. If width == 20, the heatmap
        will be squared.

    colormap : str, optional, default: "inferno"
        any Matplotlib's colormap

    Returns
    -------
    tree_style : TreeStyle class from ete3

    column_end : int, the number of columns after the tree. If you want to
        plot anything else alongside the tree, the column number should be
        equal to this value.
    """

    leaves_id = tree_ete3.get_leaf_names()
    nb_leaves = len(leaves_id)

    # allow to chose among Matplotlib's colormaps
    ProfileFace.get_color_gradient = _get_color_gradient

    # Check that sequences in the similarity matrix are ordered as in the
    # tree leaves and keep only sequences that are present in the tree
    sequences = pd.DataFrame(index=sector_id, data={"seq": sector_seq})
    reordered_sequences = sequences.loc[leaves_id, "seq"].values

    if matrix_type == 'identity':
        matrix = compute_seq_identity(reordered_sequences)
        # FIX to zero values appearing black in the heatmap whatever the cmap
        matrix[matrix == 0] = 0.00000001
    elif matrix_type == 'similarity':
        matrix = compute_seq_similarity(reordered_sequences)
    elif matrix_type == 'norm_similarity':
        matrix = compute_normalized_seq_similarity(reordered_sequences)

    # Add heatmap profile to each leaf
    for i, lf in enumerate(tree_ete3.iter_leaves()):
        lf.add_features(profile=matrix[i])
        lf.add_features(deviation=[0 for x in range(matrix.shape[0])])
        lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5,
                                width=(nb_leaves*width), height=20,
                                style='heatmap',
                                colorscheme=colormap),
                    column=column_start, position="aligned")
    column_start += nb_leaves*width
    return tree_style, column_start
