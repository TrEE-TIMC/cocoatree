"""Module to visualize phylogenetic trees along with sectors"""

# User provided file:
# - phylogenetic tree in newick format
# - multiple sequence alignment used to generate the tree in fasta format
# - annotation table in csv format

# Import necessary packages
from ete3 import ProfileFace, TreeStyle, NodeStyle, TextFace, \
    add_face_to_node, SeqMotifFace, RectFace
from pandas.api.types import is_numeric_dtype  # type: ignore
import numpy as np
from PyQt5 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

from .msa import filter_seq_id
from .msa import compute_seq_identity


def _annot_to_color(attribute, tree, df_annot, cmap='jet'):
    """
    Reads in the attributes specified by the user in the annotation csv file
    and attributes a color palette for each.

    Arguments
    ---------
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
    color_dict = _get_color_palette(list(df_annot[attribute].unique()), cmap)
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
        tree_ete3, df_annot, sector_id, sector_seq, meta_data,
        fig_title='',
        t_sector_seq=True,
        t_sector_heatmap=True,
        colormap='inferno'
        ):
    """
    Update ete3 tree with sector info and attributes
    and return tree_style for further visualization.

    Arguments
    ---------
    tree_ete3 : ete3's tree object,
            as imported by io.load_tree_ete3()

    annot_file : pandas dataframe of the annotation file

    sector_id : list of sector identifiers, as imported by io.load_msa()
            the ids must match with the tree's leaves id

    sector_seq : corresponding list of sector sequences to display,
            as imported by io.load_msa()

    meta_data : tuple of annotations to display
                 (from annotation file's header)

    fig_title : figure title (str)

    t_sector_seq : boolean,
        whether to show the sequences of the sector

    t_sector_heatmap : boolean,
        whether to add a heatmap of the identity matrix between sector
        sequences

    Returns
    -------
    tree_style : TreeStyle class from ete3
    """

    leaves_id = tree_ete3.get_leaf_names()
    nb_leaves = len(leaves_id)

    tree_style = TreeStyle()
    tree_style.layout_fn = []

    # Add bootstrap support NodeStyle
    boot_style = NodeStyle()
    boot_style["fgcolor"] = "darkred"
    boot_style["size"] = 10
    empty_style = NodeStyle()
    empty_style["size"] = 0
    for node in tree_ete3.traverse():
        if node.support >= 95:
            node.set_style(boot_style)
        else:
            node.set_style(empty_style)

    column_layout = 0
    col_legend_rectface = 0
    if meta_data:  # if no meta_data, do nothing
        # add column associated with attributes
        def layout_attribute(node, column=column_layout):
            if node.is_leaf():
                name = node.name
                rect_faces = [None for i in range(len(meta_data))]
                for i in range(len(meta_data)):
                    colors, _ = _annot_to_color(meta_data[i],
                                                tree_ete3,
                                                df_annot)

                    rect_faces[i] = RectFace(50, 20,
                                             fgcolor=colors[name],
                                             bgcolor=colors[name],
                                             label=meta_data[i])
                    rect_faces[i].margin_left = 15
                    rect_faces[i].margin_right = 15
                    if i == len(meta_data) - 1:
                        rect_faces[i].margin_right = 30
                    add_face_to_node(rect_faces[i], node, column=column,
                                     position='aligned')
                    column += 1

        tree_style.layout_fn.append(layout_attribute)

        # Add legend
        legend_face = [None for i in range(len(meta_data))]
        for i in range(len(meta_data)):
            _, col_dict = _annot_to_color(meta_data[i], tree_ete3,
                                          df_annot)
            tree_style.legend.add_face(TextFace(meta_data[i],
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
    column_layout += len(meta_data)

    if t_sector_seq:
        sector_dict = {sector_id[i]: str(sector_seq[i])
                       for i in range(len(sector_id))}

        def layout_SeqMotifFace(node, column=column_layout):
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
        column_layout += 1

    if t_sector_heatmap:
        # allow to chose among Matplotlib's colormaps
        ProfileFace.get_color_gradient = _get_color_gradient
        # Check that sequences in the similarity matrix are ordered as in the
        # tree leaves and keep only sequences that are present in the tree
        reorder_msa = filter_seq_id(sector_id, sector_seq, leaves_id)
        id_mat = compute_seq_identity(reorder_msa[2])

        count = 0
        # Add heatmap profile to each leaf
        for lf in tree_ete3.iter_leaves():
            lf.add_features(profile=id_mat[count])
            count += 1
            lf.add_features(deviation=[0 for x in range(id_mat.shape[0])])
            lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5,
                                    width=(nb_leaves*20), height=20,
                                    style='heatmap',
                                    colorscheme=colormap),
                        column=column_layout, position="aligned")

    # Add title
    tree_style.title.add_face(TextFace(fig_title, fsize=20), column=0)

    return tree_style
