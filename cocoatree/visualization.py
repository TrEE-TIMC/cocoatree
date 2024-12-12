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
from Bio import AlignIO
from PyQt5 import QtGui
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt

from .msa import filter_seq_id
from .statistics.pairwise import compute_seq_identity


def _annot_to_color(attribute, tree, df_annot, cmap='jet'):
    """
    Reads in the attributes specified by the user in the annotation csv file
    and attributes a color palette for each.

    Arguments
    ---------
    tree : ete3's tree object,
            as imported by io.load_tree()

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


def _get_sector_seq(sector_fasta, seq_id):
    """
    Get the amino acid sequences of the sector (fasta format)

    Arguments
    ---------
    sector_fasta : list of sequences as str as imported by io.load_msa()

    seq_id : list of sequence identifiers, as imported by io.load_msa()

    Returns
    -------
    sector : list of sequences

    sector_length : number of residues in the sector
    """
    sector = {}
    sector_length = len(sector_fasta[0])

    for i in range(0, len(sector_fasta)):
        sector[seq_id[i]] = str(sector_fasta[i])

    return sector, sector_length


# DON'T REVIEW YET, I STILL HAVE THINGS TO CHECK
# This function is very long, I was thinking maybe I should write smaller
# functions that create the different layouts and then they are each called
# by a wrapping function but I have to check with ete3 whether it is feasible
def plot_coev_along_phylogeny(tree, df_annot, sector_fasta, seq_id, attributes,
                              fig_title, rectface=True, seqmotif=True,
                              heatmap=True, colormap='inferno'):
    """
    Wrapping function that draws the phylogenetic tree along with specified
    sector characteristics.

    Arguments
    ---------
    tree : ete3's tree object,
            as imported by io.load_tree()

    annot_file : pandas dataframe of the annotation file

    sector_fasta : list of sequences to display, as imported by io.load_msa()

    seq_id : list of sequence identifiers, as imported by io.load_msa()
            the sequence names must match with the tree's sequence names

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

    id_lst = tree.get_leaf_names()
    nb_seq = len(id_lst)

    ts = TreeStyle()
    ts.layout_fn = []

    # Add bootstrap support NodeStyle
    boot_style = NodeStyle()
    boot_style["fgcolor"] = "darkred"
    boot_style["size"] = 10
    empty_style = NodeStyle()
    empty_style["size"] = 0
    for node in tree.traverse():
        if node.support >= 95:
            node.set_style(boot_style)
        else:
            node.set_style(empty_style)

    col_rectface = 0
    col_legend_rectface = 0
    if rectface:
        # Case with only one attribute to plot
        if isinstance(attributes, str):
            attribute_colors, col_dict = _annot_to_color(attributes, tree,
                                                         df_annot)

            def layout_RectFace(node):
                if node.is_leaf():
                    name = node.name
                    square = RectFace(50, 20, fgcolor=attribute_colors[name],
                                      bgcolor=attribute_colors[name])
                    square.margin_left = 10
                    add_face_to_node(square, node, column=col_rectface,
                                     position='aligned')
            ts.layout_fn.append(layout_RectFace)
            # Add legend
            ts.legend.add_face(TextFace(attributes, fsize=10, bold=True),
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
        # Case with several attributes to plot
        # Doesn't work because it writes a new layout function over the
        # previous one
        elif isinstance(attributes, list):
            for att in attributes:
                attribute_colors, col_dict = _annot_to_color(att, tree,
                                                             df_annot)

                def layout_RectFace(node):
                    if node.is_leaf():
                        name = node.name
                        square = RectFace(50, 20,
                                          fgcolor=attribute_colors[name],
                                          bgcolor=attribute_colors[name])
                        square.margin_left = 10
                        add_face_to_node(square, node, column=col_rectface,
                                         position='aligned')
                ts.layout_fn.append(layout_RectFace)
                col_rectface += 1
                # Add legend
                ts.legend.add_face(TextFace(att, fsize=10, bold=True),
                                   column=col_legend_rectface)
                # otherwise text is not in front of RectFace
                ts.legend.add_face(TextFace(""),
                                   column=col_legend_rectface + 1)
                for gene in col_dict.keys():
                    legend_face = RectFace(50, 20, fgcolor=col_dict[gene],
                                           bgcolor=col_dict[gene])
                    legend_face.margin_right = 5
                    ts.legend.add_face(legend_face, column=col_legend_rectface)
                    ts.legend.add_face(TextFace(gene, fsize=10),
                                       column=col_legend_rectface + 1)
                col_legend_rectface += 2

    col_seqmotif = 0
    if seqmotif:
        sector_seq, sector_length = _get_sector_seq(sector_fasta, seq_id)
        if rectface:
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
    if heatmap:
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
        for lf in tree.iter_leaves():
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

    return tree.show(tree_style=ts)
