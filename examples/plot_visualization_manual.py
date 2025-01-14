from cocoatree.io import load_MSA, load_tree
from cocoatree.msa import filter_seq_id
from cocoatree.statistics.pairwise import compute_seq_identity
from cocoatree.visualization import _annot_to_color, _get_sector_seq, \
    _get_color_gradient
import pandas as pd
from ete3 import TreeStyle, NodeStyle, RectFace, TextFace, add_face_to_node, \
    SeqMotifFace, ProfileFace

tree_file = 'data/random_seq_tree_kpicsg_ete3.treefile'
tree, id_lst = load_tree(tree_file)
annot_file = 'data/random_seq_tree_annotations.csv'
df_annot = pd.read_csv(annot_file)
sector_fasta = 'data/SCA_sector_1_core_tagged_seq.fasta'
sector_id, sector_seq1 = load_MSA(sector_fasta, 'fasta')

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


# RectFace
# First attribute
attributes = 'HMM_annotation'
col_HMM_annot = 0
col_legend_HMM_annot = 0
attribute_colors, col_dict = _annot_to_color(attributes, tree, df_annot)


def layout_HMM_annot(node):
    if node.is_leaf():
        name = node.name
        square = RectFace(50, 20, fgcolor=attribute_colors[name],
                          bgcolor=attribute_colors[name])
        square.margin_left = 10
        add_face_to_node(square, node, column=col_HMM_annot,
                         position='aligned')


ts.layout_fn.append(layout_HMM_annot)

# Add legend of first attribute
ts.legend.add_face(TextFace(attributes, fsize=10, bold=True),
                   column=col_legend_HMM_annot)
# otherwise text is not in front of RectFace
ts.legend.add_face(TextFace(""), column=col_legend_HMM_annot + 1)
for gene in col_dict.keys():
    legend_face = RectFace(50, 20, fgcolor=col_dict[gene],
                           bgcolor=col_dict[gene])
    legend_face.margin_right = 5
    ts.legend.add_face(legend_face, column=col_legend_HMM_annot)
    ts.legend.add_face(TextFace(gene, fsize=10),
                       column=col_legend_HMM_annot + 1)

# Second attribute
attributes = 'superkingdom'
col_superkingdom = col_HMM_annot + 1
col_legend_superkingdom = col_legend_HMM_annot + 2

attribute_colors, col_dict = _annot_to_color(attributes, tree, df_annot)


def layout_superkingdom(node):
    if node.is_leaf():
        name = node.name
        square = RectFace(50, 20, fgcolor=attribute_colors[name],
                          bgcolor=attribute_colors[name])
        square.margin_left = 10
        add_face_to_node(square, node, column=col_superkingdom,
                         position='aligned')


ts.layout_fn.append(layout_superkingdom)

# Add legend of second attribute
ts.legend.add_face(TextFace(attributes, fsize=10, bold=True),
                   column=col_legend_superkingdom)
ts.legend.add_face(TextFace(""), column=col_legend_superkingdom + 1)
for val in col_dict.keys():
    legend_face = RectFace(50, 20, fgcolor=col_dict[val],
                           bgcolor=col_dict[val])
    legend_face.margin_right = 5
    ts.legend.add_face(legend_face, column=col_legend_superkingdom)
    ts.legend.add_face(TextFace(val, fsize=10),
                       column=col_legend_superkingdom + 1)

# Third attribute: class
attributes = 'class'
col_class = col_superkingdom + 1
col_legend_class = col_legend_superkingdom + 2
attribute_colors, col_dict = _annot_to_color(attributes, tree, df_annot)


def layout_class(node):
    if node.is_leaf():
        name = node.name
        square = RectFace(50, 20, fgcolor=attribute_colors[name],
                          bgcolor=attribute_colors[name])
        square.amrgin_left = 10
        add_face_to_node(square, node, column=col_class, position='aligned')


ts.layout_fn.append(layout_class)

# Add legend of third attribute
ts.legend.add_face(TextFace(attributes, fsize=10, bold=True),
                   column=col_legend_class)
ts.legend.add_face(TextFace(""), column=col_legend_class + 1)
for val in col_dict.keys():
    legend_face = RectFace(50, 20, fgcolor=col_dict[val],
                           bgcolor=col_dict[val])
    legend_face.margin_right = 5
    ts.legend.add_face(legend_face, column=col_legend_class)
    ts.legend.add_face(TextFace(val, fsize=10), column=col_legend_class + 1)

# SeqMotifFace
sector_seq, sector_length = _get_sector_seq(sector_seq1, sector_id)
col_seqmotif = col_class + 1


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

# Heatmap
# allow to chose among Matplotlib's colormaps
ProfileFace.get_color_gradient = _get_color_gradient
# Check that sequences in the similarity matrix are ordered as in the
# tree leaves and keep only sequences that are present in the tree
reorder_msa = filter_seq_id(sector_seq1, sector_id, id_lst)
id_mat = compute_seq_identity(reorder_msa[2])

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
                            colorscheme='inferno'),
                column=col_heatmap, position="aligned")

ts.title.add_face(TextFace('Test HMM annot', fsize=20), column=0)

tree.show(tree_style=ts)
