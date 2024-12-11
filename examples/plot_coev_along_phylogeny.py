"""
==============================
Plot final figure
==============================

A small example that shows how to use the function
plot_coev_along_phylogeny of the visualization module.

"""

from cocoatree.visualization import plot_coev_along_phylogeny

annot_file = 'data/random_seq_tree_annotation.csv'
tree_file = 'data/random_seq_tree_kpicsg_ete3.treefile'
sector_fasta = 'data/SCA_sector_1_core_tagged_seq.fasta'

plot_coev_along_phylogeny(tree_file=tree_file, annot_file=annot_file,
                          sector_fasta=sector_fasta,
                          attributes=('HMM_annotation'),
                          fig_title='Test HMM annot',
                          rectface=True,
                          seqmotif=True,
                          heatmap=True,
                          colormap='inferno')
