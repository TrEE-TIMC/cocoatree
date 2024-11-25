from Bio import AlignIO
from .msa import _clean_msa
from ete3 import Tree


def load_MSA(file_path, format, clean=True, verbose=False):
    """Read in a multiple sequence alignment (MSA)

    Arguments
    ---------
    file_path : path to the alignment file

    format : format of the alignment file (e.g. 'fasta', 'phylip', etc.)

    verbose : boolean,
            whether to print informations about the MSA

    Returns
    -------
    seq_id : list of sequence identifiers

    sequences : list of sequences as strings
    """

    alignment = AlignIO.read(file_path, format)

    if clean:
        alignment = _clean_msa(alignment)

    seq_id = list()
    sequences = list()
    for record in alignment:
        seq_id.append(record.id)
        sequences.append(str(record.seq))

    if verbose:
        print('Number of sequences: %i' % len(alignment))
        print('Alignment of length: %i' % len(alignment[0]))

    return seq_id, sequences


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
