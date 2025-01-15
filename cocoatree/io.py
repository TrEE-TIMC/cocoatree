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


def load_tree_ete3(file_path):
    """
    From the loading of a Newick tree, generate a ete3.Tree object

    Arguments
    ---------
    file_path : path to the Newick file

    Returns
    -------
    tree_ete3 : `ete3.Tree` object

    """
    tree_ete3 = Tree(file_path, format=0)
    return tree_ete3


def export_fasta(sequences, seq_id, outpath):
    """
    Function to export intermediate files in fasta format

    Arguments
    ---------
    sequences : list of sequences as strings (as imported by load_MSA)

    seq_id : list of sequences identifiers (as imported by load_MSA)

    outpath : path to the output file
    """

    # Add checks to see if the path exists?
    Nseq = len(sequences)
    with open(outpath, 'w') as outfile:
        for record in range(0, Nseq):
            outfile.write('>' + str(seq_id[record]) + '\n')
            outfile.write(str(sequences[record]) + '\n')
