import numpy as np
from Bio import AlignIO


def MSA(filename, frmt, clean=False, verbose=True):
    """Import and formatting of multiple sequence alignment data

    Imports multiple sequence alignment data using BioPython's AlignIO

    Arguments
    ---------

    Returns
    ----------
    alignment : Bio.Align.MultipleSeqAlignment object

    id_list : list of sequence identifiers

    seq_list : list of sequences

    binary_array : MxLx20 binary array x_a_si where s = 1,..., M labels the
                   sequences, i = 1,..., L the positions, a = 1,..., 20 the
                   amino acids. x_a_si = 1 if sequence s has aa a at position
                   i and 0 otherwise.
    """

    alignment = AlignIO.read(filename, frmt)
    #print("Alignment of length %i" % alignment.get_alignment_length())

    # Option to clean the alignment
    if clean:
        alignment = clean_msa(alignment)

    seq_list = []
    id_list = []
    for record in alignment:
        #print(record.id)
        seq_list.append(str(record.seq))
        id_list.append(record.id)
    seq_list = np.array(seq_list)
    seq_list = np.char.replace(seq_list, 'X', '-')
    seq_list = np.char.replace(seq_list, 'B', '-')
    # Number of gaps
    #seq_list[1].count("-")

    tmp = np.array([np.array([char for char in row]) for row in seq_list])
    tmp = np.char.replace(tmp, 'X', '-', count = None)
    tmp = np.char.replace(tmp, 'B', '-', count = None)
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)

    if verbose is True:
        print("Alignment of length %i" % binary_array.shape[2])
        print("Number of sequences: %i" % binary_array.shape[1])

    return [alignment, id_list, seq_list, binary_array]
