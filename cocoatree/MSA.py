import numpy as np
from Bio import AlignIO
from Bio import Seq
from .__params import lett2num


def load_MSA(filename, frmt, clean=False, verbose=True):
    """Import and formatting of multiple sequence alignment data

    Imports multiple sequence alignment data using BioPython's AlignIO

    Arguments
    ---------
    filename : path to the alignment file

    frmt : format of the alignment file (e.g. 'fasta', 'phylip', etc.)

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

    # Option to clean the alignment
    if clean:
        alignment = _clean_msa(alignment)

    seq_list = []
    id_list = []
    for record in alignment:
        seq_list.append(str(record.seq))
        id_list.append(record.id)
    seq_list = np.array(seq_list)

    tmp = np.array([np.array([char for char in row]) for row in seq_list])
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)

    if verbose:
        print("Alignment of length %i" % binary_array.shape[2])
        print("Number of sequences: %i" % binary_array.shape[1])

    return [alignment, id_list, seq_list, binary_array]



def _clean_msa(msa) :

    """
    This function compares the amino acid codes in the sequence alignment with 
    the ones in lett2num and removes unknown amino acids (such as 'X' or 'B') 
    when importing the multiple sequence alignment.
    
    msa : bioalign object
    """

    for index, record in enumerate(msa) :
        for char in record.seq : 
            if char not in lett2num.keys() :
                sequence = list(record.seq)
                sequence[record.seq.index(char)] = '-'
                sequence = "".join(sequence)
                msa[index].seq = Seq(sequence)

    return msa
