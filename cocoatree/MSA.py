import numpy as np
from Bio import AlignIO
from Bio import Seq
from .__params import lett2num


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
        alignment = _clean_msa(alignment)

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



def _clean_msa(msa) :

    """Function to remove unknown amino acids when importing the multiple sequence alignment
    
    msa = bioalign object
    """

    for index, record in enumerate(msa) :
        #print(record.id) 
        for char in record.seq : 
            if char not in lett2num.keys() :
                #print('Sequence', record.id)
                #print(char, 'at position', record.seq.index(char))
                sequence = list(record.seq)
                sequence[record.seq.index(char)] = '-'
                sequence = "".join(sequence)
                #print(sequence)
                msa[index].seq = Seq(sequence)

    return msa
