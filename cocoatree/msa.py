import numpy as np
from Bio import AlignIO
from Bio import Seq
from .__params import lett2num
import sklearn.metrics as sn


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


def _clean_msa(msa):
    """
    This function compares the amino acid codes in the sequence alignment with
    the ones in lett2num and removes unknown amino acids (such as 'X' or 'B')
    when importing the multiple sequence alignment.

    Arguments
    ---------
    msa : bioalign object
    """

    for index, record in enumerate(msa):
        for char in record.seq:
            if char not in lett2num.keys():
                sequence = list(record.seq)
                sequence[record.seq.index(char)] = '-'
                sequence = "".join(sequence)
                msa[index].seq = Seq(sequence)

    return msa


def filter_gap_pos(sequences, threshold=0.4, verbose=False):
    """Filter the sequences for overly gapped positions.

    Arguments
    ---------
    sequences : list of the MSA sequences to filter

    threshold : max proportion of gaps tolerated (default=0.4)

    Returns
    -------
    filt_seqs : list of the sequences after filter

    pos_kept : numpy.ndarray of the positions that were conserved
    """

    print("Filter MSA for overly gapped positions")

    Nseq, Npos = len(sequences), len(sequences[0])

    gaps = np.array([[int(sequences[seq][pos] == '-') for pos in range(Npos)]
                     for seq in range(Nseq)])

    freq_gap_per_pos = np.sum(gaps, axis=0) / Nseq

    pos_kept = np.where(freq_gap_per_pos <= threshold)[0]

    if verbose:
        print("Keeping %i out of %i positions" % (len(pos_kept), Npos))

    filt_seqs = ["".join([sequences[seq][pos] for pos in pos_kept])
                 for seq in range(Nseq)]

    return filt_seqs, pos_kept


def filter_gap_seq(seq_id, sequences, threshold=0.2, filtrefseq=False,
                   refseq_id=None, verbose=False):
    """
    Remove sequences with a fraction of gaps greater than a specified
    value.
    Also possibility to remove sequences with sequence identity too high
    with a given reference sequence.

    Arguments
    ---------
    seq_id : list of the MSA's sequence identifiers

    sequences : list of MSA sequences

    threshold : maximum fraction of gaps per sequence (default 0.2)

    filtrefseq : boolean, whether to filter based on a reference sequence
                 (default False)

    refseq_id : str, default = None
                identifier of the reference sequence (only if filtrefseq=True)
                If 'None', a default reference sequence is chosen

    Returns
    -------
    seq_id_kept : list of conserved sequence identifiers

    seq_kept : list of kept sequences
    """

    if verbose:
        print('Filter MSA for overly gapped sequences')

    Nseq, Npos = len(sequences), len(sequences[0])

    freq_gap_per_seq = np.array([sequences[seq].count('-') / Npos
                                 for seq in range(Nseq)])

    seq_kept_index = np.where(freq_gap_per_seq <= threshold)[0]
    if verbose:
        print('Keeping %i sequences out of %i sequences' %
              (len(seq_kept_index), Nseq))

    seq_kept = [sequences[seq] for seq in seq_kept_index]
    seq_id_kept = [seq_id[seq] for seq in seq_kept_index]

    if filtrefseq:
        if verbose:
            print('Remove sequences too similar to the reference sequence')
        seq_id_kept, seq_kept = filter_ref_seq(seq_id_kept, seq_kept,
                                               delta=0.2, refseq_id=refseq_id)

    return seq_id_kept, seq_kept


def filter_ref_seq(seq_id, sequences, delta=0.2, refseq_id=None,
                   verbose=False):
    '''
    Remove sequences r with Sr < delta, where Sr is the fractional identity
    between r and a specified reference sequence.

    Arguments
    ---------
    seq_id : list of sequence identifiers in the MSA

    sequences : list of sequences in the MSA

    delta : identity threshold (default 0.2)

    refseq_id : identifier of the reference sequence, if 'None', a reference
                sequence is computed (default 'None')

    Returns
    -------
    seq_id_kept : list of identifiers of the kept sequences

    seq_kept : list of the kept sequences
    '''

    Nseq = len(sequences)

    if refseq_id is None:
        if verbose:
            print('Choose a default reference sequence within the alignment')
        refseq_idx = choose_ref_seq(sequences)
    else:
        if verbose:
            print('Reference sequence is: %i' % refseq_id)
        refseq_idx = seq_id.index(refseq_id)

    sim_matrix = seq_identity(sequences, graphic=False)
    seq_kept_index = np.where(sim_matrix[refseq_idx] >= delta)[0]
    seq_kept = [sequences[seq] for seq in seq_kept_index]
    seq_id_kept = [seq_id[seq] for seq in seq_kept_index]

    if verbose:
        print('Keeping %i out of %i sequences' % (len(seq_kept), Nseq))

    return seq_id_kept, seq_kept


def choose_ref_seq(msa):

    """This function chooses a default reference sequence for the alignment by
    taking the sequence which has the mean pairwise sequence identity closest
    to that of the entire sequence alignment.

    Parameters
    ----------
    msa : the multiple sequence alignment as a list of sequences

    Returns
    -------
    The index of the reference sequence in the given alignment
    """

    sim_matrix = seq_identity(msa)

    global_sim = np.mean(sim_matrix)

    mean_pairwise_seq_sim = np.mean(sim_matrix, axis=0)

    meanDiff = np.abs(mean_pairwise_seq_sim - global_sim)

    ref_seq = np.argmin(meanDiff)

    return ref_seq


def seq_identity(sequences):

    """
    Computes the identity between sequences in a MSA (as Hamming's pairwise
    distance)

    Arguments
    ---------
    sequences : list of sequences

    Returns
    -------
    sim_matrix : identity matrix of shape (Nseq, Nseq)
    """

    separated_aa = np.array(
        [np.array([lett2num[char] for char in row])
         for row in sequences])

    sim_matrix = 1 - sn.DistanceMetric.get_metric(
        "hamming").pairwise(separated_aa)

    return sim_matrix
