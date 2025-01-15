from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from .__params import lett2num
import numpy as np
from .statistics.pairwise import compute_seq_identity


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

    if verbose:
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

    sim_matrix = compute_seq_identity(sequences, graphic=False)
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

    sim_matrix = compute_seq_identity(msa)

    mean_pairwise_seq_sim = np.mean(sim_matrix, axis=0)

    ref_seq = np.argmin(mean_pairwise_seq_sim)

    return ref_seq


def seq_weights(sim_matrix, threshold=0.8):
    """Each sequence s is given a weight ws = 1/Ns where Ns is the number of
    sequences with an identity to s above a specified threshold.

    Parameters
    ----------
    sim_matrix : similarity matrix (e.g. output from seq_similarity() function)

    threshold : percentage identity above which the sequences are considered
                identical (default=0.8)

    Returns
    -------
    weights : np.array of each sequence weight
    """

    weights = (1 / np.sum(sim_matrix >= threshold, axis=0))

    Nseq_eff = sum(weights)

    return weights, Nseq_eff


def filter_seq_id(seq_id, sequences, list_id):
    """
    Filter a multiple sequence alignment to keep only sequences whose
    identifiers are in a user provided list.

    Parameters
    ----------
    seq_id : list of the MSA's sequence identifiers

    sequences : list of MSA sequences

    list_id : list of sequence identifiers the user wants to keep. The
    identifiers must be in the same format as in the input MSA

    Returns
    -------
    new_msa : Bio.Align.MultipleSeqAlignment object,
            filtered msa

    id_list : list of sequence ID in the filtered MSA

    seq_list : list of sequences of the filtered MSA
    """
    new_msa = MultipleSeqAlignment([])
    new_record = SeqRecord([])
    for i in list_id:
        for ident in seq_id:
            if ident == i:
                new_record = SeqRecord(Seq(sequences[seq_id.index(ident)]),
                                       id=ident)
                new_msa.append(new_record)

    seq_list = []
    id_list = []
    for record in new_msa:
        seq_list.append(str(record.seq))
        id_list.append(record.id)
    seq_list = np.array(seq_list)

    return [new_msa, id_list, seq_list]
