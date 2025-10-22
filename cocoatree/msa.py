from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, substitution_matrices
import numpy as np
import sklearn.metrics as sn
from .__params import lett2num
from joblib import Parallel, delayed


def _clean_msa(msa):
    """
    This function compares the amino acid codes in the sequence alignment with
    the ones in lett2num and removes unknown amino acids (such as 'X' or 'B')
    when importing the multiple sequence alignment.

    Parameters
    ----------
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


def filter_sequences(sequences, sequences_id,
                     gap_threshold=0.4, seq_threshold=0.2,
                     verbose=False):
    """
    Filter sequences

    Remove (1) overly gapped positions; (2) overly gapped sequences.

    Parameters
    ----------
    sequences : list of MSA sequences to filter

    sequences_id : list of the MSA's sequence identifiers

    gap_threshold : float,
        maximum proportion of gaps tolerated per position (default=0.4)

    seq_threshold : float,
        maximum proportion of gaps tolerated per sequence (default=0.2)

    Returns
    -------
    filtered_seqs : list of the remaining sequences (written as strings)
        after applying the filters

    filtered_seqs_id : list of sequence identifiers that were kept after
        applying the filters

    remaining_pos : numpy.ndarray
        remaining positions after filtering
    """

    updated_sequences, remaining_pos = _filter_gap_pos(
        sequences, threshold=gap_threshold,
        verbose=verbose)
    filtered_seqs, filtered_seqs_id,  = _filter_gap_seq(
        updated_sequences, sequences_id,
        threshold=seq_threshold, verbose=verbose)

    return filtered_seqs, filtered_seqs_id, remaining_pos


def _filter_gap_pos(sequences, threshold=0.4, verbose=False):
    """Filter the sequences for overly gapped positions.

    Parameters
    ----------
    sequences : list of the MSA sequences to filter

    threshold : max proportion of gaps tolerated (default=0.4)

    Returns
    -------
    updated_seqs : updated_list of sequences with filtered gaps

    remaining_pos : numpy.ndarray
        remaining positions after filtering
    """

    if verbose:
        print("Filter MSA for overly gapped positions")

    Nseq, Npos = len(sequences), len(sequences[0])

    gaps = np.array([[int(sequences[seq][pos] == '-') for pos in range(Npos)]
                     for seq in range(Nseq)])

    freq_gap_per_pos = np.sum(gaps, axis=0) / Nseq

    remaining_pos = np.where(freq_gap_per_pos <= threshold)[0]

    if verbose:
        print("Keeping %i out of %i positions" % (len(remaining_pos), Npos))

    updated_seqs = ["".join([sequences[seq][pos] for pos in remaining_pos])
                    for seq in range(Nseq)]

    return updated_seqs, remaining_pos


def _filter_gap_seq(sequences, sequences_id, threshold=0.2, verbose=False):
    """
    Remove sequences with a fraction of gaps greater than a specified
    value.

    Parameters
    ----------
    sequences : list of MSA sequences

    sequences_id : list of the MSA's sequence identifiers

    threshold : maximum fraction of gaps per sequence (default 0.2)

    Returns
    -------
    filt_seqs : filtered list of sequences

    filt_seqs_id : corresponding list of sequence identifiers
    """

    if verbose:
        print('Filter MSA for overly gapped sequences')

    Nseq, Npos = len(sequences), len(sequences[0])

    freq_gap_per_seq = np.array([sequences[seq].count('-') / Npos
                                 for seq in range(Nseq)])

    filt_seqs_ix = np.where(freq_gap_per_seq <= threshold)[0]
    if verbose:
        print('Keeping %i sequences out of %i sequences' %
              (len(filt_seqs_ix), Nseq))

    filt_seqs = [sequences[seq] for seq in filt_seqs_ix]
    filt_seqs_id = [sequences_id[seq] for seq in filt_seqs_ix]

    return filt_seqs, filt_seqs_id


def filter_ref_seq(sequences, sequences_id, delta=0.2, refseq_id=None,
                   verbose=False):
    '''
    Filter the alignment based on identity with a reference sequence

    Remove sequences *r* with Sr < delta, where Sr is the fractional identity
    between the sequence *r* and a specified reference sequence.

    Parameters
    ----------
    sequences : list of sequences in the MSA

    sequences_id : list of sequence identifiers in the MSA

    delta : identity threshold (default=0.2)

    refseq_id : identifier of the reference sequence, if 'None', a reference
                sequence is choosen as the sequence that has the mean pairwise
                sequence identity closest to that of the entire sequence
                alignment (default 'None')

    Returns
    -------
    filt_seqs : filtered list of sequences

    filt_seqs_id : corresponding list of sequence identifiers
    '''

    Nseq = len(sequences)

    if refseq_id is None:
        if verbose:
            print('Choose a default reference sequence within the alignment')
        refseq_idx = _choose_ref_seq(sequences)
    else:
        if verbose:
            print('Reference sequence is: %i' % refseq_id)
        refseq_idx = sequences_id.index(refseq_id)

    sim_matrix = compute_seq_identity(sequences)
    filt_seqs_ix = np.where(sim_matrix[refseq_idx] >= delta)[0]
    filt_seqs = [sequences[seq] for seq in filt_seqs_ix]
    filt_seqs_id = [sequences_id[seq] for seq in filt_seqs_ix]

    if verbose:
        print('Keeping %i out of %i sequences' % (len(filt_seqs), Nseq))

    return filt_seqs, filt_seqs_id


def _choose_ref_seq(msa):
    """
    Determine a reference sequence for the alignment

    This function chooses a default reference sequence for the alignment by
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


def filter_seq_id(sequences, sequences_id, list_id):
    """
    Filter sequences based on list

    Filter a multiple sequence alignment to keep only sequences whose
    identifiers are in a user provided list.

    Parameters
    ----------
    sequences : list of MSA sequences

    sequences_id : list of the MSA's sequence identifiers

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
    for ident in sequences_id:
        if ident in list_id:
            new_record = SeqRecord(Seq(sequences[sequences_id.index(ident)]),
                                   id=ident)
            new_msa.append(new_record)

    seq_list = []
    id_list = []
    for record in new_msa:
        seq_list.append(str(record.seq))
        id_list.append(record.id)
    seq_list = np.array(seq_list)

    return [new_msa, id_list, seq_list]


def map_to_pdb(pdb_seq, pdb_pos, sequences, sequences_id, pdb_seq_id):
    """
    Mapping of the unfiltered MSA positions on a PDB structure.

    Parameters
    ----------
    pdb_seq: str,
        amino acid sequence of the reference PDB file

    pdb_pos: list,
        Residue positions as found in the PDB file

    sequences: list,
        List of sequences of the unfiltered MSA

    sequences_id: list,
        List of sequence identifiers in the unfiltered MSA

    pdb_seq_id: str,
        identifier of the sequence the positions are mapped onto. Should be
        included in sequences_id.

    Returns
    -------
    mapping: numpy.ndarray of shape (3, len(pdb_seq)),
        the first element is an array of the residues found in the PDB sequence
        the second element is an array of the PDB position of each amino acid
        the third element is an array of the positions of those same amino
        acids in the unfiltered MSA
    """
    msa_pos = []
    pdb_seq_idx = sequences_id.index(pdb_seq_id)
    for aa_index in range(len(sequences[pdb_seq_idx])):
        if sequences[pdb_seq_idx][aa_index] != '-':
            msa_pos.append(aa_index)

    mapping = np.array((list(pdb_seq), pdb_pos, msa_pos))

    return mapping


def compute_seq_identity(sequences):
    """
    Computes the identity between sequences in a MSA (as Hamming's pairwise
    distance)

    Parameters
    ----------
    sequences : list of sequences

    Returns
    -------
    sim_matrix : identity matrix of shape (Nseq, Nseq)
    """

    separated_aa = np.array([[lett2num[char] for char in row]
                             for row in sequences])

    sim_matrix = 1 - sn.DistanceMetric.get_metric(
        "hamming").pairwise(separated_aa)

    return sim_matrix


def compute_seq_weights(sequences, threshold=0.8, verbose_every=10000,
                        n_jobs=1, verbose_parallel=5):
    """
    Compute sequence weights

    Each sequence s is given a weight ws = 1/Ns where Ns is the number of
    sequences with an identity to s above a specified threshold.

    Parameters
    ----------
    sequences : list of sequences

    threshold : float, optional, default: 0.8
        percentage identity above which the sequences are considered identical
        (default=0.8)

    verbose_every : int
        if > 0, verbose every {verbose_every} sequences

    n_jobs : int, default=1 (no parallelization)
        the maximum number of concurrently running jobs
        (see joblib doc)

    verbose_parallel : int
        verbosity level for parallelization
        (see joblib doc)

    Returns
    -------
    weights : np.array (nseq, ) of each sequence weight

    m_eff : float
        number of effective sequences
    """
    if threshold < 0 or threshold > 1:
        raise ValueError(
            "The threshold needs to be between 0 and 1." +
            f" Value provided {threshold}")

    sequences_num = np.array([[lett2num[char] for char in row]
                              for row in sequences])

    if n_jobs == 1:
        seq_weights = []
        for iseq, seq in enumerate(sequences_num):
            if iseq % verbose_every == 0:
                print('computing weight of seq %d/%d\t' %
                      (iseq+1, len(sequences_num)), end='\r')
            sim = 1 - sn.DistanceMetric.get_metric(
                "hamming").pairwise([seq], sequences_num)
            seq_weights.append(1 / np.sum(sim >= threshold))
    else:
        def _weight_f(sequence, sequence_list):
            return 1 / np.sum(1 - sn.DistanceMetric.get_metric(
                "hamming").pairwise([sequence], sequence_list) >= threshold)

        seq_weights = Parallel(n_jobs=n_jobs, verbose=verbose_parallel)(
            delayed(_weight_f)(seq, sequences_num) for seq in sequences_num)

    seq_weights = np.array(seq_weights)
    m_eff = sum(seq_weights)

    return seq_weights, m_eff


def map_msa_positions(n_loaded_pos, remaining_pos):
    """
    Maps positions between the original and the filtered MSA

    Parameters
    ----------
    n_loaded_pos : int,
        Number of positions in the original unfiltered MSA

    remaining_pos : np.ndarray,
        array containing the indexes of positions that have been conserved
        after filtering the MSA (output from cocoatree.msa.filter_sequences)

    Returns
    -------
    original2filtered : dictionnary,
        the keys are the positions in the original MSA and the values are the
        corresponding positions in the filtered MSA. When the original position
        has been filtered, the value is set to 'None'.

    filtered2original : dictionnary,
        the keys are the positions in the filtered MSA and the values are the
        corresponding positions in the original MSA.
    """
    mapping = [
        int(val) if f else None
        for f, val in zip(
            np.isin(np.arange(n_loaded_pos), remaining_pos),
            np.isin(np.arange(n_loaded_pos), remaining_pos).cumsum()-1)]
    original2filtered = {
        i: t for i, t
        in enumerate(mapping)}

    filtered2original = dict()

    for pos in range(len(remaining_pos)):
        filtered2original[pos] = int(remaining_pos[pos])

    return original2filtered, filtered2original


def compute_seq_similarity(sequences, subst_matrix='BLOSUM62', gap_penalty=-4,
                           n_jobs=1, verbose_parallel=0):
    """
    Computes a similarity matrix using a precalculated substitution matrix.

    The similarity score for a pair of sequences is obtained as the sum of
    the substitution scores at each position of the sequence pair.

    Parameters
    ----------
    sequences : list of str,
        list of Nseq MSA sequences.
    subst_matrix : str, default='BLOSUM62'
        name of the substitution matrix.
        Type `Bio.Align.substitution_matrices.load()` to obtain a list of
        available substitution matrices.
    gap_penalty : int, default=-4
        penalty score for gaps. You can adjust this parameter to reflect
        biological assumptions (e.g., -1 for mild, -10 for harsh).
    n_jobs : int, default=1 (no parallelization)
        the maximum number of concurrently running jobs (-1 uses all
        available cores)
    verbose_parallel : int, default=0
        verbosity level for parallelization (see joblib doc)

    Returns
    -------
    similarity_matrix : np.ndarray,
        a (Nseq, Nseq) array of similarity scores.
    """
    matrix = substitution_matrices.load(subst_matrix)
    n = len(sequences)
    seq_length = len(sequences[0])

    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError("All sequences must be of equal length.")

    seq_array = np.array([list(seq) for seq in sequences])

    def score_pair(i, j):
        a_seq = seq_array[i]
        b_seq = seq_array[j]
        score = sum(
            0 if a == '-' and b == '-'
            else gap_penalty if a == '-' or b == '-'
            else matrix.get((a, b))
            for a, b in zip(a_seq, b_seq)
        )
        return i, j, score

    results = Parallel(n_jobs=n_jobs, verbose=verbose_parallel)(
        delayed(score_pair)(i, j) for i in range(n) for j in range(i, n)
    )

    similarity_matrix = np.zeros((n, n), dtype=int)
    for i, j, score in results:
        similarity_matrix[i, j] = score
        similarity_matrix[j, i] = score

    return similarity_matrix


def compute_normalized_seq_similarity(sequences, subst_matrix='BLOSUM62',
                                      gap_penalty=-4, n_jobs=1,
                                      verbose_parallel=0):
    """
    Computes a normalized similarity matrix using a precalculated substitution
    matrix.

    Each pairwise similarity score is normalized by the maximum possible
    score for the pair of sequences (i.e., the score we would obtain by
    comparing the sequence to itself).

    Parameters:
    -----------
    sequences : list of str,
        list of Nseq MSA sequences.
    subst_matrix : str, default='BLOSUM62'
        name of the substitution matrix.
        Type `Bio.Align.substitution_matrices.load()` to obtain a list of
        available substitution matrices.
    gap_penalty : int, default=-4
        Penalty score for gaps. You can adjust this parameter to reflect
        biological assumptions (e.g., -1 for mild, -10 for harsh).
    n_jobs : int, default=1 (no parallelization)
        the maximum number of concurrently running jobs (-1 uses all
        available cores)
    verbose_parallel : int, default=0
        verbosity level for parallelization (see joblib doc)

    Returns:
    --------
    similarity_matrix : np.ndarray,
        a (Nseq, Nseq) array of normalized similarity scores (0.0 to 1.0).
    """
    matrix = substitution_matrices.load(subst_matrix)
    n_seq = len(sequences)
    seq_length = len(sequences[0])

    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError("All sequences must be of equal length.")

    seq_array = np.array([list(seq) for seq in sequences])

    def score_pair(i, j):
        a_seq = seq_array[i]
        b_seq = seq_array[j]
        score = sum(
            0 if a == '-' and b == '-'
            else gap_penalty if a == '-' or b == '-'
            else matrix.get((a, b))
            for a, b in zip(a_seq, b_seq)
        )
        return i, j, score

    def max_score(i):
        seq = seq_array[i]
        return sum(
            0 if a == '-' else matrix.get((a, a), -1)
            for a in seq
        )

    # Compute maximum scores for normalization
    max_scores = Parallel(n_jobs=n_jobs, verbose=verbose_parallel)(
        delayed(max_score)(i) for i in range(n_seq)
    )

    # Compute pairwise scores
    results = Parallel(n_jobs=n_jobs, verbose=verbose_parallel)(
        delayed(score_pair)(i, j) for i in range(n_seq)
        for j in range(i, n_seq)
    )

    similarity_matrix = np.zeros((n_seq, n_seq), dtype=float)
    for i, j, score in results:
        max_possible = max(max_scores[i], max_scores[j])
        normalized = score / max_possible if max_possible != 0 else 0.0
        similarity_matrix[i, j] = normalized
        similarity_matrix[j, i] = normalized

    return similarity_matrix
