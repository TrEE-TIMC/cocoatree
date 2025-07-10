import numpy as np
from ..__params import lett2num, __freq_regularization_ref, __aa_count, __freq0
from ..msa import compute_seq_weights


def _compute_aa_freqs(sequences, seq_weights=None,
                      freq_regul=__freq_regularization_ref):
    """Computes frequencies of amino acids at each position of the alignment.

    .. math::
        f_i^a = (\\sum_s w_s x_{si}^a + \\lambda/21)/(M_{eff} + \\lambda)

    where

    .. math::

        M_{eff} = \\sum_s w_s

    represents the effective number of sequences in the alignment and *lambda*
    is a regularization parameter (pseudocount).

    Parameters
    ----------
    sequences : list of sequences as imported by load_msa()

    seq_weights : numpy 1D array, optional
            Gives more or less importance to certain sequences. If
            seq_weights=None, all sequences are attributed an equal weight
            of 1.

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    aa_freq : np.ndarray of shape (Npos, aa_count)
            frequency of amino acid *a* at position *i*
    """

    tmp = np.array([[char for char in row] for row in sequences])
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)
    # weights
    if seq_weights is None:
        seq_weights = np.ones(len(sequences))
    m_eff = np.sum(seq_weights)
    weighted_binary_array = \
        binary_array * seq_weights[np.newaxis, :, np.newaxis]
    aa_freq = (np.sum(weighted_binary_array, axis=1).T
               + freq_regul * m_eff / __aa_count) / ((1 + freq_regul) * m_eff)

    return aa_freq


def _compute_background_freqs(aa_freqs, sequences, seq_weights=None,
                              freq_regul=__freq_regularization_ref):
    """Computes (regularized) background frequencies of amino acids

    Parameters
    ----------
    aa_freqs : np.ndarray of the positional amino acid frequencies

    sequences : list of sequences for which seq_weights give weights

    seq_weights : numpy 1D array, optional
            Gives more or less importance to certain sequences.
            If seq_weights=None, all sequences are attributed an equal weight
            of 1.

    freq_regul : regularization parameter (default=__freq_regularization_ref)


    Returns
    -------
    bkgd_freqs :  np.ndarray (21, )
        A (21,) np.array containing the background amino acid frequencies
        at each position; it is computed from the mean frequency of amino acid
        *a* in all proteins in the NCBI non-redundant database
        (see Rivoire et al., https://dx.plos.org/10.1371/journal.pcbi.1004817)
    """

    # q0 : fraction of gaps in the alignment
    q0 = np.mean(aa_freqs[:, 0])
    # background_freq : correction factor on __freq0 in order to take the
    # proportion of gaps into account
    bkgd_freqs = list((1 - q0) * __freq0)
    bkgd_freqs.insert(0, q0)
    bkgd_freqs = np.array(bkgd_freqs)

    # weights
    if seq_weights is None:
        seq_weights = np.ones(len(sequences))
    m_eff = np.sum(seq_weights)

    # regularization
    bkgd_freqs = (bkgd_freqs * m_eff +
                  freq_regul * m_eff / __aa_count) / ((1 + freq_regul) * m_eff)

    return bkgd_freqs


def _compute_first_order_freqs(sequences, seq_weights=None,
                               freq_regul=__freq_regularization_ref):
    """
    Compute amino acid frequencies at each position and background frequencies

    Parameters
    ----------
    sequences : list of sequences for which seq_weights gives weights

    seq_weights : numpy 1D array, optional, default=None
        Gives more or less importance to certain sequences.
        If seq_weights=None, will compute sequence weights

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    aa_freqs : np.ndarray of the positional amino acid frequencies

    bkgd_freqs : np.ndarray (21, )
    A (21, ) np.array containing the background amino acid frequencies at each
    position. It is computed from the mean frequency of amino acid *a* in all
    proteins in the NCBI non-redundant database.

    (see Rivoire et al., https://dx.plos.org/10.1371/journal.pcbi.1004817)
    """

    if seq_weights is None:
        seq_weights, _ = compute_seq_weights(sequences)

    aa_freqs = _compute_aa_freqs(
        sequences,
        freq_regul=freq_regul,
        seq_weights=seq_weights)

    bkgd_freqs = _compute_background_freqs(
        aa_freqs,
        sequences,
        seq_weights=seq_weights,
        freq_regul=__freq_regularization_ref)

    return aa_freqs, bkgd_freqs


def compute_entropy(aa_freq):
    """Computes  Shannon's entropy for each position in the alignment

    .. math::

        H(a) = -\\sum_i f_{ia} \\log f_{ia}

    where *H(a)* is the relative entropy of amino acid *a*,
        *fia* is the frequency of amino acid *a* at position *i*

    Parameters
    ----------
    aa_freq : np.ndarray,
        amino acid frequencies per position

    Returns
    -------
    s: array of shape (N_pos)
    """

    s = -np.sum(aa_freq * np.log(aa_freq), axis=1)

    return s


def compute_conservation(sequences, seq_weights=None,
                         freq_regul=__freq_regularization_ref):
    """
    Compute the conservation of amino acid at each position.

    The conservation is computed as the relative entropy (e.g., the
    Kullback-Leibler divergence)

    .. math::

        D_i^a = f_i^a \\ln \\frac{f_i^a}{q^a} + (1 - f_i^a) \\ln \
            \\frac{1 - f_i^a}{1 - q^a}

    where :math:`f_i^a` is the observed frequency of amino acid `a` at
        position i`, :math:`q^a` is the background expectation

    :math:`D_i^a` indicates how unlikely the observed frequencies of amino
    acid `a` at position `i` would be if `a` occurred randomly with
    probability :math:`q^a`.

    Parameters
    ----------
    sequences : list of sequences

    seq_weights : ndarray (nseq), optional, default: None
        if None, will compute sequence weights

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    Di : np.ndarray (npos,)
        where each entry corresponds to the conservation at this position in
        the sequences.

    """

    aa_freqs, bkgd_freqs = _compute_first_order_freqs(sequences, seq_weights,
                                                      freq_regul)

    _, Di = _compute_rel_entropy(aa_freqs, bkgd_freqs)

    return Di


def _compute_rel_entropy(aa_freqs, bkgd_freqs):
    """Compute the relative entropy

    Also know as the Kullback-Leibler divergence

    .. math::

        D_i^a = f_i^a \\ln \\frac{f_i^a}{q^a} + (1 - f_i^a) \\ln \
            \\frac{1 - f_i^a}{1 - q^a}

    where f_i^a is the observed frequency of amino acid *a* at position *i*,
        q^a is the background expectation

    D_i^a is known as the Kullback-Leibler relative entropy (Cover and Thomas,
    2012) and indicates how unlikely the observed frequencies of amino acid
    *a* at position *i* would be if *a* occurred randomly with probability q^a.

    Parameters
    ----------
    aa_freqs: np.ndarray,
        amino acid frequencies per position

    bck_freq: np.ndarray,
        background frequenvies of amino acids

    returns
    -------
    Dia: np.ndarray,
        relative entropy of aa_freq given the background distribution of amino
        acids. Indicates how unlikely the observed frequency of amino acid *a*
        at position *i* would be if a occurred randomly with probability
        background_freq

    Di: np.ndarray,
        overall conservation of position *i* taking all amino acids into
        account
    """

    Dia = aa_freqs * np.log(aa_freqs / bkgd_freqs) + \
        (1 - aa_freqs) * np.log((1 - aa_freqs) / (1 - bkgd_freqs))

    # sum on all amino acid at each position
    Di = np.sum(aa_freqs * np.log(aa_freqs / bkgd_freqs), axis=1)

    return Dia, Di
