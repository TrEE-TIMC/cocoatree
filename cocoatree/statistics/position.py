import numpy as np
from ..__params import lett2num, __freq0


def aa_freq_at_pos(sequences, lambda_coef=0.03, weights=None):
    """Computes frequencies of aminoacids at each position of the alignment.

    .. math::
        f_i^a = (1 - \\lambda) \\sum_s w_s \\frac{x_{si}^a}{M^a} +\
             \\frac{\\lambda}{21}

    Arguments
    ----------
    sequences : list of sequences as imported by load_msa()

    lambda_coef : regularization parameter lambda (default=0.03)

    weights : numpy 1D array, optional
            Gives more or less importance to certain sequences.
            If weights=None, all sequences are attributed an equal weight of 1.

    Returns
    -------
    aa_freq : np.ndarray of shape (Npos, aa_count)
            frequency of amino acid *a* at position *i*
    """

    separated_aa = np.array([[char for char in row] for row in sequences])
    N_seq, N_pos = separated_aa.shape
    if weights is None:
        weights = np.ones(N_seq)
    if lambda_coef > 0:
        Neff_seq = np.sum(weights)
    else:
        Neff_seq = N_seq
    N_pos = separated_aa.shape[1]

    separated_aa_num = []
    # Convert amino acids to numerical values
    for seq in range(N_seq):
        num_aa = []
        for residue in separated_aa[seq]:
            code = lett2num[residue]
            num_aa.append(code)

        num_aa = np.array(num_aa)
        separated_aa_num.append(num_aa)
    # array of the amino acids as numericals
    separated_aa_num = np.array(separated_aa_num)

    # np.bincount : Count number of occurrences of each value in array of
    # non-negative ints.
    aa_freq = []
    # consider gaps as a 21st amino acid
    aa_count = 21
    for pos in range(N_pos):
        tmp = np.bincount(separated_aa_num[:, pos], weights=weights,
                          minlength=aa_count) / Neff_seq
        if lambda_coef >= 0:
            tmp = (1 - lambda_coef) * tmp + lambda_coef / aa_count
        aa_freq.append(tmp)
    aa_freq = np.array(aa_freq)

    return aa_freq


def compute_background_frequencies(aa_freq, lambda_coef=0.03):
    """Computes regularized background frequencies of amino acids

    Arguments
    ---------
    aa_freq : np.ndarray of the positional amino acid frequencies

    lambda_coef : regularization parameter lambda (default=0.03)

    Returns
    -------
    background_freq : np.array of the background frequencies
    """

    # q0 : fraction of gaps in the alignment
    q0 = np.mean(aa_freq[:, 0])
    # background_freq : correction factor on __freq0 in order to take the
    # proportion of gaps into account
    background_freq = list((1 - q0) * __freq0)
    background_freq.insert(0, q0)
    background_freq = np.array(background_freq)

    # consider gaps as 21st amino acid
    aa_count = 21
    if lambda_coef > 0:
        background_freq = (1 - lambda_coef) * background_freq + \
            lambda_coef / aa_count

    return background_freq


def compute_entropy(aa_freq):
    """Computes  Shannon's entropy for each position in the alignment

    .. math::

        H(a) = -\\sum_i f_{ia} \\log f_{ia}

    where *H(a)* is the relative entropy of amino acid *a*,
        *fia* is the frequency of amino acid *a* at position *i*

    Arguments
    ----------
    aa_freq: np.ndarray,
        amino acid frequencies per position

    Returns
    -------
    s: array of shape (N_pos)
    """

    s = -np.sum(aa_freq * np.log(aa_freq), axis=1)

    return s


def compute_rel_entropy(aa_freq, background_freq):
    """Kullback-Leibler relative entropy

    .. math::

        D_i^a = f_i^a \\ln \\frac{f_i^a}{q^a} + (1 - f_i^a) \\ln \
            \\frac{1 - f_i^a}{1 - q^a}

    where f_i^a is the observed frequency of amino acid *a* at position *i*,
        q^a is the background expectation
    D_i^a is known as the Kullback-Leibler relative entropy (Cover and Thomas,
    2012) and indicates how unlikely the observed frequencies of amino acid
    *a* at position *i* would be if *a* occurred randomly with probability q^a.

    Arguments
    ----------
    aa_freq: np.ndarray,
        amino acid frequencies per position

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

    Dia = aa_freq * np.log(aa_freq / background_freq) + \
        (1 - aa_freq) * np.log((1 - aa_freq) / (1 - background_freq))

    # sum on all amino acid at each position
    Di = np.sum(aa_freq * np.log(aa_freq / background_freq), axis=1)

    return Dia, Di
