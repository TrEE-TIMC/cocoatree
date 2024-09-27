import numpy as np
from ._params import lett2num, __freq0


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


def aa_freq_at_pos(sequences, lbda=0.03, aa_count=21, weights=None):
    """Computes frequencies of aminoacids at each position of the alignment.
    .. math ::

    Arguments
    ----------
    sequences : list of sequences as imported by load_msa()

    lbda: regularization parameter lambda (default=0.03)

    aa_count: int, optional
            The number of amino acids to consider. Defaults to 21 to account
            for unknown (X) amino acids which are considered as gaps.

    weights: numpy 1D array, optional
            Gives more or less importance to certain sequences.
            If weights=None, all sequences are attributed an equal weight of 1.

    Returns
    -------
    fia : np.ndarray of shape (Npos, aa_count)
    """

    separated_aa = np.array(
        [np.array([char for char in row]) for row in sequences])
    N_seq = separated_aa.shape[0]
    if weights is None:
        weights = np.ones(N_seq)
    if lbda > 0:
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
    fia = []
    for pos in range(N_pos):
        tmp = np.bincount(separated_aa_num[:, pos], weights=weights,
                          minlength=21) / Neff_seq
        if lbda >= 0:
            tmp = (1 - lbda) * tmp + lbda / aa_count
        fia.append(tmp)
    fia = np.array(fia)

    return fia


def background_freq(fia, lbda=0.03, aa_count=21):
    """Computes regularized background frequencies of amino acids

    Arguments
    ---------
    fia : np.ndarray of the positional amino acid frequencies

    lbda : regularization parameter lambda (default=0.03)

    aa_count : int
            The number of amino acids to consider. Defaults to 21 to account
            for unknown amino acids which are considered as gaps.

    Returns
    -------
    qa : np.array of the background frequencies
    """

    # q0 : fraction of gaps in the alignment
    q0 = np.mean(fia[:, 0])
    qa = list((1 - q0) * __freq0)
    qa.insert(0, q0)
    qa = np.array(qa)

    if lbda > 0:
        qa = (1 - lbda) * qa + lbda / aa_count

    return qa


def aa_joint_freq(sequences, weights, lbda=0.03):
    """Computes the joint frequencies of each pair of amino acids in a MSA

    Arguments
    ----------
    sequences : list of sequences as imported by load_MSA()

    weights : sequence weights as calculated by seq_weights()

    lbda : regularization parameter lambda (default=0.03)

    Returns
    -------
    joint_freqs : amino acid joint frequencies

    joint_freqs_ind : amino acid joint frequencies if independent
    """

    # Convert sequences to binary format
    tmp = np.array([np.array([char for char in row]) for row in sequences])
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)

    aa_count, seq_nb, seq_length = binary_array.shape

    # Joint frequencies
    joint_freqs = np.zeros((seq_length, seq_length, aa_count, aa_count))

    # Frequencies if AAs are present independently at positions i & j
    joint_freqs_ind = np.zeros((seq_length, seq_length, aa_count, aa_count))

    # Adding weights
    weighted_binary_array = binary_array * weights[np.newaxis, :, np.newaxis]

    m_eff = np.sum(weights)
    simple_freq = weighted_binary_array / m_eff
    # Sum on the number of sequences
    simple_freq = np.sum(simple_freq, axis=1)

    simple_freq = (1 - lbda**2) * simple_freq + (lbda / aa_count)**2

    joint_freq_aibj = np.tensordot(weighted_binary_array, binary_array,
                                   axes=([1], [1])) / m_eff

    joint_freqs = joint_freq_aibj.transpose(1, 3, 0, 2)

    joint_freqs = (1 - lbda**2) * joint_freqs + lbda**2 / (aa_count)**2

    joint_freqs_ind = np.multiply.outer(simple_freq, simple_freq)

    joint_freqs_ind = np.moveaxis(joint_freqs_ind, [0, 1, 2, 3], [2, 0, 3, 1])

    return joint_freqs, joint_freqs_ind


def compute_sca_matrix(joint_freqs, joint_freqs_ind, fia, qa):
    """Compute SCA coevolution matrix

    Arguments
    ----------
    joint_freqs : amino acid joint frequencies

    joint_freqs_ind : amino acid joint frequencies if independent

    fia : frequency of amino acid *a* at position *i*

    qa : background frequency of amino acid *a*

    Returns
    -------
    Cijab : SCA coevolution matrix

    Cij : Frobenius norm of Cijab
    """

    Cijab_raw = joint_freqs - joint_freqs_ind

    # Derivee de l'entropie relative (fonction Phi)
    fia = fia.transpose([1, 0])
    phi = np.log(
        fia * (1 - qa[:, np.newaxis]) / (
            (1 - fia) * qa[:, np.newaxis])).transpose([1, 0])
    phi = np.multiply.outer(phi, phi).transpose([0, 2, 1, 3])

    Cijab_tilde = phi * Cijab_raw
    # Frobenius norm
    Cij = np.sqrt(np.sum(Cijab_tilde ** 2, axis=(2, 3)))

    return Cijab_raw, Cij
