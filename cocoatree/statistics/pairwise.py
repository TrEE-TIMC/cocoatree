import numpy as np
import sklearn.metrics as sn
from ..__params import lett2num


def compute_seq_identity(sequences):

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

    separated_aa = np.array([[lett2num[char] for char in row]
                             for row in sequences])

    sim_matrix = 1 - sn.DistanceMetric.get_metric(
        "hamming").pairwise(separated_aa)

    return sim_matrix


def aa_joint_freq(sequences, weights, lbda=0.03):
    """Computes the joint frequencies of each pair of amino acids in a MSA

    .. math::
        f_{ij}^{ab} = (1 - \lambda) \sum_s w_s \frac{x_{si}^a x_{sj}^b}{M'} + \frac{\lambda^2}{(21)^2}

    where

    .. math::
        M' = \sum_s w_s

    represents the effective number of sequences in the alignment and *lambda*
    is a small regularization parameter (default=0.03).

    Arguments
    ----------
    sequences : list of sequences as imported by load_MSA()

    weights : ndarray of shape (Nseq)
            sequence weights as calculated by seq_weights()

    lbda : float
        regularization parameter lambda (default=0.03)

    Returns
    -------
    joint_freqs : ndarray of shape (Nseq, Nseq)
                amino acid joint frequencies

    joint_freqs_ind : ndarray of shape (Nseq, Nseq)
                amino acid joint frequencies if independent
    """

    # Convert sequences to binary format
    tmp = np.array([[char for char in row] for row in sequences])
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


def compute_sca_matrix(joint_freqs, joint_freqs_ind, aa_freq, qa):
    """Compute the SCA coevolution matrix

    .. math::
        C_{ij}^{ab} = f_{ij}^{ab} - f_i^a f_j^b

    .. math::
        \tilde{C_{ij}} = \sqrt{sum_{a,b} \tilde{(C_{ij}^{ab})^2}}

    Arguments
    ----------
    joint_freqs : amino acid joint frequencies

    joint_freqs_ind : amino acid joint frequencies if independent

    aa_freq : frequency of amino acid *a* at position *i*

    qa : background frequency of amino acid *a*

    Returns
    -------
    Cijab : SCA coevolution matrix

    Cij : Frobenius norm of Cijab
    """

    Cijab_raw = joint_freqs - joint_freqs_ind

    # Derivee de l'entropie relative (fonction Phi)
    aa_freq = aa_freq.transpose([1, 0])
    phi = np.log(
        aa_freq * (1 - qa[:, np.newaxis]) / (
            (1 - aa_freq) * qa[:, np.newaxis])).transpose([1, 0])
    phi = np.multiply.outer(phi, phi).transpose([0, 2, 1, 3])

    Cijab_tilde = phi * Cijab_raw
    # Frobenius norm
    Cij = np.sqrt(np.sum(Cijab_tilde ** 2, axis=(2, 3)))

    return Cijab_raw, Cij
