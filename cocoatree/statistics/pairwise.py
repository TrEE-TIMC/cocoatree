import numpy as np
from ..__params import lett2num
from ..msa import compute_sequences_weights
from .position import _compute_aa_freq_at_pos


def _aa_joint_freq(sequences, weights, lambda_coef=0.03):
    """Computes the joint frequencies of each pair of amino acids in a MSA

    .. math::

        f_{ij}^{ab} = (1 - \\lambda) \\sum_s w_s \\frac{x_{si}^a x_{sj}^b}{M'}\
            + \\frac{\\lambda^2}{(21)^2}

    where

    .. math::

        M' = \\sum_s w_s

    represents the effective number of sequences in the alignment and *lambda*
    is a small regularization parameter (default=0.03).

    Arguments
    ----------
    sequences : list of sequences as imported by load_MSA()

    weights : ndarray of shape (Nseq)
            sequence weights as calculated by seq_weights()

    lambda_coef : float
        regularization parameter lambda (default=0.03)

    Returns
    -------
    joint_freqs : ndarray of shape (Nseq, Nseq)
                amino acid joint frequencies
    """

    # Convert sequences to binary format
    tmp = np.array([[char for char in row] for row in sequences])
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)

    aa_count, seq_nb, seq_length = binary_array.shape

    # Joint frequencies
    joint_freqs = np.zeros((seq_length, seq_length, aa_count, aa_count))

    # Frequencies if AAs are present independently at positions i & j

    # Adding weights
    weighted_binary_array = binary_array * weights[np.newaxis, :, np.newaxis]

    m_eff = np.sum(weights)
    simple_freq = weighted_binary_array / m_eff
    # Sum on the number of sequences
    simple_freq = np.sum(simple_freq, axis=1)

    simple_freq = (1 - lambda_coef**2) * simple_freq +\
        (lambda_coef / aa_count)**2

    joint_freq_aibj = np.tensordot(weighted_binary_array, binary_array,
                                   axes=([1], [1])) / m_eff

    joint_freqs = joint_freq_aibj.transpose(1, 3, 0, 2)

    joint_freqs = (1 - lambda_coef**2) * joint_freqs +\
        lambda_coef**2 / (aa_count)**2

    return joint_freqs


def compute_sca_matrix(joint_freqs, aa_freq, background_freq):
    """Compute the SCA coevolution matrix

    .. math::

        C_{ij}^{ab} = f_{ij}^{ab} - f_i^a f_j^b

    .. math::

        \\tilde{C_{ij}} = \\sqrt{sum_{a,b} \\tilde{(C_{ij}^{ab})^2}}

    Arguments
    ----------
    joint_freqs : amino acid joint frequencies

    aa_freq : frequency of amino acid *a* at position *i*

    background_freq : background frequency of amino acid *a*

    Returns
    -------
    Cijab : SCA coevolution matrix

    Cij : Frobenius norm of Cijab
    """
    joint_freqs_ind = np.multiply.outer(aa_freq, aa_freq)
    joint_freqs_ind = np.moveaxis(joint_freqs_ind, [0, 1, 2, 3],  [0, 2, 1, 3])

    Cijab_raw = joint_freqs - joint_freqs_ind

    # Derivee de l'entropie relative (fonction Phi)
    aa_freq = aa_freq.transpose([1, 0])
    phi = np.log(
        aa_freq * (1 - background_freq[:, np.newaxis]) / (
            (1 - aa_freq) * background_freq[:, np.newaxis])).transpose([1, 0])
    phi = np.multiply.outer(phi, phi).transpose([0, 2, 1, 3])

    Cijab_tilde = phi * Cijab_raw
    # Frobenius norm
    Cij = np.sqrt(np.sum(Cijab_tilde ** 2, axis=(2, 3)))

    return Cijab_raw, Cij


def compute_apc(MIij):
    """
    Computes the average product correction (APC) as described in Dunn et
    al. (2008).

    .. math::

        APC(a, b) = \\frac{MI(a, \\bar{x}) MI(b, \\bar{x}){\\overline{MI}}

    where :math:`MI(a, \\bar{x})` is the mean mutual information of column *a*
    and :math:`\\overline{MI}` is the overall mean mutual information

    The corrected mutual information is then:

    .. math::

        MIp(a, b) = MI(a, b) - APC(a, b)

    Arguments
    ----------
    MIij : np.ndarray,
        the mutual information matrix

    Returns
    -------
    APC_ij : np.ndarray,
        the average product correction (APC) matrix

    MIp : np.ndarray,
        the APC corrected mutual information matrix
    """

    n = MIij.shape[0]
    m = n - 1
    # Replace the matrix diagonal by 0
    np.fill_diagonal(MIij, 0)

    MI_colmean = (1/m) * np.sum(MIij, axis=0)
    MI_colmean = np.multiply.outer(MI_colmean, MI_colmean)

    MI_overmean = (2/(m*n)) * np.sum(np.tril(MIij))

    APC_ij = MI_colmean / MI_overmean

    MIp = MIij - APC_ij

    return APC_ij, MIp


def compute_entropy_correction(coevolution_matrix, s):

    """
    Computes the entropy correction according to Vorberg et al. (2018)

    .. math::

        C_{ij}^{EC} = C_{ij} - \\alpha s_{i}^{\\frac{1}{2}} \
            s_{j}^{\\frac{1}{2}}

    where :math:`\\alpha` is a coefficient determining the strength of the
    correction:

    .. math::

        \\alpha = \\frac{\\sum_{i \\neq j}^{L} c_ij \
        s_{i}^{\\frac{1}{2}}}{\\sum_{i \\neq j}^{L} s_i s_j}

    Arguments
    ---------
    coevolution_matrix : square matrix of shape (Nseq, Nseq)

    s : entropy computed for every position of the MSA

    Returns
    -------
    a square matrix of shape (Nseq, Nseq)
    """

    s_prod = np.multiply.outer(s, s)
    no_diag_eye = (1 - np.eye(s_prod.shape[0]))
    alpha = np.sum(
        (no_diag_eye * np.sqrt(s_prod) * coevolution_matrix) / np.sum(
            (no_diag_eye * s_prod)))

    return coevolution_matrix - alpha * np.sqrt(s_prod)


def compute_mutual_information_matrix(sequences, seq_weights=None,
                                      pseudo_count_val=0.03,
                                      normalize=True):
    """Compute the mutual information matrix

    .. math::

        I(X, Y) = \\sum_{x,y} p(x, y) \\log \\frac{p(x, y)}{p(x)p(y)}

    Arguments
    ----------
    sequences : list of sequences

    seq_weights : ndarray (nseq), optional, default: None
        if None, will compute sequence weights

    pseudo_count_val : float, default : 0.03
        Pseudo count value, to add to expected frequences (in order to have
        non-zero elements)

    normalize : boolean, default : True
        Whether to normalize the mutual information by the entropy.

    Returns
    -------
    mi_matrix : np.ndarray of shape (nseq, nseq)
        the matrix of mutual information
    """
    weights, _ = compute_sequences_weights(sequences)
    joint_freqs = _aa_joint_freq(
        sequences, weights, lambda_coef=pseudo_count_val)

    ind_freqs = _compute_aa_freq_at_pos(
        sequences, lambda_coef=pseudo_count_val,
        weights=weights)

    joint_freqs_ind = np.multiply.outer(ind_freqs, ind_freqs)
    joint_freqs_ind = np.moveaxis(joint_freqs_ind, [0, 1, 2, 3], [0, 2, 1, 3])

    mi_matrix = np.sum(
        joint_freqs * np.log(
            joint_freqs / joint_freqs_ind),
        axis=(2, 3))

    if normalize:
        joint_entropy = -np.sum(joint_freqs * np.log(joint_freqs), axis=(2, 3))
        mi_matrix /= joint_entropy

    return mi_matrix
