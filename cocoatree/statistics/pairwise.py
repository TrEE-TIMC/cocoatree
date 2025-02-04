import numpy as np
from ..__params import lett2num, __freq_regularization_ref, __aa_count
from ..msa import compute_seq_weights
from .position import _compute_first_order_freqs


def _compute_aa_joint_freqs(sequences, seq_weights=None,
                            freq_regul=__freq_regularization_ref):
    """Computes the joint frequencies of each pair of amino acids in a MSA

    .. math::

        f_{ij}^{ab} = (\\sum_s w_s x_{si}^a x_{sj}^b +
            \\lambda/(21)^2)/(M_{eff} + \\lambda)

    where

    .. math::

        M_{eff} = \\sum_s w_s

    represents the effective number of sequences in the alignment and *lambda*
    is a regularization parameter (pseudocount).

    Arguments
    ----------
    sequences : list of sequences as imported by load_MSA()

    seq_weights : numpy 1D array, optional
            Gives more or less importance to certain sequences. If
            seq_weights=None, all sequences are attributed an equal weighti
            of 1.

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    aa_joint_freqs : np.ndarray of shape (Npos, Npos, aa_count, aa_count)
        joint frequency of amino acids `a` and `b`
        at respective positions `i` and `j`
    """

    # Convert sequences to binary format
    tmp = np.array([[char for char in row] for row in sequences])
    binary_array = np.array([tmp == aa for aa in lett2num.keys()]).astype(int)

    # Adding weights
    if seq_weights is None:
        seq_weights = np.ones(len(sequences))
    weighted_binary_array = binary_array * \
        seq_weights[np.newaxis, :, np.newaxis]
    # number of effective sequences
    m_eff = np.sum(seq_weights)

    # Joint frequencies
    aa_joint_freqs = np.tensordot(weighted_binary_array, binary_array,
                                  axes=([1], [1])).transpose(1, 3, 0, 2)
    aa_joint_freqs = (aa_joint_freqs + freq_regul * m_eff / __aa_count ** 2)\
        / ((1 + freq_regul) * m_eff)
    return aa_joint_freqs


def _compute_aa_product_freqs(aa_freqs_1, aa_freqs_2):
    """Computes the product of frequencies

    (joint frequencies if residues are independent)

    Arguments
    ----------
    aa_freqs_1 : frequency of amino acid *a* at position *i* (set 1)

    aa_freqs_2 : frequency of amino acid *a* at position *i* (set 2)

    Returns
    -------
    aa_prod_freqs : np.ndarray of shape (Npos, Npos, aa_count, aa_count)
        product of frequency of amino acids *a* and $b$
        at respective positions *i* and *j*
    """

    aa_product_freqs = np.multiply.outer(aa_freqs_1, aa_freqs_2)
    aa_product_freqs = np.moveaxis(aa_product_freqs,
                                   [0, 1, 2, 3],
                                   [0, 2, 1, 3])

    return aa_product_freqs


def _compute_second_order_freqs(sequences, seq_weights=None,
                                freq_regul=__freq_regularization_ref):
    """
    balabla
    """

    # joint frequencies
    aa_joint_freqs = _compute_aa_joint_freqs(sequences,
                                             seq_weights=seq_weights,
                                             freq_regul=freq_regul)

    aa_freqs, _ = _compute_first_order_freqs(
        sequences, seq_weights=seq_weights, freq_regul=freq_regul)

    # joint frequencies if independence (product of frequencies)
    aa_product_freqs = _compute_aa_product_freqs(aa_freqs, aa_freqs)

    return aa_joint_freqs, aa_product_freqs


def compute_sca_matrix(sequences, seq_weights=None,
                       freq_regul=__freq_regularization_ref):
    """Compute the SCA coevolution matrix

    .. math::

        C_{ij}^{ab} = f_{ij}^{ab} - f_i^a f_j^b

    .. math::

        \\tilde{C_{ij}} = \\sqrt{sum_{a,b} \\tilde{(C_{ij}^{ab})^2}}

    Arguments
    ----------
    sequences : list of sequences

    seq_weights : ndarray (nseq), optional, default: None
        if None, will compute sequence weights

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    SCA_matrix : SCA coevolution matrix
    """

    # computing frequencies
    if seq_weights is None:
        seq_weights, _ = compute_seq_weights(sequences)
    aa_joint_freqs, aa_product_freqs = _compute_second_order_freqs(
        sequences, seq_weights=seq_weights, freq_regul=freq_regul)

    # Cijab
    Cijab = aa_joint_freqs - aa_product_freqs

    # derivative of relative entropy
    aa_freqs, bkgd_freqs = _compute_first_order_freqs(
        sequences, seq_weights=seq_weights, freq_regul=freq_regul)
    aa_freqs = aa_freqs.transpose([1, 0])
    phi = np.log(
        aa_freqs * (1 - bkgd_freqs[:, np.newaxis]) / (
            (1 - aa_freqs) *
            bkgd_freqs[:, np.newaxis])).transpose([1, 0])
    phi = np.multiply.outer(phi, phi).transpose([0, 2, 1, 3])

    # applying sca positional weights
    Cijab_tilde = phi * Cijab

    # Frobenius norm
    SCA_matrix = np.sqrt(np.sum(Cijab_tilde ** 2, axis=(2, 3)))

    return SCA_matrix


def compute_mutual_information_matrix(sequences, seq_weights=None,
                                      freq_regul=__freq_regularization_ref,
                                      normalize=True):
    """Compute the mutual information matrix

    .. math::

        I(X, Y) = \\sum_{x,y} p(x, y) \\log \\frac{p(x, y)}{p(x)p(y)}

    Arguments
    ----------
    sequences : list of sequences

    seq_weights : ndarray (nseq), optional, default: None
        if None, will compute sequence weights

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    normalize : boolean, default : True
        Whether to normalize the mutual information by the entropy.

    Returns
    -------
    mi_matrix : np.ndarray of shape (nseq, nseq)
        the matrix of mutual information
    """

    # computing frequencies
    if seq_weights is None:
        seq_weights, _ = compute_seq_weights(sequences)
    aa_joint_freqs, aa_product_freqs = _compute_second_order_freqs(
        sequences, seq_weights=seq_weights,
        freq_regul=freq_regul)

    # mutual information
    mi_matrix = np.sum(
        aa_joint_freqs * np.log(aa_joint_freqs / aa_product_freqs),
        axis=(2, 3))

    if normalize:
        joint_entropy = -np.sum(aa_joint_freqs * np.log(aa_joint_freqs),
                                axis=(2, 3))
        mi_matrix /= joint_entropy

    return mi_matrix


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
