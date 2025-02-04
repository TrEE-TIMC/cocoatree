import numpy as np
from .__params import __freq_regularization_ref
from cocoatree.randomize import _randomize_seqs_conserving_col_compo
from cocoatree.msa import compute_seq_weights
from cocoatree.statistics.pairwise import compute_sca_matrix
from cocoatree.pysca import _compute_ica, _icList


def extract_independent_components(sequences, coevo_matrix, method='pySCA',
                                   n_components=None, nrandom_pySCA=10,
                                   learnrate_ICA=0.1, nb_iterations_ICA=100000,
                                   freq_regul=__freq_regularization_ref,
                                   verbose_random_iter=True):
    """
    The current method is fully applicable to SCA analysis. For other metrics,
    we set n_components = 3 (to improve)

    Returns
    -------
    idpt_components : ndarray
        corresponding to a list of independent components
    """

    if n_components is None:
        if method == 'pySCA':
            n_components = _compute_n_components_as_pySCA(
                sequences, coevo_matrix,
                nrandom=nrandom_pySCA, freq_regul=freq_regul,
                verbose_random_iter=verbose_random_iter)
        else:
            n_components = 3

    V, S, Vt = np.linalg.svd(coevo_matrix)
    Vica, _ = _compute_ica(V, n_components,
                          learnrate=learnrate_ICA,
                          iterations=nb_iterations_ICA)

    idpt_components = Vica.T

    return idpt_components


def _compute_n_components_as_pySCA(sequences, coevo_matrix,
                                   seq_weights=None,
                                   nrandom=10,
                                   freq_regul=__freq_regularization_ref,
                                   verbose_random_iter=True):
    """
    Given the eigenvalues of the coevolution matrix, and the
    eigenvalues for the set of randomized matrices, return
    the number of significant eigenmodes as those above the average second
    eigenvalue plus 2 standard deviations.
    Based on S1 text of Rivoire et al. (2016)

    Rem: it concerns only SCA metrics
    For other merics (MI, adding corrections) this should be adapted    
    """

    if seq_weights is None:
        seq_weights, m_eff = compute_seq_weights(sequences)
    else:
        m_eff = np.sum(seq_weights)

    second_eigen_values_random = []
    for irand in range(nrandom):
        if verbose_random_iter:
            print('%d/%d randomized msa (to compute number of\
                  significant components) '
                  % (irand+1, nrandom), end='\r')
        rand_sequences = _randomize_seqs_conserving_col_compo(sequences)

        seq_weights_rand, m_eff_rand = compute_seq_weights(rand_sequences)

        # to get the correct m_eff
        seq_weights_rand = seq_weights_rand / m_eff_rand * m_eff

        SCA_rand = compute_sca_matrix(rand_sequences, seq_weights_rand,
                                      freq_regul=freq_regul)

        _, S, _ = np.linalg.svd(SCA_rand)
        second_eigen_values_random.append(S[1])

    mean_second_ev_rand = np.mean(second_eigen_values_random)
    std_second_ev_rand = np.std(second_eigen_values_random)

    _, S_input, _ = np.linalg.svd(coevo_matrix)

    n_components = len(S_input[S_input > mean_second_ev_rand +
                               2 * std_second_ev_rand])

    return n_components


def extract_principal_components(coevo_matrix):
    """
    """

    _, _, principal_components = np.linalg.svd(coevo_matrix)

    return principal_components


def extract_sectors(idpt_components, coevo_matrix):
    """
    """
    Vica = idpt_components.T
    _, sector_sizes, sorted_pos, _, _, _ = _icList(
        Vica, len(idpt_components), coevo_matrix)

    sectors = [[sorted_pos[i] for i in range(sector_sizes[0])]]
    ref_index = sector_sizes[0]
    for isize in range(1, len(sector_sizes)):
        sectors.append([sorted_pos[i]
                        for i in range(ref_index,
                                       ref_index + sector_sizes[isize])])
        ref_index += sector_sizes[isize]
    return sectors


def substract_first_principal_component(coevo_matrix):
    """
    In the sector literature (and data analysis), this corresponds
    to removing global correlations (from e.g. phylogenetic effects)
    """
    U, S, Vt = np.linalg.svd(coevo_matrix)
    S[0] = 0
    coevo_matrix_sub = np.maximum(np.linalg.multi_dot([U, np.diag(S), Vt]),0)

    return coevo_matrix_sub

