from . import position
from . import pairwise
from .. import msa


def compute_all_frequencies(sequences,
                            seq_weights=None,
                            pseudo_counts_val=0.03):
    """
    Compute frequencies on sequences


    Parameters
    ----------
    sequences : list of sequences

    seq_weights : {None, np.ndarray (n_seq)}
        if None, will re-compute the sequence weights.

    pseudo_counts_val : float, default : 0.03
        pseudo counts value in order to avoid null entries in the estimation
        of the amino acid frequencies.

    Returns
    -------
    aa_frequencies : np.ndarray (nseq, 21)
        A (nseq, 21) ndarray containing the amino acid frequencies at each
        positions.

    background_frequencies : np.ndarray (21, )
        A (21,) vector containing the average use of amino acid in the list of
        provided sequences.

    pairwise_frequencies : np.ndarray (nseq, nseq, 21, 21)
        An ndarray containing the pairwise frequencies for each pairs of
        positions in the list of provided sequences.
    """
    if seq_weights is None:
        seq_weights, _ = msa.compute_sequences_weights(sequences)

    aa_frequencies = position._compute_aa_freq_at_pos(
        sequences,
        lambda_coef=pseudo_counts_val,
        weights=seq_weights)

    background_frequencies = \
        position._compute_regularized_background_frequencies(
            aa_frequencies)

    pairwise_frequencies = pairwise._aa_joint_freq(
        sequences,
        weights=seq_weights,
        lambda_coef=pseudo_counts_val)

    return aa_frequencies, background_frequencies, pairwise_frequencies
