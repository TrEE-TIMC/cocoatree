from . import position
from . import pairwise
from .. import msa
from ..__params import __freq_regularization_ref


def compute_all_frequencies(sequences,
                            seq_weights=None,
                            freq_regul=__freq_regularization_ref):
    """
    Compute frequencies on sequences


    Parameters
    ----------
    sequences : list of sequences

    seq_weights : {None, np.ndarray (n_seq)}
        if None, will re-compute the sequence weights.

    freq_regul : regularization parameter (default=__freq_regularization_ref)

    Returns
    -------
    aa_freqs : np.ndarray (nseq, 21)
        A (nseq, 21) ndarray containing the amino acid frequencies at each
        positions.

    bkgd_freqs :  np.ndarray (21, )
        A (21,) np.array containing the background amino acid frequencies
        at each position; it is computed from the mean frequency of amino acid
        a in all proteins in the NCBI non-redundant database
        (see Rivoire et al., https://dx.plos.org/10.1371/journal.pcbi.1004817)

    aa_joint_freqs : np.ndarray (nseq, nseq, 21, 21)
        An ndarray containing the pairwise joint frequencies of amino acids
        for each pair of positions in the list of provided sequences.
    """
    if seq_weights is None:
        seq_weights, _ = msa.compute_seq_weights(sequences)

    aa_freqs = position._compute_aa_freqs(
        sequences,
        freq_regul=freq_regul,
        seq_weights=seq_weights)

    bkgd_freqs = position._compute_background_freqs(
        aa_freqs,
        sequences,
        seq_weights=seq_weights,
        freq_regul=__freq_regularization_ref)

    aa_joint_freqs = pairwise._compute_aa_joint_freqs(
        sequences,
        seq_weights=seq_weights,
        freq_regul=freq_regul)

    return aa_freqs, bkgd_freqs, aa_joint_freqs
