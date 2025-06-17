"""Module to perform randomization of alignments"""

import numpy as np


def _randomize_seqs_conserving_col_compo(sequences=[], seed=None):
    """
    Randomize the list of sequenecs (MSA) so that the content of each
    colum is overall conserved (conservation of aa frequencies)

    Parameters
    ----------
    sequences : list of sequences (MSA)

    seed : int
        to generate exact same list of random numbers
        (mostly for testing )

    Returns
    -------
    rand_seqs : list of sequences where the columns have been shuffled
    """

    seq_array = np.array([list(seq) for seq in sequences])
    T = seq_array.T
    rng = np.random.default_rng(seed)
    rand_seq_array = np.array([rng.permutation(T[i]) for i in range(len(T))]).T
    rand_seqs = [''.join(seq) for seq in rand_seq_array]

    return rand_seqs
