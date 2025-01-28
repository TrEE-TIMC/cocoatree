import numpy as np


def compute_seq_weights(sim_matrix, threshold=0.8):
    """
    Compute sequence weights

    Each sequence s is given a weight ws = 1/Ns where Ns is the number of
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
