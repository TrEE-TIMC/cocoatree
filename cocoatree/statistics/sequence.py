import sklearn.metrics as sn
from ..__params import lett2num
import numpy as np


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


def compute_sequences_weights(sequences, threshold=0.8):
    """
    Compute sequence weights

    Each sequence s is given a weight ws = 1/Ns where Ns is the number of
    sequences with an identity to s above a specified threshold.

    Parameters
    ----------
    sequences : list of sequences

    threshold : float, optional, default: 0.8

        percentage identity above which the sequences are considered identical
        (default=0.8)

    Returns
    -------
    weights : np.array (nseq, ) of each sequence weight

    n_seq_effective : float
        the number of effective sequences
    """
    if threshold < 0 or threshold > 1:
        raise ValueError(
            "The threshold needs to be between 0 and 1." +
            f" Value provided {threshold}")

    sim_matrix = compute_seq_identity(sequences)
    weights = (1 / np.sum(sim_matrix >= threshold, axis=0))

    n_seq_effective = sum(weights)

    return weights, n_seq_effective
