"""
Generate a simulated MSA using model by (Volberg et al, 2018)
"""


import numpy as np
from akasthesia.coevolution import Alignment
import scipy.special as sp

import types

from joblib import Parallel, delayed

aa = '-ACDEFGHIKLMNPQRSTVWY'


# Gibbs sampling
def conditional_prob(xt, res, x_single, x_pair):
    """Calculate conditional probability from v and w

    .. math::
        p(x^{t+1} = a|x^t, v, w) = \\frac{1}{Z}\\exp (v_i(a) + \\sum_{j\\neq i} w_{i, j}(a, x^t_j))

    Parameters
    ----------
        xt : numpy.ndarray
            the current sequence/state
        res : int
            the residue to calculate the conditional prob
        x_single, x_pair : numpy.ndarray
            the statistical model

    Return:
        numpy.ndarray :
            Conditional probability for 20 amino aicds
    """
    exclude_i = np.delete(np.arange(len(xt)), res)
    pot = (x_single[res] + np.sum(x_pair[res, exclude_i, :, xt[exclude_i]], axis=0))
    return np.exp((pot).T - sp.logsumexp(pot))


def gibbs_step(seq, x_single, x_pair, rng=None):
    """
    Single Gibbs sampling step

    Update the sequence according to the current state (Markov Chain) and the model

    Parameters
    ----------
        seq: numpy.ndarray
            the current sequence/state
        x_single, x_pair: numpy.ndarray
            the statistical model
        rng: numpy.random.Generator, optional
            the random generator to generate the output

    Return:
        numpy.ndarray :
            The updated sequence

    """
    if rng is None:
        rng = np.random.default_rng()
    len_seq = seq.shape[0]
    _seq = np.copy(seq)
    for a in range(len_seq):
        cond_prob = conditional_prob(_seq, a, x_single, x_pair)
        _seq[a] = rng.choice(20, p=cond_prob)
    return _seq


def gibbs_sampling(init_seq, n_seq, x_single, x_pair, n_steps, burnin=100, seed=42):
    """Gibbs sampling process

    Return a simulated MSA in numerical represenations, according to
    the model v and w using Gibbs sampling process

    The probability of encounting a sequence x:
    
    .. math::
        p(x|v,w) \\propto (\\sum_i v_i + \\sum_{i,j} w_{i,j})

    Parameters
    ----------
        init_seq: numpy.ndarray
            initial sequences
        n_seq: int
            number of desired sequences
        x_single, x_pair: numpy.ndarray
            the statistical model
        n_steps: int
            number of gibbs steps between new accepted sequence
        burnin: int
            number of burnin (throwaway) Gibbs steps before start sampling

    Return:
        generator:
            Result of simulation
    """
    seq = init_seq
    rng = np.random.default_rng(seed)
    for i in range(burnin):
        seq = gibbs_step(seq, x_single, x_pair, rng)
    for i in range(n_seq):
        for __ in range(n_steps):
            seq = gibbs_step(seq, x_single, x_pair, rng)
        yield seq


# Support functions
def num_to_aa(num_seqs: np.ndarray):
    """Return a amino acid represenation of the MSA from numerical

    Parameters
    ----------
        num_seqs: numpy.ndarray
            Numerical representation of the MSA

    Returns:
        numpy.ndarray :
            Character representation of the MSA
    """
    N = num_seqs.shape[0]
    seqs = []
    for i in range(N):
        seqs.append([aa[_+1] for _ in num_seqs[i]])
    return np.array(seqs)


def to_Alignment(seqs: np.ndarray):
    """Generate a set of headers for Cocoa's Alignment objects

    Parameters
    ----------
        seqs: numpy.ndarray (NxL)
            amino acid representation of the MSA

    Returns:
        Alignment :
            Alignment object of the MSAs
    """
    if isinstance(seqs, types.GeneratorType):
        seqs = np.array(list(seqs))
    if seqs.dtype is np.dtype(np.int_):
        seqs = num_to_aa(seqs)
    N = seqs.shape[0]
    headers = []
    for i in range(N):
        headers.append(" ".join(["Generated sequence No. ", str(i)]))
    return Alignment(headers, seqs, 1)


def recompute_v(background_freq, x_pair, z=1, lmbd=1):
    """
    Reformulation of the mathematic model

    Take the positional frequency of amino acids at positions and readjust v (x_single) so that the
    outcome frequency reflect the input frequency

    .. math::
        v_i(a) = log(f_i(a)) + log Z - \\lambda \\sum_c \\sum_{j \\neq i} w_{ij}(a, c)

    Parameters
    ----------
        background_freq : numpy.ndarray
            Desired outcome frequency of amino acid at residues
        x_pair : numpy.ndarray
            Pairwise statistical potential
        z : int, optional
            A constant to adjust the value of z; default: 1
        lmbd : int, optional
            A constant to adjust the impact of x_pair; default: 1

    Returns
    -------
        numpy.ndarray :
            Recomputed v (single residue potential)
    """
    return (np.tile(np.log(z), (20, 1)).T
            + np.log(background_freq)
            - lmbd*np.sum(x_pair, axis=(1, 3)))
