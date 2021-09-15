import numpy as np
from akasthesia.coevolution import Alignment
import scipy.special as sp

import types

from joblib import Parallel, delayed

aa = '-ACDEFGHIKLMNPQRSTVWY'


# Gibbs sampling
def conditional_prob(xt, res, x_single, x_pair):
    """Calculate conditional probability from v and w"""
    exclude_i = np.delete(np.arange(len(xt)), res)
    pot = (x_single[res] + np.sum(x_pair[res, exclude_i, :, xt[exclude_i]], axis=0))
    return np.exp((pot).T - sp.logsumexp(pot))


def gibbs_step(seq, x_single, x_pair, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    len_seq = seq.shape[0]
    _seq = np.copy(seq)
    for a in range(len_seq):
        cond_prob = conditional_prob(_seq, a, x_single, x_pair)
        _seq[a] = rng.choice(20, p=cond_prob)
    return _seq


def gibbs_sampling(init_seq, n_seq, x_single, x_pair, burnin=100, seed=42):
    """Gibbs sampling process

    Return a simulated MSA in numerical represenations, according to
    the model v and w using Gibbs sampling process

    Probability of encountering a sequence is (Volberg 2018)

    .. math:
        p(x|v,w) ~ exp(sum{i}(v_i(x_i)) + sum{ij}(w_ij(x_i, x_j)))

    Parameters
    ----------
    init_seq : initial sequences
    n_seq : int
        number of desired sequences
    x_single, x_pair : np array
        the statistical model
    n_steps : int
        number of gibbs steps per new accepted point

    Results
    -------
    Generator of sequences
        Result of simulation
    """
    seq = init_seq
    rng = np.random.default_rng(seed)
    # burn in phase
    for i in range(burnin):
        seq = gibbs_step(seq, x_single, x_pair, rng)
    # Sampling phase
    for i in range(n_seq):
        seq = gibbs_step(seq, x_single, x_pair, rng)
        yield seq


def gibbs_sampling_old(init_seq, n_seq, x_single, x_pair, n_steps, seed=42):
    """Gibbs sampling process
    #### OLD VERSION - checkout the new version above
    Return a simulated MSA in numerical represenations, according to
    the model v and w using Gibbs sampling process

    Parameters
    ----------
    init_seq : initial sequences
    n_seq : int
        number of desired sequences
    x_single, x_pair : np array
        the statistical model
    n_steps : int
        number of gibbs steps per new accepted point

    Results
    -------
    Generator of sequences
        Result of simulation
    """
    seq = init_seq
    rng = np.random.default_rng(seed)
    for i in range(n_seq):
        for __ in range(n_steps):
            seq = gibbs_step(seq, x_single, x_pair, rng)
        yield seq


def gibbs_sampling_p(init_seq, N, v, w, T, burnin=10, n_threads=2):
    """Parallelize version of Gibbs sampling
    UNDER DEVELOPMENT - DO NOT USE
    """
    return np.vstack(Parallel(n_jobs=n_threads)
                             (delayed(gibbs_sampling)
                              (init_seq, int(N/n_threads), v, w, T,
                              burnin=burnin) for _ in range(n_threads)))


# Support functions
def num_to_aa(num_seqs: np.ndarray):
    """Return a amino acid represenation of the MSA from numerical"""
    N = num_seqs.shape[0]
    seqs = []
    for i in range(N):
        seqs.append([aa[_+1] for _ in num_seqs[i]])
    return np.array(seqs)


def to_Alignment(seqs: np.ndarray):
    """Generate a set of headers for Cocoa's Alignment objects

    Parameters:
        seqs : nparray (NxL)
            amino acid representation of the MSA

    Returns
    -------
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
