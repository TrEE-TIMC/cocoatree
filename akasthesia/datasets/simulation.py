import numpy as np
from akasthesia.coevolution import Alignment
import scipy.special as sp

import types

from joblib import Parallel, delayed

aa = '-ACDEFGHIKLMNPQRSTVWY'


# Gibbs sampling
def conditional_prob(xt, res, v, w):
    """Calculate conditional probability from v and w"""
    j = np.delete(np.arange(len(xt)), res)
    pot = (v[res] + np.sum(w[res, j, :, xt[j]], axis=0))
    return np.exp((pot).T - sp.logsumexp(pot))


def gibbs_step(seq, v, w, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    len_seq = seq.shape[0]
    _seq = np.copy(seq)
    for a in range(len_seq):
        cond_prob = conditional_prob(seq, a, v, w)
        _seq[a] = rng.choice(20, p=cond_prob)
    return _seq


def gibbs_sampling(init_seq, n_seq, v, w, n_steps, seed=42):
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
    v, w : np array
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
        for _ in range(n_steps):
            seq = gibbs_step(seq, v, w, rng)
        yield seq


def gibbs_sampling_old(init_seq, n_seq, v, w, n_steps, burnin=10):
    """Gibbs sampling process
    #### OLD VERSION - checkout the new version above
    Return a simulated MSA in numerical represenations, according to
    the model v and w using Gibbs sampling process

    Parameters
    ----------
    init_seq : initial sequences
    N : int
        number of desired sequences
    v, w : np array
        the statistical model
    T : int
        number of gibbs steps

    Results
    -------
    np array (LxN)
        Result of simulation
    """
    len_seq = v.shape[0]
    seqs = np.empty(shape=[n_seq, len_seq]).astype('int64')
    seqs[0] = init_seq
    # # burn in phase:
    # for _ in range(burnin):
    #     seqs[0] = gibb_step(seqs[0], v, w)

    for i in range(n_seq-1):
        for _ in range(n_steps):
            seqs[i] = gibbs_step(seqs[i], v, w)
        seqs[i+1] = seqs[i]
    return seqs


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


# Generate parameters
def generate_v(L, type_of_v, exclusion=[], seed=42):
    # Define v
    rng = np.random.default_rng(seed)
    if type_of_v == 'simple':
        v = np.ones([L, 20])
    else:
        choices = np.delete(np.arange(L), exclusion)
        v_highrange = [8,  9]
        v_lowrange = [-0.02,  0.05]
        high_v_residue = rng.choice(choices, size=20, replace=False)
        low_v_residue = [i for i in range(L) if i not in high_v_residue]
        v = np.ones([L, 20])
        for i in high_v_residue:
            v[i] = 0
            for j in np.random.choice(range(0, 20), 3):
                v[i, j] = np.random.uniform(v_highrange[0], v_highrange[1])
        for i in low_v_residue:
            v[i] = v[i]*np.random.uniform(v_lowrange[0], v_lowrange[1])
    return v


def generate_topology(nodes,  topo,  n_edges=None):
    # edge nomenclature: [i, a, j, b]: amino acid a at res i to aa b at j
    edges = []
    r = len(nodes)
    if topo == 'line':
        for i in range(r-1):
            if nodes[i][0] < nodes[i+1][0]:
                edges.append(nodes[i]+nodes[i+1])

    elif topo == 'circle':
        for i in range(r-1):
            edges.append(nodes[i]+nodes[i+1])
        edges.append(nodes[-1]+nodes[0])

    elif topo == 'fully-connected':
        for i in range(r):
            for j in range(i+1, r):
                if nodes[i][0] != nodes[j][0]:
                    edges.append(nodes[i]+nodes[j])

    elif topo == 'randomly-connected':
        # 50% of fully connected edges are removed
        if n_edges is None:
            print('number of edges is required for random connection')
            return None
        for i in range(r):
            for j in range(i+1, r):
                edges.append(nodes[i]+nodes[j])
        edges = [edges[retain] for retain in
                 [i for i in np.random.default_rng().choice(range(len(edges)),
                  size=n_edges, replace=False)]]
    return edges


def generate_w(L, edges):
    # edge nomenclature: [i, a, j, b]: amino acid a at res i to aa b at j
    w = np.zeros([L, L, 20, 20])
    for edge in edges:
        if edge[0] == edge[2]:
            continue
        w[edge[0], edge[2], edge[1], edge[3]] = 5  # a@i to b@j
        w[edge[2], edge[0], edge[3], edge[1]] = 5  # b@j to a@b
    return w


def generate_v_and_w(L, v_type, edges):
    # Create custom binary matrices
    exclusion = []
    for edge in edges:
        exclusion += [edge[0], edge[2]]
    x_single = generate_v(L, v_type, exclusion)
    x_pair = generate_w(L, edges)
    return x_single, x_pair
