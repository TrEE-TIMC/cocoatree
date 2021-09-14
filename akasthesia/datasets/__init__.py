from akasthesia.coevolution import Alignment
import numpy as np

from . import simulation

# import os


# def s1a(threshold=.2):
#     """S1A Serine Protease protein family

#     Parameters
#     ------
#     threshold : float
#         Default .2

#     Returns
#     ------
#     Alignment :
#         Alignment object of S1A protein family
#     """
#     return Alignment.from_file(os.getcwd()+'/data/s1Ahalabi_1470.an', threshold=threshold)


def simple_simulated():
    """Generate a simple simulated dataset

    Returns:
    --------
    Alignment :
        Alignment object of the generated dataset
    """
    N = 100
    L = 75
    init_seqs = np.random.randint(20, size=[N, L])
    v, w = simulation.generate_v_and_w(L, 'simple', [])
    seqs = simulation.gibbs_sampling(init_seqs, v, w, 100)
    return simulation.to_Alignment(simulation.num_to_aa(seqs, N), N)


def simulated(N, v, w, T, seed=42, burnin=0):
    """Generate simulated dataset using Gibbs sampling (Volberg et al, 2018)

    Take the model in the form of v (Lx20) and w (LxLx20x20)
    and generate a MSA data set. The probability of encountering
    a sequence is

    .. math:
        p(x|v,w) ~ exp(sum{i}(v_i(x_i)) + sum{ij}(w_ij(x_i, x_j)))


    Parameters
    -------
    N : int
        Number of sequences in the output MSA
    v, w : np array
        the model as descibed
    t : int
        the number of Gibbs steps

    Returns
    -------
    Alignment :
        Alignment object of the generated MSA

    """
    L = v.shape[0]
    init_seqs = np.random.randint(20, size=L)
    seqs = simulation.gibbs_sampling(init_seqs, N+burnin, v, w, T, seed=seed)
    for i in range(burnin):
        seqs.__next__()
    return simulation.to_Alignment(seqs)
