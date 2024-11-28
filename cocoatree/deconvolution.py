import scipy.linalg as sp
import numpy as np
import sys


def eigen_decomp(mat):
    """
    Performs eigenvalue decomposition of a square matrix

    Arguments
    ---------
    mat : square coevolution matrix of shape (Npos, Npos)

    Returns
    -------
    eigenvalues : ndarray of shape (Npos)

    eigenvectors : ndarray of shape (Npos, Npos)
    """

    eigenvalues, eigenvectors = sp.eigh(mat)
    eigenvalues = np.sort(eigenvalues)
    eigenvectors = eigenvectors[:, np.arange(eigenvalues.shape[0] - 1, -1, -1)]

    return eigenvalues, eigenvectors


def _basicICA(x, r0, Niter, tolerance=1e-15):
    """
    Basic ICA algorithm, based on work by Bell & Sejnowski (infomax). The input
    data should preferentially be sphered, i.e., x.T.dot(x) = 1
    Source: https://github.com/ranganathanlab/pySCA/

    Arguments
    ---------
    x : LxM input matrix where L = # features and M = # samples

    r : learning rate / relaxation parameter (e.g. r=.0001)

    Niter : number of iterations (e.g. 1000)

    Returns
    -------
    w : unmixing matrix

    change: record of incremental changes during the iterations.

    **Note** r and Niter should be adjusted to achieve convergence, which
    should be assessed by visualizing 'change' with plot(range(iter), change)
    **Example**::
      [w, change] = basicICA(x, r, Niter)
    """

    [L, M] = x.shape
    w = np.eye(L)
    change = list()
    r = r0 / M
    with np.errstate(over="raise"):
        try:
            for _ in range(Niter):
                w_old = np.copy(w)
                u = w.dot(x)
                w += r * (
                    M * np.eye(L) + (1.0 - 2.0 / (1.0 + np.exp(-u))).dot(u.T)
                ).dot(w)
                delta = (w - w_old).ravel()
                val = delta.dot(delta.T)
                change.append(val)
                if np.isclose(val, 0, atol=tolerance):
                    break
                if _ == Niter - 1:
                    print("basicICA failed to converge: " + str(val))
        except FloatingPointError as e:
            sys.exit("Error: basicICA " + str(e))
    return [w, change]


def compute_ica(V, kmax=6, learnrate=0.1, iterations=10000):
    """
    ICA rotation (using _basicICA) with default parameters and normalization of
    outputs.
    Basic ICA algorithm, based on work by Bell & Sejnowski (infomax). The input
    data should preferentially be sphered, i.e., x.T.dot(x) = 1

    Source: https://github.com/ranganathanlab/pySCA/

    Arguments
    ---------
    V : ndarray,
        eigenvectors obtained after matrix decomposition

    kmax : integer,
        number of independent components to retrieve

    learnrate : integer,
        learning rate / relaxation parameter

    iterations : integer,
        number of iterations

    **Note** r and Niter should be adjusted to achieve convergence, which
    should be assessed by visualizing 'change' with plot(range(iter), change)

    Returns
    -------
    Vica : ndarray,
        contributions along each independent components

    W : ndarray of shape (kmax, kmax),
        unmixing matrix

    **Example**::
       Vica, W = rotICA(V, kmax=6, learnrate=.0001, iterations=10000)
    """

    V1 = V[:, :kmax].T
    [W, changes] = _basicICA(V1, learnrate, iterations)
    Vica = (W.dot(V1)).T
    for n in range(kmax):
        imax = abs(Vica[:, n]).argmax()
        Vica[:, n] = (
            np.sign(Vica[imax, n]) * Vica[:, n] / np.linalg.norm(Vica[:, n])
        )
    return Vica, W


def chooseKpos(eigenvalues, rand_eigenvalues):
    """
    Given the eigenvalues of the coevolution matrix (Lsca), and the
    eigenvalues for the set of randomized matrices (Lrand), return
    the number of significant eigenmodes.
    Based on Rivoire et al. (2016)

    Arguments
    ---------
    eigenvalues : ndarray of shape (Npos),
        list of eigenvalues of the coevolution matrix based on the real MSA

    rand_eigenvalues : ndarray of shape (Nrep, Npos), where Nrep is the number
        of randomization iterations
        eigenvalues of all the iteartions of randomized alignments

    Returns
    -------
    kpos : integer,
        number of significant eigenmodes
    """

    kpos = eigenvalues[eigenvalues > (rand_eigenvalues.mean() +
                                      (3 * rand_eigenvalues.std()))].shape[0]

    return kpos
