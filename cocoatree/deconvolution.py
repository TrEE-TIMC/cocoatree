import scipy.linalg as sp
import numpy as np


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
    eigenvectors = eigenvectors[:, np.arange(eigenvalues - 1, -1, -1)]

    return eigenvalues, eigenvectors
