import scipy.linalg as sp
import numpy as np
import sys
from scipy.stats import t, scoreatpercentile


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


<<<<<<< HEAD
def choose_num_components(eigenvalues, rand_eigenvalues):
=======
class Unit:
    """
    A class for units (sectors, sequence families, etc.)

    Attributes
    ----------

    name :  string describing the unit (ex: 'firmicutes')
    items : set of member items (ex: indices for all firmicutes
            sequence in an alignment)
    col :   color code associated to the unit (for plotting)
    vect :  an additional vector describing the member items (ex: a list
            of sequence weights)
    """

    def __init__(self):
        self.name = ""
        self.items = set()
        self.col = 0
        self.vect = 0


def icList(Vpica, kpos, Csca, p_cut=0.95):
    """
    Produces a list of positions contributing to each independent component
    (IC) above a defined statistical cutoff (p_cut, the cutoff on the CDF of
    the t-distribution fit to the histogram of each IC). Any position above the
    cutoff on more than one IC are assigned to one IC based on which group of
    positions to which it shows a higher degree of coevolution. Additionally
    returns the numeric value of the cutoff for each IC, and the pdf fit, which
    can be used for plotting/evaluation.

    **Example**::

      icList, icsize, sortedpos, cutoff, pd = icList(Vsca, Lsca, Lrand)
    """

    # do the PDF/CDF fit, and assign cutoffs
    Npos = len(Vpica)
    cutoff = list()
    scaled_pdf = list()
    all_fits = list()
    for k in range(kpos):
        pd = t.fit(Vpica[:, k])
        all_fits.append(pd)
        iqr = scoreatpercentile(Vpica[:, k], 75) - scoreatpercentile(
            Vpica[:, k], 25
        )
        binwidth = 2 * iqr * (len(Vpica[:, k]) ** (-0.33))
        nbins = round((max(Vpica[:, k]) - min(Vpica[:, k])) / binwidth)
        h_params = np.histogram(Vpica[:, k], int(nbins))
        x_dist = np.linspace(min(h_params[1]), max(h_params[1]), num=100)
        area_hist = Npos * (h_params[1][2] - h_params[1][1])
        scaled_pdf.append(area_hist * (t.pdf(x_dist, pd[0], pd[1], pd[2])))
        cd = t.cdf(x_dist, pd[0], pd[1], pd[2])
        tmp = scaled_pdf[k].argmax()
        if abs(max(Vpica[:, k])) > abs(min(Vpica[:, k])):
            tail = cd[tmp: len(cd)]
        else:
            cd = 1 - cd
            tail = cd[0:tmp]
        diff = abs(tail - p_cut)
        x_pos = diff.argmin()
        cutoff.append(x_dist[x_pos + tmp])

    # select the positions with significant contributions to each IC
    ic_init = list()
    for k in range(kpos):
        ic_init.append([i for i in range(Npos) if Vpica[i, k] > cutoff[k]])

    # construct the sorted, non-redundant iclist
    sortedpos = list()
    icsize = list()
    ics = list()
    icpos_tmp = list()
    Csca_nodiag = Csca.copy()
    for i in range(Npos):
        Csca_nodiag[i, i] = 0
    for k in range(kpos):
        icpos_tmp = list(ic_init[k])
        for kprime in [kp for kp in range(kpos) if kp != k]:
            tmp = [v for v in icpos_tmp if v in ic_init[kprime]]
            for i in tmp:
                remsec = np.linalg.norm(
                    Csca_nodiag[i, ic_init[k]]
                ) < np.linalg.norm(Csca_nodiag[i, ic_init[kprime]])
                if remsec:
                    icpos_tmp.remove(i)
        sortedpos += sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        icsize.append(len(icpos_tmp))
        s = Unit()
        s.items = sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        s.col = k / kpos
        s.vect = -Vpica[s.items, k]
        ics.append(s)
    return ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits


def chooseKpos(eigenvalues, rand_eigenvalues):
>>>>>>> origin/main
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
    n_component : integer,
        number of significant eigenmodes
    """

    n_component = eigenvalues[eigenvalues >
                              (rand_eigenvalues.mean() +
                               (3 * rand_eigenvalues.std()))].shape[0]

    return n_component
