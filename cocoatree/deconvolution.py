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


def icList(Vica, n_component, Cij, p_cut=0.95):
    """
    Produces a list of positions contributing to each independent component
    (IC) above a defined statistical cutoff (p_cut, the cutoff on the CDF of
    the t-distribution fit to the histogram of each IC). Any position above the
    cutoff on more than one IC are assigned to one IC based on which group of
    positions to which it shows a higher degree of coevolution. Additionally
    returns the numeric value of the cutoff for each IC, and the pdf fit, which
    can be used for plotting/evaluation.
    Based on Rivoire et al. (2016): \
        https://doi.org/10.1371/journal.pcbi.1004817

    Arguments
    ---------
    Vica : ndarray,
        independent components

    n_component : int,
        number of independent components chosen

    Cij : numpy.ndarray,
        coevolution matrix

    p_cut : int,
        cutoff on the CDF of the t-distribution fit to the histogran of each IC

    Returns
    -------
    selected_res : list of cocoatree.deconvolution.Unit,
        positions of the selected residues for each independent component.
        Beware that if the alignment used for the analysis has been filtered,
        those are the positions on the filtered alignment and not on the
        original alignment, a mapping of the positions may be needed.

    ic_size : list,
        number of selected residues for each component.

    sorted_pos : list,
        positions of the residues sorted by decreasing contribution for each
        component.

    cutoff : list,
        numeric value of the cutoff for each component.

    scaled_pdf : list of np.ndarrays,
        scaled probability distribution function for each component.

    all_fits : list,
        t-distribution fits for each component.

    **Example**::
        selected_res, ic_size, sorted_pos, cutoff, scaled_pdf, all_fits = \
            icList(Vica, n_component, Cij, p_cut=0.95)
    """

    # do the PDF/CDF fit, and assign cutoffs
    Npos = len(Vica)
    cutoff = list()
    scaled_pdf = list()
    all_fits = list()
    for k in range(n_component):
        pd = t.fit(Vica[:, k])
        all_fits.append(pd)
        iqr = scoreatpercentile(Vica[:, k], 75) - scoreatpercentile(
            Vica[:, k], 25
        )
        binwidth = 2 * iqr * (len(Vica[:, k]) ** (-0.33))
        nbins = round((max(Vica[:, k]) - min(Vica[:, k])) / binwidth)
        h_params = np.histogram(Vica[:, k], int(nbins))
        x_dist = np.linspace(min(h_params[1]), max(h_params[1]), num=100)
        area_hist = Npos * (h_params[1][2] - h_params[1][1])
        scaled_pdf.append(area_hist * (t.pdf(x_dist, pd[0], pd[1], pd[2])))
        cd = t.cdf(x_dist, pd[0], pd[1], pd[2])
        tmp = scaled_pdf[k].argmax()
        if abs(max(Vica[:, k])) > abs(min(Vica[:, k])):
            tail = cd[tmp: len(cd)]
        else:
            cd = 1 - cd
            tail = cd[0:tmp]
        diff = abs(tail - p_cut)
        x_pos = diff.argmin()
        cutoff.append(x_dist[x_pos + tmp])

    # select the positions with significant contributions to each IC
    ic_init = list()
    for k in range(n_component):
        ic_init.append([i for i in range(Npos) if Vica[i, k] > cutoff[k]])

    # construct the sorted, non-redundant iclist
    sorted_pos = list()
    ic_size = list()
    selected_res = list()
    icpos_tmp = list()
    Cij_nodiag = Cij.copy()
    for i in range(Npos):
        Cij_nodiag[i, i] = 0
    for k in range(n_component):
        icpos_tmp = list(ic_init[k])
        for kprime in [kp for kp in range(n_component) if kp != k]:
            tmp = [v for v in icpos_tmp if v in ic_init[kprime]]
            for i in tmp:
                remsec = np.linalg.norm(
                    Cij_nodiag[i, ic_init[k]]
                ) < np.linalg.norm(Cij_nodiag[i, ic_init[kprime]])
                if remsec:
                    icpos_tmp.remove(i)
        sorted_pos += sorted(icpos_tmp, key=lambda i: -Vica[i, k])
        ic_size.append(len(icpos_tmp))
        s = Unit()
        s.items = sorted(icpos_tmp, key=lambda i: -Vica[i, k])
        s.col = k / n_component
        s.vect = -Vica[s.items, k]
        selected_res.append(s)
    return selected_res, ic_size, sorted_pos, cutoff, scaled_pdf, all_fits


def choose_num_components(eigenvalues, rand_eigenvalues):
    """
    Given the eigenvalues of the coevolution matrix (eigenvalues), and the
    eigenvalues for the set of randomized matrices (rand_eigenvalues), return
    the number of significant eigenmodes as those above the average first
    eigenvalue plus 3 standard deviations.
    Based on S1 text of Rivoire et al. (2016)

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
                              (rand_eigenvalues[:, -2].mean() +
                               (3 * rand_eigenvalues[:, -2].std()))].shape[0]

    return n_component
