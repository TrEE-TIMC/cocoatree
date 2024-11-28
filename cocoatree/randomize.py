"""Module to perform randomization of alignments"""


import numpy as np
from __params import lett2num
from statistics.position import aa_freq_at_pos, background_freq
from statistics.pairwise import aa_joint_freq, compute_sca_matrix
from deconvolution import eigen_decomp


def _random_aln(fia, Nseq):
    """
    Generate a random alignment with Nseq sequences and amino acid frequency
    at each position fia

    Arguments
    ---------
    fia : frequency of amino acid *a* at position *i*

    Nseq : number of sequences

    Returns
    -------
    msa_str : Random alignment as a list of string sequences

    binarray : Random alignment as a binary array
    """

    Npos = fia.shape[0]
    msa_rand = np.zeros((Nseq, Npos), dtype=int)
    for i in range(Npos):
        Maa = np.random.multinomial(Nseq, fia[i, :])
        col = np.array([], dtype=int)
        for aa, M in enumerate(Maa):
            col = np.append(col, np.tile(aa, M))
        np.random.shuffle(col)
        msa_rand[:, i] = col

    binarray = np.array(
        [msa_rand == aa for aa in lett2num.values()]).astype(int)

    # create a MSA as a list of string sequences otherwise we have to modify
    # all the functions that compute frequencies (and potentially others)
    msa_str = []
    for seq in range(len(msa_rand)):
        seq_str = str()
        for aa in range(len(msa_rand[seq])):
            seq_str = seq_str + list(lett2num.keys())[msa_rand[seq][aa]]
        msa_str.append(seq_str)

    return msa_str, binarray


def randomization(sequences, Nrep, weights=1, lambda_coef=0.03, kmax=6,
                  metric='SCA', correction=None):
    """
    Randomize the alignment while preserving the frequencies of amino acids at
    each position and compute the resulting spectrum of coevolution matrix.

    Arguments
    ---------
    sequences : multiple sequence alignment

    Nrep : number of iterations of randomization

    Naa : number of amino acids to consider, default = 21

    weights : vector of sequence weights (default = assume equal weights)

    lbda : pseudo-counts frequency

    kmax : number of eigenvectors to keep for each randomized iteration

    metric : coevolution metric to use, 'MI' or 'SCA' (default = 'SCA')

    correction : correction to apply on wthe coevolution matrix, 'APC',
            'entropy' or None (default = None)

    Returns
    -------
    vect_rand :

    val_rand :
    """

    Nseq, Npos = len(sequences), len(sequences[0])

    # Create a vector of sequence weights = 1 if equal weighting
    if isinstance(weights, int) and weights == 1:
        weights = np.ones(Nseq)

    fia = aa_freq_at_pos(sequences, lambda_coef=lambda_coef, weights=weights)
    qa = background_freq(fia, lambda_coef=lambda_coef)

    # initialize for eigenvectors
    vect_rand = np.zeros((Nrep, Npos, kmax))
    # initialize for eigenvalues
    val_rand = np.zeros((Nrep, Npos))
    for rep in range(Nrep):
        msa_random = _random_aln(fia, Nseq)[0]
        fijab, fijab_ind = aa_joint_freq(msa_random, weights,
                                         lambda_coef=lambda_coef)
        # Compute coevolution matrix for the randomized alignment
        if metric == 'SCA':
            Coev_rand = compute_sca_matrix(fijab, fijab_ind, fia, qa)[1]
        # elif metric == 'MI':
        #    Coev_rand = compute_mi_matrix(fijab, fijab_ind)
        # add a correction
        # if correction == 'APC':
        #    Coev_rand = APC(Coev_rand)
        # if correction == 'entropy':
        #    Coev_rand = entropy_correction(Coev_rand)

        eig_val, eig_vec = eigen_decomp(Coev_rand)
        vect_rand[rep, :, :] = eig_vec[:, :kmax]
        val_rand[rep, :] = eig_val

    return vect_rand, val_rand
