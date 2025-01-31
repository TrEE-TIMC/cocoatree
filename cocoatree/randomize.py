"""Module to perform randomization of alignments"""


import numpy as np
from .__params import lett2num, __pseudo_count_ref
from .statistics.position import _compute_aa_freqs, \
    _compute_background_freqs
from .statistics.pairwise import _compute_aa_joint_freqs, compute_sca_matrix
from .deconvolution import eigen_decomp
from sklearn.utils import check_random_state


def _random_aln(fia, n_seq, random_state):
    """
    Generate a random alignment with n_seq sequences and amino acid frequency
    at each position fia

    Arguments
    ---------
    fia : frequency of amino acid *a* at position *i*

    n_seq : number of sequences

    random_state : np.random.RandomState

    Returns
    -------
    msa_str : Random alignment as a list of string sequences

    binarray : Random alignment as a binary array
    """

    n_pos = fia.shape[0]
    msa_rand = np.zeros((n_seq, n_pos), dtype=int)
    for i in range(n_pos):
        Maa = random_state.multinomial(n_seq, fia[i, :])
        col = np.array([], dtype=int)
        for aa, M in enumerate(Maa):
            col = np.append(col, np.tile(aa, M))
        random_state.shuffle(col)
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


def randomization(sequences, n_rep, seq_weights=1,
                  pseudo_count=__pseudo_count_ref,
                  kmax=6,
                  random_state=None):
    """
    Randomize the alignment while preserving the frequencies of amino acids at
    each position and compute the resulting spectrum of coevolution matrix.

    Arguments
    ---------
    sequences : multiple sequence alignment

    n_rep : int
        number of iterations of randomization

    weights : {1, np.array (n_seq)}, default : 1
        weights provided to :fun:`cocoatree.statistics.pairwise.aa_joint_freq`

            - If int, assumes equal weights on all sequences
            - If, vector of sequence weights for each sequence

    lambda_coef : float, optional, default: 0.03
        pseudo-counts

    kmax : int, optional, default: 6
        number of eigenvectors to keep for each randomized iteration

    random_state : int or RandomState instance, default=None
        Determines random number generation for centroid initialization. Pass
        an int for reproducible output across multiple function calls.

    Returns
    -------
    vect_rand :

    val_rand :
    """
    random_state = check_random_state(random_state)

    n_seq, n_pos = len(sequences), len(sequences[0])

    # Create a vector of sequence weights = 1 if equal weighting
    if isinstance(seq_weights, int) and seq_weights == 1:
        seq_weights = np.ones(n_seq)

    aa_freqs = _compute_aa_freqs(sequences,
                                 seq_weights=seq_weights,
                                 pseudo_count=pseudo_count)
    bkgd_freqs = _compute_background_freqs(aa_freqs,
                                           sequences,
                                           seq_weights=seq_weights,
                                           pseudo_count=pseudo_count)

    # initialize for eigenvectors
    vect_rand = np.zeros((n_rep, n_pos, kmax))
    # initialize for eigenvalues
    val_rand = np.zeros((n_rep, n_pos))
    for rep in range(n_rep):
        rand_sequences = _random_aln(aa_freqs, n_seq,
                                     random_state=random_state)[0]
        aa_joint_freqs = _compute_aa_joint_freqs(rand_sequences, seq_weights,
                                                 pseudo_count=pseudo_count)
        # Compute coevolution matrix for the randomized alignment
        sca_rand = compute_sca_matrix(aa_joint_freqs, aa_freqs,
                                      bkgd_freqs)[1]

        eig_val, eig_vec = eigen_decomp(sca_rand)
        vect_rand[rep, :, :] = eig_vec[:, :kmax]
        val_rand[rep, :] = eig_val

    return vect_rand, val_rand
