import numpy as np
from numpy import testing as npt

from datasets import simulation


def test_cond_prob_simple():
    
    # Generate a null model
    v = np.zeros([75, 20])
    w = np.zeros([75, 75, 20, 20])
    # Calculate the conditional probability
    cond_prob = simulation.conditional_prob(np.random.randint(20, size=75), 0, v, w)

    # Check if is probability: sum up to 1
    npt.assert_almost_equal(
        np.sum(cond_prob), 1
    )

    # Check if uniform: all close to .05
    npt.assert_almost_equal(
        cond_prob, np.tile(0.05, (20))
        )


def test_cond_prob_v():
    v = np.zeros([75, 20])
    w = np.zeros([75, 75, 20, 20])

    # Pick 3 significant amino acid: the first 2 have 25% probability, the last one have 50%
    choose = np.random.choice(20, size=3)
    v[:, choose[-1]] = 8  # 50%
    v[:, choose[:-1]] = 8 - np.log(4) + np.log(2)  # 25% each

    # Generate conditional probability
    cond_prob = simulation.conditional_prob(np.random.randint(20, size=75), 0, v, w)
    # Check if is probability: sum up to 1
    npt.assert_almost_equal(
        np.sum(cond_prob), 1
    )
    # Check if distribution is correct (up to 2 decimal)
    npt.assert_allclose(
        cond_prob[choose], np.array([.25, .25, .5]),
        .01,
        .02
    )


def test_cond_prob_w():
    pass


def test_gibbs_sampling_simple():
    v = np.zeros([75, 20])
    w = np.zeros([75, 75, 20, 20])
    msa = np.array([_ for _ in simulation.gibbs_sampling(np.random.randint(20, size=75),
                                                         1000, v, w, 1)])
    # check size of the output
    npt.assert_equal(msa.shape, (1000, 75))

    # check entropy: entropy should be close to each other (+-0.01)
    alg = simulation.to_Alignment(msa)
    npt.assert_allclose(alg.vorberg_entropy(), np.random.permutation(alg.vorberg_entropy()),
                        0.01, 0.02)

    # check overall frequencies: all should be relatively close to .05
    freq1 = np.sum([alg.aa_freq_at_pos(alg.filtered_alignment()[:, i], lbda=0)
                    for i in np.arange(75)], axis=0)[1:]/75
    npt.assert_allclose(freq1, np.tile(0.05, 20), .01, .02)
