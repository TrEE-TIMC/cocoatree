import numpy as np
from numpy import testing as npt

from datasets import simulation

v = np.zeros([75, 20])
w = np.zeros([75, 75, 20, 20])


def test_cond_prob_simple():
    cond_prob = simulation.conditional_prob(np.random.randint(20, size=75), 0, v, w)
    npt.assert_almost_equal(
        cond_prob, np.tile(0.05, (20))
        )
    npt.assert_almost_equal(
        np.sum(cond_prob), 1
    )


def test_gibbs_sampling():
    msa = np.array([_ for _ in simulation.gibbs_sampling(np.random.randint(20, size=75),
                                                         1000, v, w, 1)])
    npt.assert_equal(msa.shape, (1000, 75))
