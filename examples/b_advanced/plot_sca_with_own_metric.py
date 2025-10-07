"""
========================================
Running Cocoatree with your own metric
========================================

It is possible to perform Cocoatree's SCA pipeline by defining it's own
metric. In this example, we show how to use a callable to compute the
coevolution metric.
"""
import numpy as np

import cocoatree.datasets as c_data
import cocoatree
from sklearn.metrics import pairwise
from sklearn.preprocessing import OneHotEncoder

# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset. In this case, we can directly load the S1
# serine protease dataset provided in :mod:`cocoatree`. To work on your on
# dataset, you can use the `cocoatree.io.load_msa` function.

serine_dataset = c_data.load_S1A_serine_proteases()
loaded_seqs = serine_dataset["alignment"]
loaded_seqs_id = serine_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)


def compute_correlation(sequences, seq_weights=None, freq_regul=0.08):
    seqs = np.array([[i for i in s] for s in sequences])
    X = OneHotEncoder().fit_transform(seqs.T)
    return pairwise.pairwise_kernels(X)


coevol_matrix, results = cocoatree.perform_sca(
    loaded_seqs_id, loaded_seqs, n_components=3,
    coevolution_metric=compute_correlation)
