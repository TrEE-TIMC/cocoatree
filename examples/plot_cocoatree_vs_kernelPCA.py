"""
==============================
SVD versus sklearn's kernelPCA
==============================


Small example showing that cocoatree's extract principal components should be
done on a centered co-evolution matrix
"""

# %%
# Import necessary
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.msa import filter_sequences

from cocoatree.statistics.position import compute_conservation
from cocoatree.statistics.pairwise import compute_sca_matrix
from cocoatree.statistics.pairwise import compute_mutual_information_matrix
from cocoatree.statistics.pairwise import compute_apc
import cocoatree.deconvolution as c_deconv
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import KernelCenterer

import matplotlib.pyplot as plt

import numpy as np

# %%
# Load the dataset and compute coevolution matrices
# --------------------------------------------------
#
# We start by importing the dataset.

serine_dataset = load_S1A_serine_proteases()
seq_id = serine_dataset["sequence_ids"]
sequences = serine_dataset["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)

seq_kept, seq_id_kept, pos_kept = filter_sequences(sequences, seq_id)
sca_matrix = compute_sca_matrix(seq_kept)
sca_matrix_centered = KernelCenterer().fit_transform(sca_matrix)

sca_principal_components = c_deconv.extract_principal_components(sca_matrix)
scac_principal_components = c_deconv.extract_principal_components(
    sca_matrix_centered)

pca = KernelPCA(kernel="precomputed", n_components=3)
sca_sklearn = pca.fit_transform(sca_matrix)
scac_sklearn = pca.fit_transform(sca_matrix_centered)

fig, axes = plt.subplots(figsize=(10, 4), ncols=3, tight_layout=True)
ax = axes[0]
ax.scatter(sca_sklearn[:, 0], sca_principal_components[0])
ax.set_xlabel("SCA matrix with sklearn's kernelPCA - PC1")
ax.set_ylabel("SCA matrix with SVD (C1)")


ax = axes[1]
ax.scatter(sca_sklearn[:, 0], scac_sklearn[:, 0])
ax.set_xlabel("SCA centered with sklearn's kernelPCA- PC1")
ax.set_ylabel("SCA without centering with sklearn's kernelPCA - PC1")


ax = axes[2]
ax.scatter(sca_sklearn[:, 0], scac_principal_components[0])
ax.set_xlabel("SCA with sklearn's kernelPCA - PC1")
ax.set_ylabel("SCA centered with SVD - PC1")

