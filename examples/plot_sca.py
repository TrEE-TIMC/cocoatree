"""
=====================================
Running SCA on the S1A protein family
=====================================

This example showcases how to run the whole SCA pipeline on the S1A protein
family. This reproduces the results of Halabi et al, 2009.
"""

from akasthesia.coevolution import Alignment, CoevolutionMatrix
import numpy as np
import matplotlib.pyplot as plt

# Reading from files
alg = Alignment.from_file("data/s1Ahalabi_1470.an")
print("Loaded " + str(alg.seq_count()) + " sequence of size " + str(alg.seq_len()))
print("After trimming:", str(np.shape(alg.filtered_alignment())[0]),
      "sequences of length", str(np.shape(alg.filtered_alignment())[1]))

###############################################################################
# First: Calculate Weights - Sequence Weight and Position Weights

sca_weights = alg.sca_seq_weights()
pos_weights = alg.position_weights_sca(sca_weights)

###############################################################################
#  Coevolution matrix - weighted

weighted_cij = alg.sca_coevolution(seq_weights=sca_weights,
                                   pos_weights=pos_weights)

figure, ax = plt.subplots()
figure.set_size_inches(20, 10)
im = ax.matshow(weighted_cij.matrix)
ax.set_yticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
ax.set_xticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
cbar = ax.figure.colorbar(im, ax=ax)

plt.show()

###############################################################################
#  Perform ICA

vica, w = CoevolutionMatrix.ICA(weighted_cij.matrix)

figure, ax = plt.subplots(2, 3)
figure.set_size_inches(15, 10)
ax[0, 0].scatter(x=vica[:, 0], y=vica[:, 1])
ax[0, 1].scatter(x=vica[:, 2], y=vica[:, 3])
ax[0, 2].scatter(x=vica[:, 4], y=vica[:, 5])
ax[1, 0].scatter(x=vica[:, 1], y=vica[:, 2])
ax[1, 1].scatter(x=vica[:, 3], y=vica[:, 4])
ax[1, 2].scatter(x=vica[:, 0], y=vica[:, 5])

###############################################################################
#  IC Histogram

figure, axes = plt.subplots(2, 3)
figure.set_size_inches(15, 10)

i = 0
for ax1 in axes:
    ax1[0].set_ylabel("Count")
    for ax in ax1:
        ax.hist(vica[:, i], bins=20, rwidth=0.95)
        i += 1
        ax.set_xlabel("IC" + str(i))
