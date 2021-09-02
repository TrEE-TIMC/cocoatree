"""
=============================================
Running SCA on the S1A serine protease family
=============================================

This example showcases how to run the whole SCA pipeline on the S1A protein
family. This reproduces the results of Halabi et al, 2009.
"""

import itertools
import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import quantile_transform
from akasthesia.coevolution import Alignment, CoevolutionMatrix

###############################################################################
# Loading the S1A serine data
###############################################################################
#
# This example is based on the serine protein family, studied by Halabi et al,
# and then further by Rivoire at al.
alg = Alignment.from_file('data/s1Ahalabi_1470.an')

n_sequences = alg.seq_count()
n_residues = alg.filtered_seq_len()
print("Loaded %d sequence of with %d residues" % (n_sequences, alg.seq_len()))

###############################################################################
# Calculate Weights - Sequence Weight and Position Weights
###############################################################################
#
# The first step of the pipeline consists in evaluating two sets of weights.
# First, what we denote by sequence weight or organism weight: these weights
# consists in correcting for the "phylogeny effect," ie the effect of an
# uneven sampling of sequences due to the overreprentation of some organisms
# in datasets (such as E. coli). The second set of weights corresponds to what
# we denote by the "residue weights" or "position weights." These weights aim
# at correcting for the fact that some residues are much more conserved than
# others.

seq_weights = alg.sca_seq_weights()
pos_weights = alg.position_weights_sca(seq_weights)

###############################################################################
#  Estimating the coevolution matrix
###############################################################################
#
# Then, we estimate the co-evolution matrix from the alignment, corrected with
# these sequence and residue weights. The goal here is to estimate robustly
# whether pairs of residues co-evolve.

weighted_cij = alg.sca_coevolution(seq_weights=seq_weights,
                                   pos_weights=pos_weights)

figure, ax = plt.subplots()
im = ax.matshow(weighted_cij.matrix)
ax.set_yticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
ax.set_xticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
cbar = ax.figure.colorbar(im, ax=ax)

###############################################################################
# Statistical analysis using Independant Component Analysis (ICA)
###############################################################################
#
# The last step of the analysis consists in performing somewhat of an ICA on
# the coevolution matrix estimated in the previous steps.

n_components = 6
vica, w = CoevolutionMatrix.ICA(weighted_cij.matrix, kmax=n_components)

# Set some colors based on the projects on the second component
colors = quantile_transform(
    vica[:, 1].reshape(-1, 1), n_quantiles=n_residues).flatten()
fig, axes = plt.subplots(n_components-1, n_components-1, tight_layout=True)
for i, j in itertools.product(range(n_components), range(n_components)):
    if i >= j:
        if j >= 1 and (i < n_components-1):
            fig.delaxes(axes[i, j-1])
        continue
    ax = axes[i, j-1]
    ax.scatter(vica[:, j], vica[:, i], marker=".", c=colors)
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines["top"].set_linewidth(0)
    ax.spines["right"].set_linewidth(0)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.tick_params(axis='both', which='major', labelsize="xx-small")
    ax.tick_params(axis='both', which='minor', labelsize="xx-small")

    ax.text(
        xmin - (xmax - xmin) * 0.05, 0, "IC%d" % (i+1),
        fontweight="bold", rotation=90, horizontalalignment="right",
        fontsize="x-small", verticalalignment="center")
    ax.text(
        0, ymin - (ymax - ymin) * 0.05,
        "IC%d" % (j+1), fontsize="x-small", fontweight="bold",
        horizontalalignment="center", verticalalignment="top")

###############################################################################
# We can see from the results of the ICA that some components clearly captures
# clusters of residues. These likely correspond to sectors.
