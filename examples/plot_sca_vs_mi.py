"""
=============================
Mutual information versus SCA
=============================

In this example, we are comparing the results of the co-evolution analysis on
serine proteases using SCA and the mutual information.
"""

from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.msa import filter_sequences

from cocoatree.statistics import compute_all_frequencies
from cocoatree.statistics.pairwise import compute_sca_matrix
from cocoatree.statistics.pairwise import compute_mutual_information_matrix

import matplotlib.pyplot as plt


# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset.

serine_dataset = load_S1A_serine_proteases()
seq_id = serine_dataset["sequence_ids"]
sequences = serine_dataset["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)

# %%
# Filtering of the multiple sequence alignment
# --------------------------------------------
#
# We are going to filter and clean the MSA
seq_id_kept, seq_kept, pos_kept = filter_sequences(seq_id, sequences)

# %%
# Compute the SCA matrix
# ----------------------

aa_freq, background_frequencies, joint_frequencies = compute_all_frequencies(
    seq_kept, seq_weights=None)
Cijab_raw, sca = compute_sca_matrix(joint_freqs=joint_frequencies,
                                    aa_freq=aa_freq,
                                    background_freq=background_frequencies)

# %%
# Compute the Mutual information matrix
# -------------------------------------
normalized_mi = compute_mutual_information_matrix(seq_kept)
mi = compute_mutual_information_matrix(seq_kept, normalize=False)

# %%
# Compare MI versus SCA


fig, axes = plt.subplots(ncols=3, figsize=(8, 4), tight_layout=True)
ax = axes[0]
ax.matshow(sca)
ax.set_title("SCA matrix")

ax = axes[1]
ax.matshow(normalized_mi)
ax.set_title("Normalized Mutual Information")

ax = axes[2]
ax.matshow(mi)
ax.set_title("Mutual Information")
