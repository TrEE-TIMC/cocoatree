"""
=============================
Mutual information versus SCA
=============================

In this example, we are comparing the results of the co-evolution analysis on
serine proteases using SCA and the mutual information.
"""

from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.msa import filter_gap_seq, filter_gap_pos
from cocoatree.statistics.position import aa_freq_at_pos, \
    compute_background_frequencies

from cocoatree.statistics.sequence import compute_seq_weights
from cocoatree.statistics.pairwise import aa_joint_freq, compute_sca_matrix, \
    compute_seq_identity
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
filt_seqs, pos_kept = filter_gap_pos(sequences, threshold=0.4)
seq_id_kept, seq_kept = filter_gap_seq(seq_id, filt_seqs, threshold=0.2,
                                       filtrefseq=False)

# %%
# Compute the SCA matrix
# ----------------------

sim_matrix = compute_seq_identity(seq_kept)
weights, n_eff_seq = compute_seq_weights(sim_matrix)
aa_freq = aa_freq_at_pos(seq_kept, lambda_coef=0.03, weights=weights)
background_frequencies = compute_background_frequencies(aa_freq)
fijab = aa_joint_freq(seq_kept, weights=weights, lambda_coef=0.03)
Cijab_raw, sca = compute_sca_matrix(joint_freqs=fijab,
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
