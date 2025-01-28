"""
============================================================
Perform full SCA analysis on the S1A serine protease dataset
============================================================

This example shows the full process to perform a complete coevolution
analysis in order to detect protein sectors from data importation, MSA
filtering, computation of positional and joint amino acid frequencies,
and computation of the SCA coevolution matrix.

The matrix is then decomposed into principal components and independent
component analysis is performed.

In the end, we export a fasta file of the residues contributing to the first
independent component.

"""

# Author: Margaux Jullien <margaux.jullien@univ-grenoble-alpes.fr>
#         Nelle Varoquaux <nelle.varoquaux@univ-grenoble-alpes.fr>
#         Ivan Junier <ivan.junier@univ-grenoble-alpes.fr>
# License: TBD

# %%
# Import necessary
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.io import export_fasta, load_pdb, export_sector_for_pymol
from cocoatree.msa import filter_gap_seq, filter_gap_pos, seq_weights, \
    map_to_pdb
from cocoatree.statistics.position import aa_freq_at_pos, \
    compute_background_frequencies, compute_rel_entropy
from cocoatree.statistics.pairwise import aa_joint_freq, compute_sca_matrix, \
    compute_seq_identity
from cocoatree.deconvolution import eigen_decomp, compute_ica, \
    choose_num_components, extract_positions_from_IC
from cocoatree.randomize import randomization
import matplotlib.pyplot as plt
import numpy as np

# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset. In this case, we can directly load the S1
# serine protease dataset provided in :mod:`cocoatree`. To work on your on
# dataset, you can use the :fun:`cocoatree.io.load_msa` function.

serine_dataset = load_S1A_serine_proteases()
seq_id = serine_dataset["sequence_ids"]
sequences = serine_dataset["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)

print(f"The loaded MSA has {n_seq} sequences and {n_pos} positions.")
# %%
# MSA filtering
# -------------
#
# We are going to clean a bit the loaded MSA by filtering both sequences and
# positions.
#
# Filter overly gapped positions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

filt_seqs, pos_kept = filter_gap_pos(sequences, threshold=0.4)
n_pos_kept = len(pos_kept)
print(f"After filtering, we have {n_pos_kept} remaining positions.")

# %%
# Filter overly gapped sequences
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

seq_id_kept, seq_kept = filter_gap_seq(seq_id, filt_seqs, threshold=0.2,
                                       filtrefseq=False)
print(f"After filtering, we have {len(seq_kept)} remaining sequences.")

# %%
# Compute the matrix of pairwise sequence identity
# ------------------------------------------------

sim_matrix = compute_seq_identity(seq_kept)

fig, ax = plt.subplots()
m = ax.imshow(sim_matrix, vmin=0, vmax=1, cmap='inferno')
ax.set_xlabel("sequences", fontsize=10)
ax.set_ylabel("sequences", fontsize=10)
ax.set_title('Matrix of pairwise sequence identity', fontweight="bold")
cb = fig.colorbar(m)
cb.set_label("Pairwise sequence identity", fontweight="bold")

# %%
# Compute sequence weights
weights, n_eff_seq = seq_weights(sim_matrix)
print(f"Number of effective sequences {n_eff_seq}")

# %%
# compute allele frequencies
aa_freq = aa_freq_at_pos(seq_kept, lambda_coef=0.03, weights=weights)

# %%
# Compute background frequencies
background_frequencies = compute_background_frequencies(aa_freq)

# %%
# Compute joint allele frequencies
fijab, fijab_ind = aa_joint_freq(seq_kept, weights=weights, lambda_coef=0.03)

# %%
# Compute conservation along the MSA
# ----------------------------------
Dia, Di = compute_rel_entropy(aa_freq, background_frequencies)

fig, ax = plt.subplots(figsize=(9, 4))
xvals = [i+1 for i in range(len(Di))]
xticks = [0, 50, 100, 150, 200, 250]
ax.bar(xvals, Di, color='k')
plt.tick_params(labelsize=11)
ax.set_xticks(xticks)
ax.set_xlabel('Residue position', fontsize=14)
ax.set_ylabel('Di', fontsize=14)

# %%
# Compute the SCA coevolution matrix
# ----------------------------------

Cijab_raw, Cij = compute_sca_matrix(joint_freqs=fijab,
                                    joint_freqs_ind=fijab_ind,
                                    aa_freq=aa_freq,
                                    background_freq=background_frequencies)

fig, ax = plt.subplots()
im = ax.imshow(Cij, vmin=0, vmax=1.4, cmap='inferno')

ax.set_xlabel('Residue', fontsize=10)
ax.set_ylabel(None)
ax.set_title('Coevolution matrix')
fig.colorbar(im, shrink=0.7)

# %%
# Decomposition of the matrix into principal components
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

eigenvalues, eigenvectors = eigen_decomp(Cij)

# %%
# Plot distribution of eigenvalues
fig, ax = plt.subplots()
ax.hist(eigenvalues, bins=100, color="black")
ax.set_ylabel('Number', fontweight="bold")
ax.set_xlabel('Eigenvalue', fontweight="bold")

# %%
# Select number of significant components
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# We use a randomization strategy in order to select the number of components.
# The function :fun:cocoatree.randomize.randomization runs the full SCA
# pipeline on randomized version of the MSA. Make sure that weights and lambda
# coefficient are set the same way as the when performing the analysis on the
# real dataset.
v_rand, l_rand = randomization(seq_kept, n_rep=10,
                               weights=weights, lambda_coef=0.03, kmax=10)
n_components = choose_num_components(eigenvalues, l_rand)
print('n_components = ' + str(n_components))

hist0, bins = np.histogram(l_rand.flatten(), bins=n_pos_kept,
                           range=(0, eigenvalues.max()))
hist1, bins = np.histogram(eigenvalues, bins=n_pos_kept,
                           range=(0, eigenvalues.max()))

fig, ax = plt.subplots()
ax.bar(bins[:-1], hist1, np.diff(bins), color='k')
ax.plot(bins[:-1], hist0/10, 'r', linewidth=3)
ax.set_xlabel('Eigenvalues', fontweight="bold")
ax.set_ylabel('Numbers', fontweight="bold")

# %%
# Independent component analysis (ICA)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

independent_components, W = compute_ica(
    eigenvectors, kmax=n_components, learnrate=0.1,
    iterations=100000)

# Plot results
if n_components % 2 != 0:
    print('Uneven number of axes, discard the last one for visual \
          representation')
    n_components -= 1

pairs = [[x, x+1] for x in range(0, n_components, 2)]
ncols = len(pairs)
plt.rcParams['figure.figsize'] = 14, 8
fig, axes = plt.subplots(nrows=2, ncols=len(pairs), tight_layout=True)
for k, [k1, k2] in enumerate(pairs):
    ax = axes[0, k]
    ax.plot(eigenvectors[:, k1], eigenvectors[:, k2], 'ok')
    ax.set_xlabel("eigenvector %i" % (k1+1), fontsize=16)
    ax.set_ylabel("eigenvector %i" % (k2+1), fontsize=16)

    ax = axes[1, k]
    ax.plot(independent_components[:, k1], independent_components[:, k2], 'ok')
    ax.set_xlabel("independent component %i" % (k1+1), fontsize=16)
    ax.set_ylabel("independent component %i" % (k2+1), fontsize=16)

# %%
# Select residues that significantly contribute to each independent component
ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits = \
    extract_positions_from_IC(independent_components, n_components, Cij,
                              p_cut=0.95)

print(f"Sizes of the {n_components} ICs: {icsize}")

# %%
# Plot coevolution within and between the sectors
# Each white square corresponds to a sector, with the residues ordered in
# decreasing contribution to the independent component associated from top to
# bottom and from left to right.
fig, ax = plt.subplots(tight_layout=True)
im = ax.imshow(Cij[np.ix_(sortedpos, sortedpos)], vmin=0, vmax=2,
               interpolation='none', aspect='equal',
               extent=[0, sum(icsize), 0, sum(icsize)], cmap='inferno')
cb = fig.colorbar(im)
cb.set_label("Coevolution measure")

line_index = 0
for i in range(n_components):
    ax.plot([line_index + icsize[i], line_index + icsize[i]],
            [0, sum(icsize)], 'w', linewidth=2)
    ax.plot([0, sum(icsize)], [sum(icsize) - line_index,
                               sum(icsize) - line_index], 'w', linewidth=2)
    line_index += icsize[i]

# %%
# Export fasta files of the sectors for all the sequences
# Those fasta can then be used for visualization along a phylogenetic tree
# as implemented in the cocoatree.visualization module

sector_1_pos = list(pos_kept[ics[0].items])
sector_1 = []
for sequence in range(len(seq_id)):
    seq = ''
    for pos in sector_1_pos:
        seq += sequences[sequence][pos]
    sector_1.append(seq)

export_fasta(sector_1, seq_id, 'sector_1.fasta')

# %
# Export files necessary for Pymol visualization
# Load PDB file of rat's trypsin
pdb_seq, pdb_pos = load_pdb('data/3TGI.pdb', pdb_id='TRIPSIN', chain='E')
# Map PDB positions on the MSA sequence corresponding to rat's trypsin:
# seq_id='14719441'
pdb_mapping = map_to_pdb(pdb_seq, pdb_pos, sequences, seq_id,
                         ref_seq_id='14719441')
# Export lists of the first sector positions and each residue's contribution
# to the independent component to use for visualization on Pymol.
# The residues are ordered in the list by decreasing contribution score (the
# first residue in the list is the highest contributing)
export_sector_for_pymol(pdb_mapping, independent_components, axis=0,
                        sector_pos=sector_1_pos,
                        ics=ics,
                        outpath='color_sector_1_pymol.npy')
