"""
===========================================================
Perform full SCA analysis on the Rhomboid proteases dataset
===========================================================

In this example, we aim to reproduce the coevolution analysis performed
in Mihaljevic and Urban, 2020 (doi:
https://doi-org.insb.bib.cnrs.fr/10.1016/j.str.2020.07.015).

Using the SCA metric, they identified two protein sectors in the rhomboid
protein superfamily and successfully performed experimental sector
grafting which converted rhomboid architectural and proteolytic properties.

Let's see if we can identify those sectors using COCOA-Tree!
"""

# Author: Margaux Jullien <margaux.jullien@univ-grenoble-alpes.fr>
#         Nelle Varoquaux <nelle.varoquaux@univ-grenoble-alpes.fr>
#         Ivan Junier <ivan.junier@univ-grenoble-alpes.fr>
# License: TBD

# %%
# Import necessary
import cocoatree.datasets as c_data
# import cocoatree.io as c_io
import cocoatree.msa as c_msa
import cocoatree.statistics.position as c_pos
import cocoatree.statistics.pairwise as c_pw
import cocoatree.deconvolution as c_deconv

import matplotlib.pyplot as plt
import numpy as np

# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset. In this case, we can directly load the S1
# serine protease dataset provided in :mod:`cocoatree`. To work on your own
# dataset, you can use the `cocoatree.io.load_msa` function.

rhomboid_dataset = c_data.load_rhomboid_proteases()
loaded_seqs = rhomboid_dataset["alignment"]
loaded_seqs_id = rhomboid_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)

print(f"The loaded MSA has {n_loaded_seqs} sequences and {n_loaded_pos} \
      positions.")

# %%
# MSA filtering
# -------------
#
# We clean the loaded MSA by filtering both sequences and positions.
# As the authors provide the already processed alignment, this step doesn't
# modify the alignment.

sequences, sequences_id, positions = c_msa.filter_sequences(
    loaded_seqs, loaded_seqs_id, gap_threshold=0.4, seq_threshold=0.2)
n_pos = len(positions)
print(f"After filtering, we have {n_pos} remaining positions.")
print(f"After filtering, we have {len(sequences)} remaining sequences.")

# %%
# Compute the matrix of pairwise sequence identity
# ------------------------------------------------

identity_matrix = c_msa.compute_seq_identity(sequences)

fig, ax = plt.subplots()
m = ax.imshow(identity_matrix, vmin=0, vmax=1, cmap='inferno')
ax.set_xlabel("sequences", fontsize=10)
ax.set_ylabel("sequences", fontsize=10)
ax.set_title('Matrix of pairwise sequence identity', fontweight="bold")
cb = fig.colorbar(m)
cb.set_label("Pairwise sequence identity", fontweight="bold")

# %%
# Compute sequence weights
seq_weights, m_eff = c_pos.compute_seq_weights(sequences)
print('Number of effective sequences %d' %
      np.round(m_eff))

# %%
# Compute conservation along the MSA
# ----------------------------------
#
# We compute the conservation of residues over each position of the alignment.
Di = c_pos.compute_conservation(sequences, seq_weights)

fig, ax = plt.subplots(figsize=(9, 4))
xvals = [i+1 for i in range(len(Di))]
xticks = [0, 50, 100, 150]
ax.bar(xvals, Di, color='k')
plt.tick_params(labelsize=11)
ax.set_xticks(xticks)
ax.set_xlabel('Residue positions', fontsize=14)
ax.set_ylabel('Di', fontsize=14)

# %%
# Compute the SCA coevolution matrix
# ----------------------------------

SCA_matrix = c_pw.compute_sca_matrix(sequences, seq_weights)

fig, ax = plt.subplots()
im = ax.imshow(SCA_matrix, vmin=0, vmax=1.4, cmap='inferno')
ax.set_xlabel('Residue positions', fontsize=10)
ax.set_ylabel(None)
ax.set_title('SCA matrix', fontweight='bold')
fig.colorbar(im, shrink=0.7)

# %%
# Extraction of principal components (PCA analysis)
# and of independent components (ICA analysis)
# (this can take some time because of randomization)
# -------------------------------------------------

n_components = 9
principal_components = c_deconv.extract_principal_components(
    SCA_matrix)
idpt_components = c_deconv.extract_independent_components(
    SCA_matrix,
    n_components=n_components)


# %%
# Plot components
# ---------------

n_components_to_plot = n_components
if n_components_to_plot % 2 != 0:
    print('Odd number of components: the last one is discarded for \
          visualization')
    n_components_to_plot -= 1
print(n_components_to_plot)

pairs = [[x, x+1] for x in range(0, n_components_to_plot, 2)]
ncols = len(pairs)
plt.rcParams['figure.figsize'] = 14, 8
fig, axes = plt.subplots(nrows=2, ncols=len(pairs), tight_layout=True)
for k, [k1, k2] in enumerate(pairs):
    ax = axes[0, k]
    ax.plot(principal_components[k1], principal_components[k2], 'ok')
    ax.set_xlabel("PC %i" % (k1+1), fontsize=16)
    ax.set_ylabel("PC %i" % (k2+1), fontsize=16)

    ax = axes[1, k]
    ax.plot(idpt_components[k1], idpt_components[k2], 'ok')
    ax.set_xlabel("IC %i" % (k1+1), fontsize=16)
    ax.set_ylabel("IC %i" % (k2+1), fontsize=16)

# %%
# Extract sectors
# ---------------

sectors = c_deconv.extract_sectors(idpt_components, SCA_matrix)

print('Sector positions on (filtered) sequences:')
for isec, sec in enumerate(sectors):
    print('sector %d: %s' % (isec+1, sec))

# %%
# Plot coevolution within and between the sectors
# Each white square corresponds to a sector, with the residues ordered in
# decreasing contribution to the independent component associated from top to
# bottom and from left to right.
fig, ax = plt.subplots(tight_layout=True)

sector_sizes = [len(sec) for sec in sectors]
cumul_sizes = sum(sector_sizes)
sorted_pos = [s for sec in sectors for s in sec]

im = ax.imshow(SCA_matrix[np.ix_(sorted_pos, sorted_pos)],
               vmin=0, vmax=2,
               interpolation='none', aspect='equal',
               extent=[0, cumul_sizes, 0, cumul_sizes], cmap='inferno')
ax.set_title('SCA matrix, sorted according to sectors')
cb = fig.colorbar(im)
cb.set_label("coevolution level")

line_index = 0
for i in range(n_components):
    ax.plot([line_index + sector_sizes[i], line_index + sector_sizes[i]],
            [0, cumul_sizes], 'w', linewidth=2)
    ax.plot([0, cumul_sizes], [cumul_sizes - line_index,
                               cumul_sizes - line_index], 'w', linewidth=2)
    line_index += sector_sizes[i]

# %%
# Do the same but for the SCA matrix where the global correlation mode has
# been removed and, hence, such that sectors are better highlighted.
# See Rivoire et al., PLOSCB, 2016

# removing global model (ngm = no global mode),
# i.e., removing first principal component
SCA_matrix_ngm = c_deconv.substract_first_principal_component(SCA_matrix)

# plotting
fig, ax = plt.subplots(tight_layout=True)
im = ax.imshow(SCA_matrix_ngm[np.ix_(sorted_pos, sorted_pos)], vmin=0, vmax=1,
               interpolation='none', aspect='equal',
               extent=[0, sum(sector_sizes), 0, sum(sector_sizes)],
               cmap='inferno')
ax.set_title('SCA matrix without global mode, sorted according to sectors')
cb = fig.colorbar(im)
cb.set_label("coevolution level")

line_index = 0
for i in range(n_components):
    ax.plot([line_index + sector_sizes[i], line_index + sector_sizes[i]],
            [0, sum(sector_sizes)], 'w', linewidth=2)
    ax.plot([0, sum(sector_sizes)], [sum(sector_sizes) - line_index,
                                     sum(sector_sizes) - line_index],
            'w', linewidth=2)
    line_index += sector_sizes[i]
