"""
=========================
Perform full SCA analysis
=========================

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
from cocoatree.io import load_MSA, export_fasta
from cocoatree.msa import filter_gap_seq, filter_gap_pos, seq_weights
from cocoatree.statistics.position import aa_freq_at_pos, background_freq
from cocoatree.statistics.pairwise import aa_joint_freq, compute_sca_matrix, \
    compute_seq_identity
from cocoatree.deconvolution import eigen_decomp, compute_ica, chooseKpos, \
    icList
from cocoatree.randomize import randomization
import matplotlib.pyplot as plt
import numpy as np

# %%
# Import MSA
# ----------

seq_id, sequences = load_MSA("data/s1Ahalabi_1470.an", format="fasta")
Npos, Nseq = len(sequences[0]), len(sequences)

# %%
# MSA filtering
# -------------
#
# Filter overly gapped positions
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

filt_seqs, pos_kept = filter_gap_pos(sequences, threshold=0.4)
Npos_kept = len(pos_kept)

# %%
# Filter overly gapped sequences
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

seq_id_kept, seq_kept = filter_gap_seq(seq_id, filt_seqs, threshold=0.2,
                                       filtrefseq=False)

# %%
# Compute the matrix of pairwise sequence identity
# --------------------------------------------

sim_matrix = compute_seq_identity(seq_kept)

fig = plt.figure()
plt.rcParams['figure.figsize'] = 10, 10
plt.xlabel('Sequences', fontsize=10)
plt.ylabel(None)
plt.title('Matrix of pairwise sequence identity')
plt.imshow(sim_matrix, vmin=0, vmax=1, cmap='inferno')
plt.colorbar(shrink=0.7)
plt.show()

# %%
# Compute sequence weights
weights, Neff = seq_weights(sim_matrix)
# %%
# compute allele frequencies
aa_freq = aa_freq_at_pos(seq_kept, lambda_coef=0.03, weights=weights)
# %%
# Compute background frequencies
qa = background_freq(aa_freq)
# %%
# Compute joint allele frequencies
fijab, fijab_ind = aa_joint_freq(seq_kept, weights=weights, lambda_coef=0.03)

# %%
# Compute the SCA coevolution matrix
# ----------------------------------

Cijab_raw, Cij = compute_sca_matrix(joint_freqs=fijab,
                                    joint_freqs_ind=fijab_ind,
                                    aa_freq=aa_freq,
                                    qa=qa)

fig = plt.figure()
plt.rcParams['figure.figsize'] = 10, 10
plt.xlabel('Residue', fontsize=10)
plt.ylabel(None)
plt.title('Coevolution matrix')
plt.imshow(Cij, vmin=0, vmax=1.4, cmap='inferno')
plt.colorbar(shrink=0.7)
plt.show()

# %%
# Decomposition of the matrix into principal components
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

eigenvalues, eigenvectors = eigen_decomp(Cij)

# %%
# Plot distribution of eigenvalues
plt.figure()
plt.hist(eigenvalues, bins=100, color="black")
plt.ylabel('Number')
plt.xlabel('Eigenvalue')
plt.show()

# %%
# Select number of significant components
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vrand, Lrand = randomization(sequences, Nrep=10,
                             weights=weights, lbda=0.03, kmax=10, metric='SCA',
                             correction=None)
kpos = chooseKpos(eigenvalues, Lrand)
print('kpos = ' + str(kpos))

plt.rcParams['figure.figsize'] = 9, 4
hist0, bins = np.histogram(Lrand.flatten(), bins=Npos_kept,
                           range=(0, eigenvalues.max()))
hist1, bins = np.histogram(eigenvalues, bins=Npos_kept,
                           range=(0, eigenvalues.max()))
plt.bar(bins[:-1], hist1, np.diff(bins), color='k')
plt.plot(bins[:-1], hist0/10, 'r', linewidth=3)
plt.tick_params(labelsize=11)
plt.xlabel('Eigenvalues', fontsize=18)
plt.ylabel('Numbers', fontsize=18)

# %%
# Independent component analysis (ICA)
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vica, W = compute_ica(eigenvectors, kmax=kpos, learnrate=0.1,
                      iterations=100000)

# Plot results
EVs = eigenvectors
ICs = Vica
pairs = [[x, x+1] for x in range(0, kpos, 2)]
ncols = len(pairs)
plt.rcParams['figure.figsize'] = 14, 8
for k, [k1, k2] in enumerate(pairs):
    plt.subplot(2, ncols, k+1)
    plt.plot(EVs[:, k1], EVs[:, k2], 'ok')
    plt.xlabel("EV%i" % (k1+1), fontsize=16)
    plt.ylabel("EV%i" % (k2+1), fontsize=16)
    plt.subplot(2, ncols, k+1+ncols)
    plt.plot(ICs[:, k1], ICs[:, k2], 'ok')
    plt.xlabel("IC%i" % (k1+1), fontsize=16)
    plt.ylabel("IC%i" % (k2+1), fontsize=16)
plt.tight_layout()

# %%
# Select residues that significantly contribute to each independent component
ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits = icList(Vica, kpos, Cij,
                                                              p_cut=0.95)

print('Sizes of the ' + str(kpos) + ' ICs: ' + icsize)

# %%
# Plot coevolution within and between the sectors
plt.rcParams['figure.figsize'] = 10, 10
plt.subplot(121)
plt.imshow(Cij[np.ix_(sortedpos, sortedpos)], vmin=0, vmax=2,
           interpolation='none', aspect='equal',
           extent=[0, sum(icsize), 0, sum(icsize)], cmap='inferno')
plt.colorbar(shrink=0.25)
line_index = 0
for i in range(kpos):
    plt.plot([line_index + icsize[i], line_index + icsize[i]],
             [0, sum(icsize)], 'w', linewidth=2)
    plt.plot([0, sum(icsize)], [sum(icsize) - line_index,
                                sum(icsize) - line_index], 'w', linewidth=2)
    line_index += icsize[i]

# %%
# Export fasta files of the sectors for all the sequences
seq_id, sequences

sector_1_pos = list(ics[0].items)
sector_1 = []
for sequence in range(len(seq_id)):
    seq = ''
    for pos in sector_1_pos:
        seq += sequences[sequence][pos]
    sector_1.append(seq)

export_fasta(sequences, seq_id, outpath)
# outpath needs to be defined properly
