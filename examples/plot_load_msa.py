"""
Load an MSA
===========

A small example that shows how to load and filter an MSA, compute positional
and joint amino acid frequencies and finally, compute the SCA coevolution
matrix.

"""


from cocoatree.msa import load_MSA, filter_gap_seq, filter_gap_pos, \
    compute_seq_identity
from cocoatree.sca_functions import seq_weights, aa_freq_at_pos, \
    background_freq, aa_joint_freq, compute_sca_matrix

seq_id, sequences = load_MSA("data/s1Ahalabi_1470.an", format="fasta")
Npos, Nseq = len(sequences[0]), len(sequences)

# Filter overly gapped positions
filt_seqs, pos_kept = filter_gap_pos(sequences, threshold=0.4)

# Filter overly gapped sequences
seq_id_kept, seq_kept = filter_gap_seq(seq_id, filt_seqs, threshold=0.2,
                                       filtrefseq=False)

# Compute matrix of pairwise sequence identity
sim_matrix = compute_seq_identity(seq_kept, graphic=False)

# Compute sequence weights
weights, Neff = seq_weights(sim_matrix)

# compute allele frequencies
fia = aa_freq_at_pos(seq_kept, lbda=0.03, aa_count=21, weights=weights)
# Compute background frequencies
qa = background_freq(fia)
# Compute joint allele frequencies
fijab, fijab_ind = aa_joint_freq(seq_kept, weights=weights, lbda=0.03)

# Compute the SCA coevolution matrix
Cijab_raw, Cij = compute_sca_matrix(joint_freqs=fijab,
                                    joint_freqs_ind=fijab_ind,
                                    fia=fia,
                                    qa=qa)
