"""
DHFR proteases
====================

This example shows how the DHFR dataset is organized and how you can use it to
test cocoatree's features.

"""

import numpy as np

from cocoatree.datasets import load_DHFR
import cocoatree.msa as c_msa
import cocoatree.statistics.position as c_pos

# %%
# Load the dataset
# ^^^^^^^^^^^^^^^^
dataset = load_DHFR()
print(dataset.keys())

# %%
# `dataset` is a dictionary with 5 elements:
#   - `sequence_ids`: a list of the sequence identifiers for the MSA
#   - `alignment`: the aligned sequences
#   - `sector_positions`: a Numpy NpzFile of the sectors' positions as found in
# the Kalmer study
#   - `pdb_sequence`: a tuple of the pdb sequence of 3QL3 (wild-type E. coli
# DHFR)
#   - `pdb_positions`: a list of the named positions found in the PDB
#
# Access the sequence data
# ^^^^^^^^^^^^^^^^^^^^^^^^
loaded_seqs = dataset["alignment"]
loaded_seqs_id = dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)

print(f"The loaded MSA has {n_loaded_seqs} sequences and {n_loaded_pos} \
      positions.")

# %%
# Example of sequence filtering
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
sequences, sequences_id, positions = c_msa.filter_sequences(
    loaded_seqs, loaded_seqs_id, gap_threshold=0.4, seq_threshold=0.2)
n_pos = len(positions)
print(f"After filtering, we have {n_pos} remaining positions.")
print(f"After filtering, we have {len(sequences)} remaining sequences.")

seq_weights, m_eff = c_pos.compute_seq_weights(sequences)
print('Number of effective sequences %d' %
      np.round(m_eff))

# %%
# Access sector information
# ^^^^^^^^^^^^^^^^^^^^^^^^^
sectors = dataset["sector_positions"]
print("There are", len(sectors), "sectors")

for sect in sectors:
    print(sect, "is composed of:", len(sectors[sect]), "positions.")

# %%
# Access PDB information
# ^^^^^^^^^^^^^^^^^^^^^^
pdb_seq = dataset["pdb_sequence"]
pdb_pos = dataset["pdb_positions"]
print("Start of E. coli DHFR's sequence:\n", pdb_seq[:10])
print("Named positions of the first ten residues in E. coli's sequence:\n",
      pdb_pos[:10])

# %%
# References
# ----------
#
# .. [1] `Statistical Coupling Analysis Predicts Correlated Motions in
#        Dihydrofolate Reductase
#        <https://pubs.acs.org/doi/10.1021/acs.jpcb.4c04195>`_,
#        Kalmer et al., J. Phys. Chem. 2024
