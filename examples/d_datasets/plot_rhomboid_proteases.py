"""
Rhomboid proteases
====================

This example shows how the Rhomboid dataset is organized and how you can use it
to test cocoatree's features.

"""

# %%
# Load the dataset
import numpy as np
from cocoatree.datasets import load_rhomboid_proteases
import cocoatree.msa as c_msa
import cocoatree.statistics.position as c_pos

dataset = load_rhomboid_proteases()
print(dataset.keys())

# %%
# `dataset` is a dictionary with 6 elements:
#   - `sequence_ids`: a list of the sequence identifiers for the MSA
#   - `alignment`: the aligned sequences
#   - `sector_positions`: a Numpy NpzFile of the sectors' positions as found in
# the Mihaljević and Urban study
#   - `metadata`: a Pandas dataframe containing metadata associated with the
# alignment, such as taxonomic information.
#   - `pdb_sequence`: a tuple of the 2NRF pdb sequence from E. coli
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
# Note that the above sector positions are relative to E. coli's PDB sequence.

# %%
# Access metadata
# ^^^^^^^^^^^^^^^
metadata = dataset["metadata"]
print(metadata)

# %%
# Access PDB information
# ^^^^^^^^^^^^^^^^^^^^^^
pdb_seq = dataset["pdb_sequence"]
pdb_pos = dataset["pdb_positions"]
print("Start of 2NRF E. coli's sequence:\n", pdb_seq[:10])
print("Named positions of the first ten residues in E. coli's sequence:\n",
      pdb_pos[:10])


# %%
# References
# ----------
#
# .. [1] `Decoding the Functional Evolution of an Intramembrane Protease
#        Superfamily by Statistical Coupling Analysis
#        <https://doi.org/10.1016/j.str.2020.07.015>`_,
#        Mihaljević and Urban, Structure 2020
