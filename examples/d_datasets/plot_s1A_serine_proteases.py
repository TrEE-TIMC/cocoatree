"""
S1A serine proteases
====================

This example shows how the S1A serine proteases dataset is organized and how
you can use it to test cocoatree's features.

This dataset is actually composed of two MSAs:
    - the MSA from Halabi et al., (2008)
    - the MSA from Rivoire et al., (2016)

"""

# %%
# Imports
import numpy as np
from cocoatree.datasets import load_S1A_serine_proteases

import cocoatree.msa as c_msa
import cocoatree.statistics.position as c_pos

# %%
# Halabi dataset
# --------------
# Use `load_S1A_serine_dataset` function to load the Halabi dataset
halabi = load_S1A_serine_proteases(paper='halabi')

print(halabi.keys())

# %%
# `halabi` is a dictionary with 6 elements:
#   - `sequence_ids`: a list of the sequence identifiers for the MSA
#   - `alignment`: the aligned sequences
#   - `sector_positions`: a Numpy NpzFile of the sectors' positions as found in
# the Halabi study
#   - `metadata`: a Pandas dataframe containing metadata associated with the
# alignment, such as taxonomic information.
#   - `pdb_sequence`: a tuple of the pdb sequence of 3TGI (rat's trypsin)
#   - `pdb_positions`: a list of the named positions found in the PDB
#
# Access the sequence data
# ^^^^^^^^^^^^^^^^^^^^^^^^
loaded_seqs = halabi["alignment"]
loaded_seqs_id = halabi["sequence_ids"]
print(loaded_seqs[:5])
print(loaded_seqs_id[:5])
print()
print("Number of sequences:", len(loaded_seqs))

n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)
print(f"The loaded MSA has {n_loaded_seqs} sequences and {n_loaded_pos} \
      positions.")

# %%
# Filter the MSA
# ^^^^^^^^^^^^^^
sequences, sequences_id, positions = c_msa.filter_sequences(
    loaded_seqs, loaded_seqs_id, gap_threshold=0.4, seq_threshold=0.2)
n_pos = len(positions)
print(f"After filtering, we have {n_pos} remaining positions.")
print(f"After filtering, we have {len(sequences)} remaining sequences.")

# %%
# Compute the number of effective sequences
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
seq_weights, m_eff = c_pos.compute_seq_weights(sequences)
print('Number of effective sequences %d' % np.round(m_eff))
print()

# %%
# Access the metadata
# ^^^^^^^^^^^^^^^^^^^
metadata = halabi["metadata"]
print(metadata)

# %%
# The `Seq_ID` column corresponds to the sequence identifiers found in
# `loaded_seqs_id`
#
# Access PDB information
# ^^^^^^^^^^^^^^^^^^^^^^
pdb_seq = halabi["pdb_sequence"]
pdb_pos = halabi["pdb_positions"]
print("Start of 3TGI sequence:\n", pdb_seq[:10])
print("Named positions of the first ten residues in 3TGI protein sequence:\n",
      pdb_pos[:10])

# %%
# Access sector information
# ^^^^^^^^^^^^^^^^^^^^^^^^^
sectors = halabi["sector_positions"]
print(sectors)
# %%
# There are three sectors that have been identified in Halabi et al. (2008)
print("The first sector identified in Halabi et al. is composed of the \
      following positions:", sectors["sector_1"])

# %%
# To see an example on how to use this dataset for a coevolution analysis, go
# to: :ref:`sphx_glr_auto_examples_plot_full_SCA_analysis.py`
#
# For an example on how to visualize cocoatree's results with this dataset, go
# to : :ref:`sphx_glr_auto_examples_plot_sector_along_tree_and_metadata.py`


# %%
# Rivoire dataset
# ---------------
# Similarly, you can access the dataset from Rivoire et al. (2016) paper:
rivoire = load_S1A_serine_proteases(paper='rivoire')
print(rivoire.keys())

loaded_seqs = rivoire["alignment"]
loaded_seqs_id = rivoire["sequence_ids"]
print("Number of sequences:", len(loaded_seqs))

n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)
print(f"The loaded MSA has {n_loaded_seqs} sequences and {n_loaded_pos} \
      positions.")

# %%
# The dataset is organized in the same way as Halabi's, although in this case
# there are 6 sectors:
sectors = rivoire["sector_positions"]
print(sorted(sectors))

for sect in sectors:
    print(sect, "is composed of:", len(sectors[sect]), "positions.")


# %%
# References
# ----------
#
# .. [1] `Protein Sectors: Evolutionary Units of Three-Dimensional Structure
#        <https://doi-org.insb.bib.cnrs.fr/10.1016/j.cell.2009.07.038>`_,
#        Halabi et al., Cell 2008
# .. [2] `Evolution-Based Functional Decomposition of Proteins
#        <https://doi.org/10.1371/journal.pcbi.1004817>`_,
#        Rivoire et al., PLOS Computational Biology 2016
