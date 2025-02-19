"""
==============================
The simplest SCA analysis ever
==============================

This example shows the full process to perform a complete SCA analysis and
detect protein sectors from data importation, MSA filtering.
"""

import cocoatree.datasets as c_data
import cocoatree

# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset. In this case, we can directly load the S1
# serine protease dataset provided in :mod:`cocoatree`. To work on your on
# dataset, you can use the `cocoatree.io.load_msa` function.

serine_dataset = c_data.load_S1A_serine_proteases('rivoire')
loaded_seqs = serine_dataset["alignment"]
loaded_seqs_id = serine_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)


# %%
# Compute the SCA analysis
# ------------------------
#

coevol_matrix, results = cocoatree.perform_sca(loaded_seqs_id, loaded_seqs)
print(results.head())

# %%
# Select the position of the first sector

print(results.loc[results["sector_1"]].head())
