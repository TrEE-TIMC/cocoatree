"""
==============================
The simplest SCA analysis ever
==============================

This example shows the full process to perform a complete SCA analysis and
detect protein sectors from data importation, MSA filtering.
"""

import cocoatree.datasets as c_data
import cocoatree
import matplotlib.pyplot as plt

# %%
# Load the dataset
# ----------------
#
# We start by importing the dataset. In this case, we can directly load the S1
# serine protease dataset provided in :mod:`cocoatree`. To work on your on
# dataset, you can use the `cocoatree.io.load_msa` function.

serine_dataset = c_data.load_S1A_serine_proteases()
loaded_seqs = serine_dataset["alignment"]
loaded_seqs_id = serine_dataset["sequence_ids"]
n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)


# %%
# Compute the SCA analysis
# ------------------------
#

coevol_matrix, results = cocoatree.perform_sca(
    loaded_seqs_id, loaded_seqs, n_components=3)
print(results.head())

# %%
# Select the position of the first sector in the results dataframe.

print(results.loc[results["sector_1"]].head())

# %%
# Visualizing the sectors on the first and second IC
# --------------------------------------------------

fig, ax = plt.subplots()
# Start by plotting elements that are in sector 1
ax.plot(results.loc[:, "PC1"],
        results.loc[:, "PC2"],
        ".",
        c="black")

# Start by plotting elements that are in sector 1
ax.plot(results.loc[results["sector_1"], "PC1"],
        results.loc[results["sector_1"], "PC2"],
        ".",
        c="r", label="Sector 1")

ax.plot(results.loc[results["sector_2"], "PC1"],
        results.loc[results["sector_2"], "PC2"],
        ".",
        c="g", label="Sector 2")


ax.plot(results.loc[results["sector_3"], "PC1"],
        results.loc[results["sector_3"], "PC2"],
        ".",
        c="b", label="Sector 3")

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

ax.legend()
