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
# Visualizing the sectors on the first and second PC
# --------------------------------------------------

# Plotting all elements in components
fig, ax = plt.subplots()
ax.plot(results.loc[:, "PC1"],
        results.loc[:, "PC2"],
        ".", c="black")

# Plotting elements in sectors
for isec, color in zip([1, 2, 3], ['r', 'g', 'b']):
    ax.plot(results.loc[results["sector_%d" % isec], "PC1"],
            results.loc[results["sector_%d" % isec], "PC2"],
            ".", c=color, label="Sector %d" % isec)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

ax.legend()


# %%
# Visualizing the sectors on the first and second IC
# --------------------------------------------------

# Plotting all elements in components
fig, ax = plt.subplots()
ax.plot(results.loc[:, "IC1"],
        results.loc[:, "IC2"],
        ".", c="black")

# Plotting elements in sectors
for isec, color in zip([1, 2, 3], ['r', 'g', 'b']):
    ax.plot(results.loc[results["sector_%d" % isec], "IC1"],
            results.loc[results["sector_%d" % isec], "IC2"],
            ".", c=color, label="Sector %d" % isec)

ax.set_xlabel("IC1")
ax.set_ylabel("IC2")

ax.legend()
