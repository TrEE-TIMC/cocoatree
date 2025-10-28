"""
======================
A minimal SCA analysis
======================

This example provides a minimal code for identifying sectors -- more precisely,
XCoRs :cite:p:`jullien_cocoatree_2025` -- associated with a MSA.
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
# dataset, you can use the :func:`cocoatree.io.load_msa` function.
#
# For more details on the S1A serine proteases dataset, go to
# :ref:`sphx_glr_auto_examples_d_datasets_plot_s1A_serine_proteases.py`.

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
# Extracting the list of residues (i.e. positions) composing the XCoRs.

print('XCoR_1:', [])
print('XCoR_2:', [])
print('XCoR_3:', [])

# %%
# Visualizing XCoRs on independent components
# -------------------------------------------

# Plotting independent components
fig, ax = plt.subplots()
ax.plot(results.loc[:, "IC1"],
        results.loc[:, "IC2"],
        ".", c="black")

# Plotting XCoRs elements
for isec, color in zip([1, 2, 3], ['r', 'g', 'b']):
    ax.plot(results.loc[results["xcor_%d" % isec], "IC1"],
            results.loc[results["xcor_%d" % isec], "IC2"],
            ".", c=color, label="XCoR_%d" % isec)

ax.set_xlabel("IC1")
ax.set_ylabel("IC2")

ax.legend()
