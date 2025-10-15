"""
S1A serine proteases
====================

.. image:: images/Chymotrypsin_enzyme.png
    :width: 600px
    :height: 600px
    :alt: 3D structure of bovine chymotrypsin
    :align: center

This example shows how the S1A serine proteases dataset is organized and how you
can use it to test cocoatree's features.

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
# :ref:`sphx_glr_auto_examples_plot_sector_along_tree_and_metadata.py`


for paper in ["halabi", "rivoire"]:
    print(paper)
    dataset = load_S1A_serine_proteases(paper=paper)

    print("Number of sequences", len(dataset["alignment"]))
    loaded_seqs = dataset["alignment"]
    loaded_seqs_id = dataset["sequence_ids"]
    n_loaded_pos, n_loaded_seqs = len(loaded_seqs[0]), len(loaded_seqs)

    print(f"The loaded MSA has {n_loaded_seqs} sequences and {n_loaded_pos} \
        positions.")

    sequences, sequences_id, positions = c_msa.filter_sequences(
        loaded_seqs, loaded_seqs_id, gap_threshold=0.4, seq_threshold=0.2)
    n_pos = len(positions)
    print(f"After filtering, we have {n_pos} remaining positions.")
    print(f"After filtering, we have {len(sequences)} remaining sequences.")

    seq_weights, m_eff = c_pos.compute_seq_weights(sequences)
    print('Number of effective sequences %d' %
          np.round(m_eff))
    print()

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
# 
