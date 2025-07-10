"""
S1A serine proteases
====================

Load the dataset

"""

import numpy as np
from cocoatree.datasets import load_S1A_serine_proteases

import cocoatree.msa as c_msa
import cocoatree.statistics.position as c_pos


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
