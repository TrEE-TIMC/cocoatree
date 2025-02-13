"""
====================================================
Mapping original MSA, filtered MSA, PDB, and sectors
====================================================

In this example, we showcase how to create a pandas.DataFrame to map the
original MSA's positions with the PDB positions, PDB named position, the MSA
filtered positions, and the sectors.
"""

import numpy as np
from cocoatree.datasets import load_S1A_serine_proteases
from cocoatree.msa import filter_sequences
from cocoatree.msa import map_msa_positions
import pandas as pd


serine_dataset = load_S1A_serine_proteases(paper="halabi")
seq_id = serine_dataset["sequence_ids"]
sequences = serine_dataset["alignment"]
n_pos, n_seq = len(sequences[0]), len(sequences)
sectors = [
    [int(i) for i in serine_dataset["sector_positions"][key]]
    for key in serine_dataset["sector_positions"].keys()]
pdb_pos = serine_dataset["pdb_positions"]

seq_kept, seq_id_kept, pos_kept = filter_sequences(sequences, seq_id)


pos_mapping, _ = map_msa_positions(n_pos, pos_kept)


# Sectors are in the PDB referential. The sequence corresponding to the PDB is
# the first of the MSA.
is_mapped = np.array([s != "-" for s in sequences[0]])
pdb_mapping = [int(val) if f else None
               for f, val in zip(
               is_mapped, (is_mapped.cumsum()-1))]
pdb_pos_mapping = [
    pdb_pos[j]
    if i else None
    for i, j in zip(is_mapped, is_mapped.cumsum()-1)]
mapping = pd.DataFrame(
    {"original_msa_pos": np.arange(n_pos, dtype=int),
     "pdb_pos": pdb_mapping,
     "pdb_named_pos": pdb_pos_mapping,
     "pdb_pos_trunc": [int(i)
                       if i is not None and len(i) < 4 else None
                       for i in pdb_pos_mapping],
     "filtered_msa_pos": pos_mapping.values()})

mapping["sector_1"] = np.isin(mapping["pdb_pos_trunc"], sectors[0])
mapping["sector_2"] = np.isin(mapping["pdb_pos_trunc"], sectors[1])
mapping["sector_3"] = np.isin(mapping["pdb_pos_trunc"], sectors[2])
mapping = mapping.drop("pdb_pos_trunc", axis=1)

# %%
# Now print the indices of sectors 1 in the different referentials

print(mapping.loc[mapping["sector_1"]].head())
