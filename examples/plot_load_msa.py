"""
Load an MSA
===========

A small example that shows how to load an MSA

"""


from cocoatree.msa import load_MSA

seq_id, sequences = load_MSA("data/s1Ahalabi_1470.an", format="fasta")
