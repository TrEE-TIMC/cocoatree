"""
=============================================================================
Generating a simulated dataset with cocoa's datasets submodule - reformulated
=============================================================================

Reformulated version
====================

This example showcases how to generate a simulated datasets using the datasets
submodule.

In this document we use a reformulated model. The reformulated model take a positional
background frequency and w to recalculate v, which compenstated the impact of w on
the distribution of amino acid at position.
The reformulated v is

.. math::
  v_i(a) = log(f_i(a)) + log Z - \\lambda \\sum_c \\sum_{j \\neq i} w_{ij}(a, c)

"""

import numpy as np
import matplotlib.pyplot as plt

from akasthesia.datasets import simulation
import __helper


##########################################
# Specify the parameters of the simulation
##########################################

N = 100
L = 75
n_nodes = 7  # Each sectors consist of 7 nodes
n_secs = 2   # 2 sectors overall
seed = 4

######################
# Generating the model
######################
# The model is defined by 2 parameters :math:`\mathbf{v}` and :math:`\mathbf{w}`
# The probability of a certain sequence is
#
# .. math::
#   p(x|v,w) \propto (\sum_i v_i + \sum_{i,j} w_{i,j})
#
# Here we specify a positional frequency for amino acids, which are uniformal
# for all amino acid at every residue.
#
# %%

np.random.seed(seed)
edges = []
nodes = []
freq_ai = np.tile(.05, [L, 20])
for _ in range(n_secs):
    s_node = np.random.randint(L-n_nodes)
    _nodes = []
    table = [np.random.permutation(20) for i in range(n_nodes)]
    for k in range(20):
        _nodes += [[i, table[i-s_node][k]]
                   for i in range(s_node,  s_node+n_nodes)]
    nodes.append(_nodes)
    edges += __helper.generate_topology(_nodes, topo='line')

w = __helper.generate_w(L, edges=edges)
v = simulation.recompute_v(freq_ai, w, lmbd=1)


########################
# Generate the data sets
########################
# The simulation take the model v and w; and perform T Gibbs cycle
#
# The initial sequences are drawed randomly

T = 10
msa_raw = simulation.gibbs_sampling(np.random.randint(20, size=[L]), N, v, w, T, burnin=500)
alg = simulation.to_Alignment(msa_raw)

###############################################################################
# Check if the simulation is correct

assert alg.seq_count() == N
assert alg.seq_len() == L


##############################################################################
# Plot the positional entropy of the data set: sectors should the same entropy
# with the rest of the residues now (See other examples for un)

# %%
color = ['red']*L
for i in edges:
    color[i[0]] = 'blue'
    color[i[2]] = 'blue'

fig, ax = plt.subplots()
fig.set_size_inches(10, 4)
ax.bar(np.arange(L), alg.vorberg_entropy(), color=color)
ax.set_xlabel('Residues')
ax.set_ylabel('Entropy')

############################
# Plot the overall histogram

fig, ax = plt.subplots()
fig.set_size_inches(10, 3)
hist = np.histogram(alg.filtered_alignment().flatten()-1, bins=np.arange(21))[0]/(N*L)
ax.bar(np.arange(20), hist)
ax.set_xlabel("Amino Acid")
ax.set_ylabel("Frequency")

##############################################################################
# Plot the coevolution matrix of the data set: sectors should show up brighter
# compare to the rest of the residues
# %%
weighted_cij = alg.sca_coevolution()

figure, ax = plt.subplots()
im = ax.matshow(weighted_cij.matrix)
ax.set_yticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
ax.set_xticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
cbar = ax.figure.colorbar(im, ax=ax)
