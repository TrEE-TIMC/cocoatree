"""
==============================================================
Generating a simulated dataset with cocoa's datasets submodule
==============================================================

This example showcases how to generate a simulated datasets using the datasets
submodule.
"""

import numpy as np
import matplotlib.pyplot as plt

from akasthesia.datasets import simulation
import __helper


###############################################################################
# Specify the parameters of the simulation
###############################################################################

N = 100
L = 75
n_nodes = 7  # Each sectors consist of 7 nodes
n_secs = 2   # 2 sectors overall
seed = 42

###############################################################################
# Generating the model
###############################################################################
# The model is defined by 2 parameters :math:`\mathbf{v}` and :math:`\mathbf{w}`
# The probability of a certain sequence is
#
# .. math::
# p(x|v,w) \propto (\sum_i v_i + \sum_{i,j} w_{i,j})

np.random.seed(seed)
edges = []
nodes = []
for _ in range(n_secs):
    s_node = np.random.randint(L-n_nodes)
    _nodes = []
    table = [np.random.permutation(20) for i in range(n_nodes)]
    for __ in range(20):
        _nodes += [[i, table[i-s_node][__]]
                   for i in range(s_node,  s_node+n_nodes)]
    nodes.append(_nodes)
    edges += __helper.generate_topology(_nodes, topo='line')

(v, w) = __helper.generate_v_and_w(L, "simple", edges=edges)


###############################################################################
# Generate the data sets
###############################################################################
# The simulation take the model v and w; and perform T Gibbs cycle
#
# The initial sequences are drawed randomly
# %%

T = 10
msa_raw = simulation.gibbs_sampling(np.random.randint(20, size=[L]), N, v, w, T)
alg = simulation.to_Alignment(msa_raw)

###############################################################################
# Check if the simulation is correct

assert alg.seq_count() == N
assert alg.seq_len() == L


###############################################################################
# Plot the positional entropy of the data set: sectors should have lower entropy
# compare to the rest of the residues

# %%
color = ['red']*L
for i in edges:
    color[i[0]] = 'blue'
    color[i[2]] = 'blue'

fig, ax = plt.subplots()
fig.set_size_inches(10, 3)
ax.bar(np.arange(L), alg.vorberg_entropy(), color=color)
ax.set_xlabel('Residues')
ax.set_ylabel('Entropy')


###############################################################################
# Plot the coevolution matrix of the data set: sectors should show up brighter
# compare to the rest of the residues
# %%
seq_weights = alg.sca_seq_weights()
pos_weights = alg.position_weights_sca(seq_weights)

weighted_cij = alg.sca_coevolution(seq_weights=seq_weights,
                                   pos_weights=pos_weights)

figure, ax = plt.subplots()
im = ax.matshow(weighted_cij.matrix)
ax.set_yticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
ax.set_xticks(range(0, np.shape(weighted_cij.matrix)[0], 40))
cbar = ax.figure.colorbar(im, ax=ax)
###############################################################################
# Further analysis can be done on this data sets, similar to S1A
