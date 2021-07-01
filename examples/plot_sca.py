"""
======================
SCA
======================

Here's an example of running SCA pipeline on a test dataset
"""

from akasthesia.coevolution import Alignment, CoevolutionMatrix
import numpy as np
import matplotlib.pyplot as plt

def draw_ic(mat,nodes,nodes_color = 'g', legend=False, ics=[], ax = None):
    nodes_ic = mat[nodes,:]
    non_nodes_ic = np.delete(mat, nodes, axis=0)
    if ax == None:
        figure, ax = plt.subplots()
    ax.scatter(x=non_nodes_ic[:,0], y=non_nodes_ic[:,1], color='r',label = 'Normal residues');
    ax.scatter(x=nodes_ic[:,0], y=nodes_ic[:,1], color=nodes_color, label = 'Sector');
    if not ics == []:
        ax.set_xlabel('IC'+str(ics[0]))
        ax.set_ylabel('IC'+str(ics[1]))
    if legend: ax.legend()

s = [
        [27, 30, 35, 49, 50, 51, 59, 60, 61, 62, 63, 102, 128, 129, 130, 143, 154, 168, 177, 179, 180, 181, 194, 199, 204, 216, 217],
        [36, 148, 149, 158, 162, 163, 166, 169, 170, 174, 175, 178, 193, 195, 196, 200, 201, 205, 206, 207, 209],
        [32, 33, 34, 53, 57, 70, 71, 73, 79, 83, 104, 106, 111, 112, 117, 118, 183],
        [37, 38, 39, 41, 47, 56, 127, 140, 144, 182],
        [89, 91, 92, 94, 95],
        [54, 58, 101, 103, 105, 208, 213]
    ]


# Reading from files
alg = Alignment.from_file("data/s1Ahalabi_1470.an")
print("Loaded "+str(alg.seq_count())+" sequence of size "+str(alg.seq_len()))
print("After trimming:",str(np.shape(alg.filtered_alignment())[0]),"sequences of length",str(np.shape(alg.filtered_alignment())[1]))

###############################################################################
# First: Calculate Weights - Sequence Weight and Position Weights

sca_weights = alg.sca_seq_weights()
pos_weights = alg.position_weights_sca(sca_weights)

###############################################################################
# Sanity check: Plot them out

plt.plot(sca_weights, linewidth=0, marker='+');

###############################################################################
#  Coevolution matrix - weighted

weighted_cij = alg.sca_coevolution(seq_weights=sca_weights,pos_weights=pos_weights)

figure, ax = plt.subplots()
figure.set_size_inches(20,10)
im = ax.matshow(weighted_cij.matrix)
ax.set_yticks(range(0,np.shape(weighted_cij.matrix)[0],40))
ax.set_xticks(range(0,np.shape(weighted_cij.matrix)[0],40))
cbar = ax.figure.colorbar(im, ax=ax)

plt.show()

###############################################################################
#  Perform ICA

vica, w = CoevolutionMatrix.ICA(weighted_cij.matrix)

figure, ax = plt.subplots(2,3)
figure.set_size_inches(15,10)
# draw_ic(np.array([vica[:,0],vica[:,1]]).T,nodes=s[0], ax=ax[0,0])
ax[0,0].scatter(x=vica[:,0], y=vica[:,1]);
ax[0,1].scatter(x=vica[:,2], y=vica[:,3]);
ax[0,2].scatter(x=vica[:,4], y=vica[:,5]);
ax[1,0].scatter(x=vica[:,1], y=vica[:,2]);
ax[1,1].scatter(x=vica[:,3], y=vica[:,4]);
ax[1,2].scatter(x=vica[:,0], y=vica[:,5]);

###############################################################################
#  IC Histogram

figure, axes = plt.subplots(2,3)
figure.set_size_inches(15,10)

i=0
for ax1 in axes:
    ax1[0].set_ylabel("Count")
    for ax in ax1:
        ax.hist(vica[:,i],bins=20,rwidth=0.95);
        i+=1
        ax.set_xlabel("IC"+str(i))