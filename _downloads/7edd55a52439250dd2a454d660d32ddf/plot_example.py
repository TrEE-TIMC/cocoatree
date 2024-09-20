"""
======================
A magnificient example
======================

Here's an example with maplotlib and numpy
"""

import numpy as np
import matplotlib.pyplot as plt
from cocoatree import coevolution

X = np.arange(10)

###############################################################################
# I could add another code block

fig, ax = plt.subplots()
ax.plot(np.arange(10))
ax.set_xlabel("X-label")
ax.set_ylabel("Y-label")

