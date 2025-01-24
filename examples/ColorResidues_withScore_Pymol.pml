# Copy paste the following line in pymol
# @path_to_this_pmlCode

python
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
# from cocoatree.visualization import _generate_colors_from_colormaps


def _generate_colors_from_colormaps(n_colors, cmap="jet", as_hex=True):
    """
    Generate a list of n colors from colormap
    """

    colormap = plt.get_cmap(str(cmap))
    indx = np.linspace(0, 1, n_colors)
    indexed_colors = [colormap(i) for i in indx]
    if as_hex:
        indexed_colors = [colors.to_hex(i) for i in indexed_colors]
    return indexed_colors


residues_swap = list(np.load('/home/jullimar/Documents/Postdoc_TIMC/test_color_residues.npy').astype('int'))
Scores =  list(np.load('/home/jullimar/Documents/Postdoc_TIMC/test_color_residues_contributions.npy'))

nvals = len(residues_swap)
colormap = _generate_colors_from_colormaps(nvals, cmap='jet', as_hex=False)

for i in range(len(residues_swap)):
    residue = residues_swap[i]
    r = colormap[i][0]
    g = colormap[i][1]
    b = colormap[i][2]

    color_name = f"residue_{residue}"
    cmd.set_color(color_name, [r, g, b])
    cmd.color(color_name, f"resi {residue}")
    cmd.show("spheres", f"resi {residue}")
python end
