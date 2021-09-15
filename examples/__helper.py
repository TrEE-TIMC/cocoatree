import numpy as np

aa = '-ACDEFGHIKLMNPQRSTVWY'


# Generate parameters
def generate_v(L, type_of_v, exclusion=[], seed=42):
    # Define v
    rng = np.random.default_rng(seed)
    if type_of_v == 'simple':
        v = np.ones([L, 20])
    else:
        choices = np.delete(np.arange(L), exclusion)
        v_highrange = [8,  9]
        v_lowrange = [-0.02,  0.05]
        high_v_residue = rng.choice(choices, size=20, replace=False)
        low_v_residue = [i for i in range(L) if i not in high_v_residue]
        v = np.ones([L, 20])
        for i in high_v_residue:
            v[i] = 0
            for j in np.random.choice(range(0, 20), 3):
                v[i, j] = np.random.uniform(v_highrange[0], v_highrange[1])
        for i in low_v_residue:
            v[i] = v[i]*np.random.uniform(v_lowrange[0], v_lowrange[1])
    return v


def generate_topology(nodes,  topo,  n_edges=None):
    # edge nomenclature: [i, a, j, b]: amino acid a at res i to aa b at j
    edges = []
    r = len(nodes)
    if topo == 'line':
        for i in range(r-1):
            if nodes[i][0] < nodes[i+1][0]:
                edges.append(nodes[i]+nodes[i+1])

    elif topo == 'circle':
        for i in range(r-1):
            edges.append(nodes[i]+nodes[i+1])
        edges.append(nodes[-1]+nodes[0])

    elif topo == 'fully-connected':
        for i in range(r):
            for j in range(i+1, r):
                if nodes[i][0] != nodes[j][0]:
                    edges.append(nodes[i]+nodes[j])

    elif topo == 'randomly-connected':
        # 50% of fully connected edges are removed
        if n_edges is None:
            print('number of edges is required for random connection')
            return None
        for i in range(r):
            for j in range(i+1, r):
                edges.append(nodes[i]+nodes[j])
        edges = [edges[retain] for retain in
                 [i for i in np.random.default_rng().choice(range(len(edges)),
                  size=n_edges, replace=False)]]
    return edges


def generate_w(L, edges):
    # edge nomenclature: [i, a, j, b]: amino acid a at res i to aa b at j
    w = np.zeros([L, L, 20, 20])
    for edge in edges:
        if edge[0] == edge[2]:
            continue
        w[edge[0], edge[2], edge[1], edge[3]] = 5  # a@i to b@j
        w[edge[2], edge[0], edge[3], edge[1]] = 5  # b@j to a@b
    return w


def generate_v_and_w(L, v_type, edges):
    # Create custom binary matrices
    exclusion = []
    for edge in edges:
        exclusion += [edge[0], edge[2]]
    x_single = generate_v(L, v_type, exclusion)
    x_pair = generate_w(L, edges)
    return x_single, x_pair

##
#
# v = (np.tile(np.log(z), (20, 1)).T + np.log(fia) - .25*np.sum(w, axis=(1, 3)))
