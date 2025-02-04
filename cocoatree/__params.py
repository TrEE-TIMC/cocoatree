# Parameters for COCOA-Tree

import numpy as np

# in pySCA, the defautl value is 0.03 (default here)
# in mean-field DCA, the default value is 0.5
__freq_regularization_ref = 0.03

__freq0 = np.array(
    [
        0.073,
        0.025,
        0.050,
        0.061,
        0.042,
        0.072,
        0.023,
        0.053,
        0.064,
        0.089,
        0.023,
        0.043,
        0.052,
        0.040,
        0.052,
        0.073,
        0.056,
        0.063,
        0.013,
        0.033
    ]
)

lett2num = {
    '-': 0,
    'A': 1,
    'C': 2,
    'D': 3,
    'E': 4,
    'F': 5,
    'G': 6,
    'H': 7,
    'I': 8,
    'K': 9,
    'L': 10,
    'M': 11,
    'N': 12,
    'P': 13,
    'Q': 14,
    'R': 15,
    'S': 16,
    'T': 17,
    'V': 18,
    'W': 19,
    'Y': 20}

__aa_count = len(lett2num)

aatable = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
