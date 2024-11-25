import os
from ..cocoatree_io import load_MSA


def load_S1A_serine_proteases():
    """
    Load the S1A serine protease dataset

    Note that this is not the original dataset from Halabi et al,
    Cell, 2008. The dataset has been augmented with more recent data
    (hence the MSA is larger), and reprocessed.

    Returns
    -------

    a dictionnary containing :
        - `sequences_ids`: a list of strings corresponding to sequence names
        - `alignment`: a list of strings corresponding to sequences. Because it
          is an MSA, all the strings are of same length.
    """

    module_path = os.path.dirname(__file__)
    filename = os.path.join(
        module_path,
        "data/S1A_serine_proteases/alignment.fasta")
    sequence_ids, sequences = load_MSA(filename, format="fasta")

    return {"sequence_ids": sequence_ids,
            "alignment": sequences}
