import os
from ..io import load_MSA
import numpy as np
import pandas as pd


def load_S1A_serine_proteases(paper):
    """
    Load the S1A serine protease dataset

    Note that this is not the original dataset from Halabi et al,
    Cell, 2008. The dataset has been augmented with more recent data
    (hence the MSA is larger), and reprocessed.

    Parameters
    ----------
    paper: str, either 'halabi' or 'rivoire'
        whether to load the dataset from Halabi et al, Cell, 2008 or from
        Rivoire et al, PLoS Comput Biol, 2016

    Returns
    -------

    a dictionnary containing :
        - `sequences_ids`: a list of strings corresponding to sequence names
        - `alignment`: a list of strings corresponding to sequences. Because it
          is an MSA, all the strings are of same length.
        - `metadata`: a pandas dataframe containing the metadata associated
          with the alignment.
        - `sector_positions`: a dictionnary of arrays containing the residue
          positions associated to sectors, either in Halabi et al, or in
          Rivoire et al.
    """

    module_path = os.path.dirname(__file__)

    if paper == 'halabi':
        # Load the alignment used in Halabi et al, 2008
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/halabi_alignment.fasta")
        sequence_ids, sequences = load_MSA(filename, format="fasta")
        # Load the positions of the 3 sectors identified in Halabi et al, Cell,
        # 2008
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/halabi_sectors.npz")
        sectors = np.load(filename)

    if paper == 'rivoire':
        # Load the alignment used in Rivoire et al, 2016
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/rivoire_alignment.fasta")
        sequence_ids, sequences = load_MSA(filename, format="fasta")
        # Load the positions of the 6 sectors identified in Rivoire et al, PLoS
        # Comput Biol, 2016
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/rivoire_sectors.npz")
        sectors = np.load(filename)
        # Load the metadata
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/rivoire_metadata.csv")
        metadata = pd.read_csv(filename)

    return {"sequence_ids": sequence_ids,
            "alignment": sequences,
            "sector_positions": sectors,
            "metadata": metadata}


def load_rhomboid_proteases():
    """
    Load the rhomboid protease dataset

    This dataset comes from Mihaljevic & Urban, Cell, 2020
    (DOI: https://doi.org/10.1016/j.str.2020.07.015).

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
        "data/rhomboid_proteases/Data_S1_Rhomboid_MSA.fasta")
    sequence_ids, sequences = load_MSA(filename, format="fasta")

    return {"sequence_ids": sequence_ids,
            "alignment": sequences}
