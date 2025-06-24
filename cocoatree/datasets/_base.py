import os
from ..io import load_MSA, load_pdb
import numpy as np
import pandas as pd
import gzip


def load_S1A_serine_proteases(paper='rivoire'):
    """
    Load the S1A serine protease dataset

    Halabi dataset: 1470 sequences of length 832; 3 sectors identified
    Rivoire dataset : 1390 sequences of length 832 (snake sequences were
    removed for the paper's analysis); 6 sectors identified (including the
    3 from Halabi et al, 2008)

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
          positions associated to each sector, either in Halabi et al, or in
          Rivoire et al.
        - `pdb_sequence`: sequence extracted from rat's trypsin PDB structure
        - `pdb_positions`: positions extracted from rat's trypsin PDB structure
    """

    module_path = os.path.dirname(__file__)

    if paper == 'halabi':
        # Load the alignment used in Halabi et al, 2008
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/halabi_alignment.fasta")
        data = load_MSA(filename, format="fasta")
        # Load the positions of the 3 sectors identified in Halabi et al, Cell,
        # 2008
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/halabi_sectors.npz")
        sectors = np.load(filename)
        # Load the metadata
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/halabi_metadata.csv")
        metadata = pd.read_csv(filename)

    elif paper == 'rivoire':
        # Load the alignment used in Rivoire et al, 2016
        filename = os.path.join(
            module_path,
            "data/S1A_serine_proteases/rivoire_alignment.fasta")
        data = load_MSA(filename, format="fasta")
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

    else:
        raise ValueError(f"invalid paper: {paper}. Options are 'halabi' or \
                         'rivoire'")

    # Load the PDB structure
    filename = os.path.join(
        module_path,
        "data/S1A_serine_proteases/3tgi.pdb")
    pdb_sequence, pdb_positions = load_pdb(filename, '3TGI', 'E')
    data["sector_positions"] = sectors
    data["metadata"] = metadata
    data["pdb_sequence"] = pdb_sequence,
    data["pdb_positions"] = pdb_positions

    return data


def load_rhomboid_proteases():
    """
    Load the rhomboid protease dataset

    This dataset comes from Mihaljevic & Urban, Cell, 2020
    (DOI: https://doi.org/10.1016/j.str.2020.07.015).

    Returns
    -------
    a dictionnary containing :
        - `sequence_ids`: a list of strings corresponding to sequence names
        - `alignment`: a list of strings corresponding to sequences. Because it
          is an MSA, all the strings are of same length.
         - `sector_positions`: a dictionnary of arrays containing the residue
          positions associated to each sector as published in the original
          paper.
        - `pdb_sequence`: sequence extracted from E. coli's PDB structure
        - `pdb_positions`: positions extracted from E. coli's PDB structure

    """

    module_path = os.path.dirname(__file__)
    filename = os.path.join(
        module_path,
        "data/rhomboid_proteases/Data_S1_Rhomboid_MSA.fasta")
    data = load_MSA(filename, format="fasta")

    filename = os.path.join(
        module_path,
        "data/rhomboid_proteases/rhomboid_sectors.npz")
    sectors = np.load(filename)

    # Load the metadata
    filename = os.path.join(
        module_path,
        "data/rhomboid_proteases/rhomboid_Uniprot_metadata.tsv")
    metadata = pd.read_csv(filename, sep="\t")

    # Load the PDB structure
    filename = os.path.join(
        module_path,
        "data/rhomboid_proteases/2NRF.pdb")
    # Two chains: A or B
    pdb_sequence, pdb_positions = load_pdb(filename, '2NRF', 'A')

    data["sector_positions"] = sectors
    data["metadata"] = metadata
    data["pdb_sequence"] = pdb_sequence,
    data["pdb_positions"] = pdb_positions

    return data


def load_DHFR():
    """
    load the DHFR dataset

    This dataset comes from Kalmer et al, The Journal of Physical Chemistry B,
    2024 (https://pubs.acs.org/doi/10.1021/acs.jpcb.4c04195)

    Returns
    -------
    a dictionnary containing :
        - `sequence_ids`: a list of strings corresponding to sequence names
        - `alignment`: a list of strings corresponding to sequences. Because it
          is an MSA, all the strings are of same length.
         - `sector_positions`: a dictionnary of arrays containing the residue
          positions associated to each sector as published in the original
          paper.
        - `pdb_sequence`: sequence extracted from E. coli's PDB structure
        - `pdb_positions`: positions extracted from E. coli's PDB structure
    """
    module_path = os.path.dirname(__file__)

    filename = os.path.join(
        module_path,
        "data/DHFR/alignment.faa.gz")
    with gzip.open(filename, "rt") as f:
        data = load_MSA(f, format="fasta")

    filename = os.path.join(
        module_path,
        "data/DHFR/DHFR_sectors.npz")
    sectors = np.load(filename)

    # Load the PDB structure
    filename = os.path.join(
        module_path,
        "data/DHFR/3QL0.pdb")
    pdb_sequence, pdb_positions = load_pdb(filename, '3QL0', 'A')

    data["sector_positions"] = sectors
    data["pdb_sequence"] = pdb_sequence,
    data["pdb_positions"] = pdb_positions

    return data
