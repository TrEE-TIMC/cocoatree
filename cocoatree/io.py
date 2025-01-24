from Bio import AlignIO
from Bio.PDB import PDBParser
from .msa import _clean_msa
from .__params import aatable
from ete3 import Tree
import numpy as np


def load_MSA(file_path, format, clean=True, verbose=False):
    """Read in a multiple sequence alignment (MSA)

    Arguments
    ---------
    file_path : path to the alignment file

    format : format of the alignment file (e.g. 'fasta', 'phylip', etc.)

    verbose : boolean,
            whether to print informations about the MSA

    Returns
    -------
    seq_id : list of sequence identifiers

    sequences : list of sequences as strings
    """

    alignment = AlignIO.read(file_path, format)

    if clean:
        alignment = _clean_msa(alignment)

    seq_id = list()
    sequences = list()
    for record in alignment:
        seq_id.append(record.id)
        sequences.append(str(record.seq))

    if verbose:
        print('Number of sequences: %i' % len(alignment))
        print('Alignment of length: %i' % len(alignment[0]))

    return seq_id, sequences


def load_tree_ete3(file_path):
    """
    From the loading of a Newick tree, generate a ete3.Tree object

    Arguments
    ---------
    file_path : path to the Newick file

    Returns
    -------
    tree_ete3 : `ete3.Tree` object

    """
    tree_ete3 = Tree(file_path, format=0)
    return tree_ete3


def export_fasta(sequences, seq_id, outpath):
    """
    Function to export intermediate files in fasta format

    Arguments
    ---------
    sequences : list of sequences as strings (as imported by load_MSA)

    seq_id : list of sequences identifiers (as imported by load_MSA)

    outpath : path to the output file
    """

    # Add checks to see if the path exists?
    Nseq = len(sequences)
    with open(outpath, 'w') as outfile:
        for record in range(0, Nseq):
            outfile.write('>' + str(seq_id[record]) + '\n')
            outfile.write(str(sequences[record]) + '\n')


def load_pdb(path2pdb, pdb_id, chain):

    '''
    Import a PDB file and extract the associated sequence along with the
    amino acid positions

    Arguments
    ---------
    path2pdb : path to the PDB file

    pdb_id : str,
        identifier of the PDB file

    chain : str,
        name of the chain to read

    Returns
    -------
    pbd_seq : str,
        amino acid sequence of the PDB file

    pdb_pos : list,
        PDB position of each amino acid
    '''

    P = PDBParser(PERMISSIVE=1)
    structure = P.get_structure(pdb_id, path2pdb)

    # Fill up sequence and label information
    pdb_seq = ""
    pdb_pos = list()
    residues = [res for res in structure[0][chain] if res.get_id()[0] == " "]
    for res in residues:
        pdb_pos.append(str(res.get_id()[1]) + str(res.get_id()[2]).strip())
        try:
            pdb_seq += aatable[res.get_resname()]
        except BaseException as e:
            print("Error: " + str(e))
            pdb_seq += "X"

    return pdb_seq, pdb_pos


def export_sector_for_pymol(mapping, independent_components, axis, sector_pos,
                            ics, outpath):
    """
    Export numpy arrays of a sector's residue positions and their contribution
    for coloring in PyMol.

    Arguments
    ---------
    mapping : numpy.ndarray,
        mapping between the unfiltered MSA and the PDB structure, output of
        cocoatree.msa.map_to_pdb() function

    independent_components : numpy.ndarray,
        output of cocoatree.deconvolution.compute_ica() function

    axis : int,
        rank of the independent component associated with the desired sector

    sector_pos : list,
        positions of the sector's residues in the unfiltered MSA

    ics : numpy.ndarray,
        positions of the sector's residues in the filtered MSA

    outpath : str,
        path to the output file as a binary in .npy format

    Returns
    -------
    binary file in .npy format containing an array with the positions of the
    sector's residues and an array with their contribution to the independent
    component.
    """

    sector_pdb_pos = []
    for residue in sector_pos:
        index = np.where(mapping[2] == residue)[0][0]
        sector_pdb_pos.append(mapping[1][index])

    ic_contributions = []
    for residue in ics[0].items:
        ic_contributions.append(independent_components[residue, axis])

    np.save(outpath, np.array([sector_pdb_pos, ic_contributions]))
