from . import msa
from . import datasets
from . import statistics
from . import io
from . import deconvolution
from . import __params
import pandas as pd
import numpy as np


__version__ = "0.0.0a0.dev0"


def perform_sca(sequences_id, sequences,
                n_components=4,
                freq_regul=__params.__freq_regularization_ref,
                gap_threshold=0.4, seq_threshold=0.2,
                coevolution_metric="SCA", correction=None):
    """
    Perform statistical coupling analysis (SCA)

    Parameters
    ----------
    sequences : list of MSA sequences to filter

    sequences_id : list of the MSA's sequence identifiers

    n_components : int, default: 4

    gap_threshold : float [0, 1], default: 0.4
        max proportion of gaps tolerated

    seq_threshold : maximum fraction of gaps per sequence (default 0.2)

    coevolution_metric : {'SCA', 'NMI', 'MI}, optional, default: 'SCA'
        which coevolution metric to use:

        - SCA: the coevolution matrix from Rivoire et al
        - MI: the mutual information
        - NMI: the normalized mutual information

    correction : {None, 'APC', 'entropy'}, default: None
        which correction to use

    Returns
    -------
    results : pd.DataFrame with the following columns

        - original_msa_pos : the original MSA position
        - filtered_msa_pos : the position in the filtered MSA

        and  for each component:

        - PCk: the projection of the residue onto the kth principal component
        - ICk: the projeciton of the residue onto the kth independent
          component
        - sector_k: wherether the residue is found to be part of sector k

    """

    # Start by filtering sequences
    seq_kept, seq_kept_id, pos_kept = msa.filter_sequences(
        sequences, sequences_id, gap_threshold=gap_threshold,
        seq_threshold=seq_threshold)

    # Compute sequence weights. This is mostly to avoid recomputing it at
    # several step in the pipeline and thus speed things up a bit
    seq_weights, _ = msa.compute_seq_weights(seq_kept)

    # Compute co-evolution matrix
    if coevolution_metric == "SCA":
        coevol_matrix = statistics.pairwise.compute_sca_matrix(
            seq_kept,
            seq_weights=seq_weights,
            freq_regul=freq_regul)
    elif coevolution_metric == "MI":
        coevol_matrix = statistics.pairwise.compute_mutual_information_matrix(
            seq_kept, seq_weights=seq_weights, freq_regul=freq_regul,
            normalize=False)
    elif coevolution_metric == "NMI":
        coevol_matrix = statistics.pairwise.compute_mutual_information_matrix(
            seq_kept, seq_weights=seq_weights, freq_regul=freq_regul)
    else:
        raise ValueError(
            "Unknown 'coevol_metric' value. User provided"
            f"{coevolution_metric}. Options are 'SCA', 'MI', 'NMI'")

    # Compute correction on coevolution matrix
    if correction is not None:
        if correction == "APC":
            _, coevol_matrix = statistics.pairwise.compute_apc(coevol_matrix)
        elif correction == "entropy":
            entropy_aa = statistics.position.compute_conservation(
                sequences,
                seq_weights=seq_weights)
            coevol_matrix = statistics.pairwise.compute_entropy_correction(
                coevol_matrix, entropy_aa)
        else:
            raise ValueError(
                "Unknown 'correction' value. User provided"
                f"{correction}. Options are 'APC', 'entropy'")

    # Now, compute deconvolution

    principal_components = deconvolution.extract_principal_components(
        coevol_matrix)
    independent_components = deconvolution.extract_independent_components(
        coevol_matrix, n_components=n_components)
    sectors = deconvolution.extract_sectors(
        independent_components, coevol_matrix)

    # Now, map everything into a nice pandas DataFrame
    pos_mapping, _ = msa.map_msa_positions(len(sequences[0]), pos_kept)

    results = pd.DataFrame(
        {"original_msa_pos": np.arange(len(sequences[0]), dtype=int),
         "filtered_msa_pos": pos_mapping.values()})

    # Add PCA and ICA results
    for k in range(n_components):
        results.loc[~results["filtered_msa_pos"].isna(),
                    "PC%d" % (k+1)] = principal_components[k]
        results.loc[~results["filtered_msa_pos"].isna(),
                    "IC%d" % (k+1)] = independent_components[k]

        results["sector_%d" % (k+1)] = np.isin(
            results["filtered_msa_pos"],
            sectors[k])
        results.loc[~results["filtered_msa_pos"].isna(),
                    "sector_%d" % (k+1)] = np.isin(
                        results.loc[~results["filtered_msa_pos"].isna(),
                                    "filtered_msa_pos"],
                        sectors[k])

    return results
