from . import msa
from . import statistics
from . import decomposition
from . import __params

import pandas as pd
import numpy as np


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

    coevolution_metric : str or callable, optional, default: 'SCA'
        which coevolution metric to use:

        - SCA: the coevolution matrix from Rivoire et al
        - MI: the mutual information
        - NMI: the normalized mutual information
        - callable: a function that takes as arguments (1) sequences, (2)
          `seq_weights`, and `freq_regul`

    correction : {None, 'APC', 'entropy'}, default: None
        which correction to use

    Returns
    -------
    coevol_matrix : np.ndarray (n_filtered_pos, n_filtered_pos)
        coevolution matrix

    coevol_matrix_ngm : np.ndarray (n_filtered_pos, n_filtered_pos)
        coevolution matrix without global mode (ngm = no global mode)

    df : pd.DataFrame with the following columns

        - original_msa_pos : the original MSA position
        - filtered_msa_pos : the position in the filtered MSA

        and  for each component:

        - PCk: the projection of the residue onto the kth principal component
        - ICk: the projeciton of the residue onto the kth independent
          component
        - xcor_k: wherether the residue is found to be part of xcor k

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
    elif callable(coevolution_metric):
        coevol_matrix = coevolution_metric(
            seq_kept, seq_weights=seq_weights,
            freq_regul=freq_regul)
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
                seq_kept,
                seq_weights=seq_weights)
            coevol_matrix = statistics.pairwise.compute_entropy_correction(
                coevol_matrix, entropy_aa)
        else:
            raise ValueError(
                "Unknown 'correction' value. User provided"
                f"{correction}. Options are 'APC', 'entropy'")

    # Now, compute deconvolution

    principal_components = decomposition.extract_principal_components(
        coevol_matrix)
    independent_components = decomposition.extract_independent_components(
        coevol_matrix, n_components=n_components)
    xcors = decomposition.extract_xcors_from_ICs(
        independent_components, coevol_matrix)

    # Now, map everything into a nice pandas DataFrame
    pos_mapping, _ = msa.map_msa_positions(len(sequences[0]), pos_kept)

    df = pd.DataFrame(
        {"original_msa_pos": np.arange(len(sequences[0]), dtype=int),
         "filtered_msa_pos": pos_mapping.values()})
    # make filtered_msa_pos stay integer with NaN support
    df["filtered_msa_pos"] = df["filtered_msa_pos"].astype("Int64")

    # Add PCA and ICA results
    for k in range(n_components):
        df.loc[~df["filtered_msa_pos"].isna(),
               "PC%d" % (k+1)] = principal_components[k]
        df.loc[~df["filtered_msa_pos"].isna(),
               "IC%d" % (k+1)] = independent_components[k]
        df["xcor_%d" % (k+1)] = np.isin(
            df["filtered_msa_pos"], xcors[k])
        df.loc[~df["filtered_msa_pos"].isna(),
               "xcor_%d" % (k+1)] = np.isin(
                   df.loc[~df["filtered_msa_pos"].isna(),
                          "filtered_msa_pos"], xcors[k])

    coevol_matrix_ngm = decomposition.remove_global_correlations(coevol_matrix)

    return coevol_matrix, coevol_matrix_ngm, df
