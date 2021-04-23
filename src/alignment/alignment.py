from typing import Tuple, List
# from numpy.typing import ArrayLike
import numpy as np

class Alignment:
    """
    Represents a Multiple Sequence Alignment

    Attributes
    ----------


    """
    # Used to translate AAs to Ints to make use of Numpy vectorization
    __code_tonum = {
        '-': 0,
        'X': 0,  # Consider unknown AAs as gaps
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
        'Y': 20
    }

    _kws = {
        'text': '__text_rep',
        'numerical': '__num_rep',
        'num': '__num_rep'
    }

    _freq0 = np.array(
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
            0.033,
        ]
    )

    _lbda = 0.03

    def __init__(self, alignment_sequence_list: List[str]) -> None:
        """

        Parameters
        ----------
        alignment_sequence_list: list of str
            list of strings corresponding to one string per sequence
        """
        self.__raw = alignment_sequence_list
        # Split characters so that one element in each row represents one aminoacid
        self.__text_rep = np.array(
            [np.array([char for char in row]) for row in alignment_sequence_list]
        )
        self.__num_rep = np.zeros(self.__text_rep.shape, dtype=np.int64)
        for k in Alignment.__code_tonum.keys():
            self.__num_rep[self.__text_rep == k] = self.__code_tonum[k]

    def gap_frequency(self, threshold=0.2) -> Tuple[np.ndarray, float]:
        """Computes the gap frequency for a MSA

        Parameters
        ----------
        threshold: float
            Positions will be considered overly-gapped if their gap ratio is
            above this value.

        Returns
        -------
        gap_frequency: numpy 1D array
            An array of length n_pos with each value set to the gap frequency at
            the position corresponding to the index.
        background_gap_frequency: float
            Containing the gap frequency accross all positions which are below
            the desired threshold.
        """
        gap_frequency = []

        # per position
        for residue in range(self.__text_rep.shape[1]):
            currt_freq = np.count_nonzero(
                self.__text_rep[:, residue] == '-'
            ) / self.__text_rep.shape[0]

            gap_frequency.append(currt_freq)

        gap_frequency = np.array(gap_frequency)

        working_pos = gap_frequency < threshold
        working_pos = working_pos.reshape(working_pos.shape[0], -1)

        mask = np.ones((self.__text_rep.shape[0], 1), dtype="bool").dot(working_pos.T)
        bg_gap_frequency = self.__text_rep[mask]

        return gap_frequency, float(np.mean(bg_gap_frequency == '-'))

    def filtered_alignment(self, threshold=0.2) -> np.ndarray:
        """Returns the filtered alignment

        Parameters
        ----------
        threshold: float or None
            The threshold to discriminate positions to keep versus overly-gapped
            positions.

        Returns
        -------
        alignment: numpy array of int
            The filtered alignment with overly-gapped positions removed.
        """
        gap_frequency, _ = self.gap_frequency(threshold)
        return self.__num_rep[:, gap_frequency < threshold].copy()

    def similarity(self, use_filtered=True, threshold=0.2) -> np.ndarray:
        """Computes the similarity between sequences in a MSA

        Parameters
        ----------
        use_filtered: bool or None
            Make use of the filtered MSA, according to the gap frequency
            threshold (defaults to True).
        threshold: float or None
            The threshold to filter overly-gapped positions, unused if
            use_filtered is False (default 0.2).

        Returns
        -------
        similarity_matrix: numpy 2D matrix of shape (n_seq,n_seq)
            The similarity matrix between sequences: a symmetric square matrix.
        """
        src = self.filtered_alignment(threshold) if use_filtered else self.__num_rep
        n_seqs = src.shape[0]
        res = np.zeros(shape=(n_seqs, n_seqs))
        for seq_1_i in range(n_seqs):
            seq_1 = src[seq_1_i]
            for seq_2_i in range(seq_1_i, n_seqs):
                seq_2 = src[seq_2_i]
                dist = np.mean(np.array([seq_1[i] == seq_2[i] for i in range(len(seq_1))]))
                res[seq_1_i, seq_2_i] = dist
                res[seq_2_i, seq_1_i] = dist
        return res

    @staticmethod
    def aa_freq_at_pos(
        algnt_col: np.ndarray, lbda=0.03, aa_count=21, weights=None
    ) -> np.ndarray:
        """Computes frequencies of aminoacids at a specific position

        Parameters
        ----------
        algnt_col: numpy series
            The column for which frequency should be computed.
        lbda: float or None
            A regularization parameter. Defaults to 0.03
        aa_count: int or None
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category.
        weights: numpy 1D array or None
            Gives more or less importance to certain sequences.
            If provided, the length of this array should be equal to the length
            of algnt_col

        Returns
        -------
        freqs: numpy 1D array
            A numpy array of length aa_count, containing regularized frequencies
            for each amino acid. A mapping from amino-acid to indices in this
            result can be found in __code_tonum.
        """
        # AA proportion
        bins = np.bincount(algnt_col, weights=weights)

        # Pad with zeroes for AAs not present at currt position
        freqs = np.pad(bins, (0, aa_count - bins.shape[0]))
        freqs = np.array(freqs)

        # Frequency
        if weights is None:
            weights = np.ones(algnt_col.shape)

        freqs = freqs / weights.sum()

        # Regularization
        freqs = (1 - lbda) * freqs + lbda / aa_count
        return freqs

    def sca_entropy(self, aa_count=21, weights=None, threshold=0.2) -> np.ndarray:
        """Computes entropy as per SCA guidelines

        Parameters
        ----------
        aa_count: int
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category.
        weights: numpy 1D array
            Gives more or less importance to certain sequences in the MSA.
        threshold: float
            The gap frequency low-pass threshold.

        Returns
        -------
        entropy: numpy 1D array
            The entropy for each of the positions in the filtered alignment.
        """
        _, gap_bg_freq = self.gap_frequency(threshold)
        working_alignment = self.filtered_alignment(threshold)

        mult_fact = 1 - gap_bg_freq
        freq0 = Alignment._freq0 * mult_fact
        freq0g = np.insert(freq0, 0, gap_bg_freq)

        # Prevent MSA without gaps & regularize background frequencies
        # the same way data is regularized
        freq0g_n = (1 - Alignment._lbda) * freq0g + Alignment._lbda / aa_count

        freqs = np.apply_along_axis(
            Alignment.aa_freq_at_pos,
            0,
            working_alignment,
            weights=weights
        )
        rel_entropy = np.sum(freqs * np.log(freqs / freq0g_n[:,np.newaxis]), 0)

        rel_entropy_derivative = np.apply_along_axis(
            lambda x: abs(np.sum(np.log((x * (1 - freq0g_n)) / ((1 - x) * freq0g_n)))),
            0,
            freqs
        )
        return rel_entropy, rel_entropy_derivative

    def similarity_weights(self, similarity: np.ndarray, threshold=0.2) -> float:
        """Get the effective number of sequences after applying weights

        Parameters
        ----------
        similarity: numpy 2D array
            A similarity matrix of shape (n_seq,n_seq).
        threshold: float or None
            The "width" (delta) of an effective sequence: defines the threshold
            at which sequences are regrouped as an effective sequence.
            Defaults to 0.2.

        Returns
        -------
        weight: float
            Weights
        """
        return (1 / np.sum(similarity >= threshold, axis=0))

    def weight_bins(
        self, threshold_increment=.001, low_bin_count=100, weight_vector_count=5,
        use_filtered=True, threshold=0.2
    ) -> Tuple[List[float], List[List[float]]]:
        """Computes weight bins using a similarity matrix, in order to regroup
           "similar-enough" sequences in bins

        Parameters
        ----------
        threshold_increment: float or None
            The rate at which delta should be decremented when trying to build
            bins, defaults to 0.001.
        low_bin_count: int or None
            The target effective sequence count for the smallest bin.
        weight_vector_count: int or None
            The number of bins to create. Defaults to 5.
        use_filtered: bool or None
            Use a filtered version of the alignment according to gap cutoffs.
            Defaults to True.
        threshold: float or None
            The threshold over which positions are considered overly-gapped.

        Returns
        -------
        delta_list: list of float
            A list of the used deltas: the width of effective sequences groups.
        weight_list: list of list of float
            A list of list of weights: for each position in this list, there are
            n_seq weights, grouped using the corresponding delta found in
            delta_list.
        """
        similarity_matrix = self.similarity(use_filtered=use_filtered, threshold=threshold)
        bin_sizes = np.logspace(
            np.log10(low_bin_count),
            np.log10(similarity_matrix.shape[0]),
            weight_vector_count
        ).tolist()

        delta_list = []
        weight_list = []
        currt_delta = np.max(similarity_matrix)
        fnd_values = len(delta_list)
        while fnd_values != weight_vector_count:
            sim_wt = self.similarity_weights(similarity_matrix, currt_delta)
            eff_seqs = np.sum(sim_wt)
            # If the tree is very well structured, there is no point in trying
            # to find values very close to the generated bin sizes
            if eff_seqs < bin_sizes[-1 - fnd_values]:
                delta_list.append(currt_delta)
                weight_list.append(sim_wt)
                fnd_values += 1
            currt_delta -= threshold_increment
        return delta_list, weight_list

    def get_sequence_count(self) -> int:
        """Computes the number of sequences in the current instance

        Returns
        -------
        n_seq: int
            Number of sequences
        """
        return self.__num_rep.shape[0]

    def most_frequent_aa(self, filtered=True, threshold=.2) -> Tuple[np.ndarray, np.ndarray]:
        """Computes the most frequent amino acid at each position

        Parameters
        ----------
        filtered: bool or None
            Whether overly-gapped positions should be excluded. Defaults to
            True.
        threshold: float or None
            The threshold over which positions are considered overly-gapped.
            Useless if filtered is False. Defaults to 0.2.
        """
        most_freq_aa = []
        most_freq_aa_freq = []
        working_alignment = (
            self.filtered_alignment(threshold=threshold) if filtered else self.__num_rep
        )
        for i in range(working_alignment.shape[1]):
            currt_col = working_alignment[:, i]
            currt_bins = np.bincount(currt_col)
            currt_max_i = np.argmax(currt_bins)
            most_freq_aa.append(currt_max_i)
            most_freq_aa_freq.append(currt_bins[currt_max_i] / np.sum(currt_bins))
        return np.array(most_freq_aa), np.array(most_freq_aa_freq)

    def _frequency(self, algnt_col1, algnt_col2, aa1, aa2, weights) -> Tuple[float, float, float]:
        """FOR LEARNING PURPOSES ONLY: prefer the optimized version
        Computes frequencies for amino acids to build a coevolution matrix
        """
        a_i = algnt_col1 == aa1
        b_j = algnt_col2 == aa2
        fia = np.sum(a_i * weights) / np.sum(weights)
        fjb = np.sum(b_j * weights) / np.sum(weights)
        fijab = np.sum((a_i & b_j) * weights) / np.sum(weights)
        return fia, fjb, fijab

    def _slow_coevolution_matrix(
        self, weights, threshold=.2, aa_count=21
    ) -> Tuple[np.ndarray, np.ndarray]:
        """FOR LEARNING PURPOSES ONLY: prefer the optimized version
        """
        algnt = self.filtered_alignment(threshold=threshold)
        seq_size = algnt.shape[1]
        Cijab = np.zeros((seq_size, seq_size, aa_count, aa_count))
        Cij = np.zeros((seq_size, seq_size))
        for i in range(seq_size):
            for j in range(i, seq_size):
                for a in range(aa_count):
                    for b in range(aa_count):
                        fia, fjb, fijab = self._frequency(algnt[:, i], algnt[:, j], a, b, weights)
                        Cijab[i, j, a, b] = fijab - fia * fjb
                val = np.sqrt(np.sum(Cijab[i, j, :, :] ** 2))
                Cij[i, j] = val
                Cij[j, i] = val
        return Cijab, Cij

    def aminoacid_joint_frequencies(
        self, weights: np.ndarray, threshold=.2, aa_count=21, use_tensorflow=False
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Computes the joint frequencies for every amino acid at every position

        Parameters
        ----------
        weights: numpy 1D array of float
            Gives more or less importance to certain sequences in the MSA.
        threshold: float or None
            The threshold over which positions are considered overly-gapped and
            are thus filtered out.
        aa_count: int or None
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category.
        use_tensor: bool or None
            Whether to use tensordot to compute joint frequencies

        Returns
        -------
        joint_freqs: numpy 4D matrix of shape (n_pos, n_pos, aa_count, aa_count)
            For joint_freqs[i,j,a,b], the joint frequencies for of amino acids
            a and b at positions i and j respectively.

        joint_freqs_ind: numpy 4D matrix of shape (n_pos, n_pos, aa_count, aa_count)
            For joint_freqs[i,j,a,b], the joint frequencies for of amino acids
            a and b at positions i and j respectively, if a at position i and j
            at position b are independent random variables.
        """
        algnt = self.filtered_alignment(threshold=threshold)
        seq_size = algnt.shape[1]

        # Joint frequencies
        joint_freqs = np.zeros((seq_size, seq_size, aa_count, aa_count))

        # Frequencies IF AAs are present independently as positions i & j
        joint_freqs_ind = np.zeros((seq_size, seq_size, aa_count, aa_count))

        # One-hot encoding (binary array)
        binary_array = np.array([algnt == aa for aa in range(aa_count)])

        # Adding weights
        weighted_binary_array = binary_array * weights[np.newaxis, :, np.newaxis]

        m_eff = np.sum(weights)
        simple_freq = weighted_binary_array / m_eff
        simple_freq = np.sum(simple_freq, axis=1)

        if not use_tensorflow:
            joint_freq_aibj = np.tensordot(
                weighted_binary_array, binary_array, axes=([1], [1])
            )/m_eff
            joint_freqs = joint_freq_aibj.transpose(1, 3, 0, 2)
        else:
            import tensorflow as tf
            joint_freq_aibj = tf.tensordot(
                weighted_binary_array, binary_array * 1.0, axes=[[1], [1]]
            ) / m_eff
            joint_freqs = tf.transpose(joint_freq_aibj, perm=[1, 3, 0, 2])

        joint_freqs_ind = np.multiply.outer(simple_freq, simple_freq)

        joint_freqs_ind = np.moveaxis(joint_freqs_ind, [0, 1, 2, 3], [2, 0, 3, 1])

        return joint_freqs, joint_freqs_ind

    def coevolution_matrix(
        self, weights: np.ndarray, threshold=.2, aa_count=21, use_tensorflow=False
    ) -> np.ndarray:
        """Computes the coevolution matrix (the SCA-way)

        Parameters
        ----------
        weights : np.ndarray
            Sequence weights
        threshold : float, optional
            Gap-frequency low-pass value, by default .2
        aa_count : int, optional
            Number of aminoacids to consider, by default 21
        use_tensorflow : bool, optional
            Whether tensorflow should be used to compute frequency matrices, by default False

        Returns
        -------
        np.ndarray
            The weighted coevolution matrix of shape (N_pos, N_pos)
        """
        jf, jfi = self.aminoacid_joint_frequencies(
            weights, threshold, aa_count, use_tensorflow
        )
        _, w = self.sca_entropy(aa_count, weights, threshold)
        cijab = jf - jfi
        cijab_w = (cijab
                   * w[np.newaxis, :, np.newaxis, np.newaxis]
                   * w[:, np.newaxis, np.newaxis, np.newaxis])
        return np.sqrt(np.sum(cijab_w ** 2, axis=(2, 3)))

