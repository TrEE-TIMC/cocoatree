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

    # Values taken from pySCA
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
            0.033,
        ]
    )

    __lbda = 0.03

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
        threshold: float, optional
            Positions will be considered overly-gapped if their gap ratio is
            above this value, by default 0.2

        Returns
        -------
        np.ndarray of shape (n_pos)
            An array of length n_pos with each value set to the gap frequency at
            the position corresponding to the index.
        float
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
        threshold: float, optional
            The threshold to discriminate positions to keep versus overly-gapped
            positions, by default 0.2

        Returns
        -------
        np.ndarray of int
            The filtered alignment with overly-gapped positions removed.
        """
        gap_frequency, _ = self.gap_frequency(threshold)
        return self.__num_rep[:, gap_frequency < threshold].copy()

    def seq_len(self):
        """Retrieves the length of the MSA

        Returns
        -------
        int
            Length of the MSA
        """
        return self.__num_rep.shape[1]

    def seq_count(self):
        """Retrieves the number of sequences in the MSA

        Returns
        -------
        int
            Number of sequences in the MSA
        """
        return self.__num_rep.shape[0]

    def filtered_seq_len(self, threshold):
        """Retrieves the length of the filtered MSA

        Parameters
        ----------
        threshold : float
            The threshold to discriminate positions to keep versus overly-gapped
            positions.

        Returns
        -------
        int
            Length of the filtered MSA.
        """
        return self.filtered_alignment(threshold).shape[1]

    def similarity(self, use_filtered=True, threshold=0.2) -> np.ndarray:
        """Computes the similarity between sequences in a MSA

        Parameters
        ----------
        use_filtered: bool, optional
            Make use of the filtered MSA, according to the gap frequency
            threshold, by default True
        threshold: float, optional
            The threshold to filter overly-gapped positions, unused if
            use_filtered is False, by default 0.2

        Returns
        -------
        np.ndarray of shape (n_seq,n_seq)
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
        algnt_col: np.ndarray, lbda=0.03, aa_count=21, weights=None, regularization=True
    ) -> np.ndarray:
        """Computes frequencies of aminoacids at a specific position

        Parameters
        ----------
        algnt_col: np.series
            The column for which frequency should be computed.
        lbda: float, optional
            A regularization parameter, by default 0.03
        aa_count: int, optional
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category.
        weights: numpy 1D array, optional
            Gives more or less importance to certain sequences.
            If provided, the length of this array should be equal to the length
            of algnt_col, by default None

        Returns
        -------
        np.ndarray of shape (aa_count)
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

        if regularization:
            freqs = (1 - lbda) * freqs + lbda / aa_count
        return freqs

    def regularized_background_frequencies(self, threshold=0.2, aa_count=21):
        """Computes background frequencies for every amino acid in the current alignment

        Parameters
        ----------
        threshold : float, optional
            Low-pass gap frequency filter, by default 0.2
        aa_count : int, optional
            Number of amino acids to consider, by default 21

        Returns
        -------
        np.ndarray of shape (aa_count)
            The regularized frequencies for every amino acid.
        """
        _, gap_bg_freq = self.gap_frequency(threshold)
        mult_fact = 1 - gap_bg_freq
        freq0 = Alignment.__freq0 * mult_fact
        freq0g = np.insert(freq0, 0, gap_bg_freq)

        # Prevent MSA without gaps & regularize background frequencies
        # the same way data is regularized
        freq0g_n = (1 - Alignment.__lbda) * freq0g + Alignment.__lbda / aa_count
        return freq0g_n

    def sca_entropy(self, aa_count=21, weights=None, threshold=0.2) -> np.ndarray:
        """Computes entropy as per SCA guidelines

        Parameters
        ----------
        aa_count: int
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category.
        weights: np.ndarray of shape (N_seq)
            Gives more or less importance to certain sequences in the MSA.
        threshold: float
            The gap frequency low-pass threshold.

        Returns
        -------
        np.ndarray of shape (N_pos, N_aa)
            The entropy for each of the positions in the filtered alignment.
        """
        working_alignment = self.filtered_alignment(threshold)

        freq0g_n = self.regularized_background_frequencies(threshold)

        freqs = np.apply_along_axis(
            Alignment.aa_freq_at_pos,
            0,
            working_alignment,
            weights=weights
        )
        rel_entropy = (
            freqs * np.log(freqs / freq0g_n[:, np.newaxis])
            + (1 - freqs) * np.log((1 - freqs)/(1 - freq0g_n[:, np.newaxis]))
        )
        return rel_entropy.transpose([1, 0])

    def position_weights_sca(self, threshold=0.2, weights=None):
        """Computes SCA weights positions

        Parameters
        ----------
        threshold : float, optional
            The gap frequency low-pass threshold, by default 0.2
        weights : np.ndarray, optional
            Gives more or less importance to certain sequences in the MSA, by default None

        Returns
        -------
        np.ndarray of shape (N_pos, N_pos, aa_count, aa_count)
            Weights as computed in SCA (using the derivative of the relative entropy)
        """
        working_alignment = self.filtered_alignment(threshold)

        freq0g_n = self.regularized_background_frequencies(threshold)

        freqs = np.apply_along_axis(
            Alignment.aa_freq_at_pos,
            0,
            working_alignment,
            weights=weights
        )

        rel_entropy_derivative = (
            np.log(freqs * (1 - freq0g_n[:, np.newaxis]) / (1 - freqs) * freq0g_n[:, np.newaxis])
            .transpose([1, 0])
        )

        return np.multiply.outer(
            rel_entropy_derivative,
            rel_entropy_derivative
        ).transpose([0, 2, 1, 3])

    def similarity_weights(self, similarity: np.ndarray, threshold=0.2) -> float:
        """Get the effective number of sequences after applying weights

        Parameters
        ----------
        similarity: numpy 2D array of shape (N_seq, N_seq)
            A similarity matrix
        threshold: float, optional
            The "width" (delta) of an effective sequence: defines the threshold
            at which sequences are regrouped as an effective sequence.
            Defaults to 0.2.

        Returns
        -------
        float
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
        threshold_increment: float, optional
            The rate at which delta should be decremented when trying to build
            bins, defaults to 0.001.
        low_bin_count: int, optional
            The target effective sequence count for the smallest bin.
        weight_vector_count: int, optional
            The number of bins to create. Defaults to 5.
        use_filtered: bool, optional
            Use a filtered version of the alignment according to gap cutoffs.
            Defaults to True.
        threshold: float, optional
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

    def most_frequent_aa(self, filtered=True, threshold=.2) -> Tuple[np.ndarray, np.ndarray]:
        """Computes the most frequent amino acid at each position

        Parameters
        ----------
        filtered: bool, optional
            Whether overly-gapped positions should be excluded, by default True
        threshold: float, optional
            The threshold over which positions are considered overly-gapped.
            Useless if filtered is False, by default 0.2
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
        weights: np.ndarray of float
            Gives more or less importance to certain sequences in the MSA.
        threshold: float, optional
            The threshold over which positions are considered overly-gapped and
            are thus filtered out, by default 0.2
        aa_count: int, optional
            The number of amino acids to consider for the shape of the return
            value. Defaults to 21 to account for unknown (X) amino acids which
            are considered as gaps for the lack of a better category, by default 21
        use_tensor: bool, optional
            Whether to use tensordot to compute joint frequencies, by default False

        Returns
        -------
        np.ndarray of shape (n_pos, n_pos, aa_count, aa_count)
            For joint_freqs[i,j,a,b], the joint frequencies for of amino acids
            a and b at positions i and j respectively.

        np.ndarray of shape (n_pos, n_pos, aa_count, aa_count)
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

        simple_freq = (1 - self.__lbda) * simple_freq + self.__lbda / aa_count

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

        joint_freqs = (1 - self.__lbda) * joint_freqs + self.__lbda / aa_count

        joint_freqs_ind = np.multiply.outer(simple_freq, simple_freq)

        joint_freqs_ind = np.moveaxis(joint_freqs_ind, [0, 1, 2, 3], [2, 0, 3, 1])

        return joint_freqs, joint_freqs_ind

    def sca_coevolution(
        self, seq_weights=None, pos_weights=None, threshold=.2, aa_count=21, use_tensorflow=False
    ) -> np.ndarray:
        """Computes the SCA coevolution matrix

        Parameters
        ----------
        seq_weights : np.ndarray of float of shape (N_seq)
            Sequence weights
        pos_weights : np.ndarray of float of shape (N_pos, N_pos, aa_count, aa_count)
            Position weights
        threshold : float, optional
            Gap-frequency low-pass value, by default .2
        aa_count : int, optional
            Number of aminoacids to consider, by default 21
        use_tensorflow : bool, optional
            Whether tensorflow should be used to compute frequency matrices, by default False

        Returns
        -------
        np.ndarray
            The weighted SCA coevolution matrix of shape (N_pos, N_pos)
        """
        seq_weights = seq_weights
        if seq_weights is None:
            seq_weights = np.ones(self.seq_count())
        else:
            assert seq_weights.shape[0] == self.seq_count()

        pos_weights = pos_weights
        if pos_weights is None:
            filtered_seq_len = self.filtered_seq_len(threshold)
            pos_weights = np.ones((filtered_seq_len, filtered_seq_len, aa_count, aa_count))

        jf, jfi = self.aminoacid_joint_frequencies(
            seq_weights, threshold, aa_count, use_tensorflow
        )

        cijab = jf - jfi
        cijab_w = cijab * pos_weights
        return np.sqrt(np.sum(cijab_w ** 2, axis=(2, 3)))

    def mutual_information(
        self, seq_weights=None, pos_weights=None, threshold=.2, aa_count=21, use_tensorflow=False
    ) -> np.ndarray:
        """Computes a mutual-information coevolution matrix

        Parameters
        ----------
        seq_weights : np.ndarray
            Sequence weights
        pos_weights : np.ndarray of float of shape (N_pos, N_pos, aa_count, aa_count)
            Position weights
        threshold : float, optional
            Gap-frequency low-pass value, by default .2
        aa_count : int, optional
            Number of aminoacids to consider, by default 21
        use_tensorflow : bool, optional
            Whether tensorflow should be used to compute frequency matrices, by default False

        Returns
        -------
        np.ndarray
            The weighted mutual-information coevolution matrix of shape (N_pos, N_pos)
        """

        seq_weights = seq_weights
        if seq_weights is None:
            seq_weights = np.ones(self.seq_count())
        else:
            assert seq_weights.shape[0] == self.seq_count()

        pos_weights = pos_weights
        if pos_weights is None:
            filtered_seq_len = self.filtered_seq_len(threshold)
            pos_weights = np.ones((filtered_seq_len, filtered_seq_len, aa_count, aa_count))

        jf, jfi = self.aminoacid_joint_frequencies(
            seq_weights,
            threshold,
            aa_count,
            use_tensorflow
        )

        mat_ijab = (jf
                    * np.log(jf / jfi)
                    * pos_weights)

        return np.sum(mat_ijab, axis=(2, 3))

    def vorberg_entropy(self, threshold=0.2):
        """vorberg_entropy [summary]

        Parameters
        ----------
        threshold : float, optional
            Gap-frequency low-pass value, by default .2

        Returns
        -------
        np.ndarray
            Entropy values for each position in the filtered alignment
        """
        working_alignment = self.filtered_alignment(threshold)
        f = np.apply_along_axis(
            Alignment.aa_freq_at_pos,
            0,
            working_alignment
        ).transpose(1, 0)
        s = - np.sum(f * np.log(f), axis=1)
        return s

    @staticmethod
    def basic_ICA(x, r, Niter):
        """ Basic ICA algorithm, based on work by Bell & Sejnowski (infomax). The input data should preferentially be sphered, i.e., x.T.dot(x) = 1
        Source: https://github.com/ranganathanlab/pySCA/

        Parameters
        ----------
        -  `x` = LxM input matrix where L = # features and M = # samples
        -  `r` = learning rate / relaxation parameter (e.g. r=.0001)
        -  `Niter` =  number of iterations (e.g. 1000)

        **Returns:**
        -  `w` = unmixing matrix
        -  `change` = record of incremental changes during the iterations.

        **Note:** r and Niter should be adjusted to achieve convergence, which should be assessed by visualizing 'change' with plot(range(iter) ,change)

        **Example:**
        >>> [w, change] = basicICA(x, r, Niter)

        """
        [L, M] = x.shape
        w = np.eye(L)
        change = list()
        for _ in range(Niter):
            w_old = np.copy(w)
            u = w.dot(x)
            w += r*(M*np.eye(L)+(1-2*(1./(1+np.exp(-u)))).dot(u.T)).dot(w)
            delta = (w-w_old).ravel()
            change.append(delta.dot(delta.T))
        return [w, change]

    @staticmethod
    def rot_ICA(V, kmax=6, learnrate=.0001, iterations=10000):
        """ ICA rotation (using basicICA) with default parameters and normalization of
        outputs.
        Source: https://github.com/ranganathanlab/pySCA/

        :Example:
        >>> Vica, W = rotICA(V, kmax=6, learnrate=.0001, iterations=10000)
        """
        V1 = V[:,:kmax].T
        [W, changes_s] = Alignment.basic_ICA(V1, learnrate, iterations)
        Vica = (W.dot(V1)).T
        for n in range(kmax):
            imax = abs(Vica[:, n]).argmax()
            Vica[:, n] = np.sign(Vica[imax, n])*Vica[:, n]/np.linalg.norm(Vica[:, n])
        return Vica, W

    @staticmethod
    def ICA(matrix, kmax=6, learn_rate=.0001, iterations=10000):
        v, _, _ = np.linalg.svd(matrix)
        return Alignment.rot_ICA(v, kmax, learnrate=learn_rate, iterations=iterations)
