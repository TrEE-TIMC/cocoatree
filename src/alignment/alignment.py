import numpy as np

class Alignment:
    # Used to translate AAs to Ints to make use of Numpy vectorization
    _code_tonum = {
        '-':0,
        'X':0, # Consider unknown AAs as gaps
        'A':1,
        'C':2,
        'D':3,
        'E':4,
        'F':5,
        'G':6,
        'H':7,
        'I':8,
        'K':9,
        'L':10,
        'M':11,
        'N':12,
        'P':13,
        'Q':14,
        'R':15,
        'S':16,
        'T':17,
        'V':18,
        'W':19,
        'Y':20
    }

    _kws = {
        'text': '_text_rep',
        'numerical': '_num_rep',
        'num': '_num_rep'
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

    def __init__(self, alignment_array) -> None:
        """
        alignment_array: list of strings corresponding to one string per sequence
        """
        self._raw = alignment_array
        # Split characters so that one element in each row represents one aminoacid
        self._text_rep = np.array([np.array([char for char in row]) for row in alignment_array])
        self._num_rep = np.zeros(self._text_rep.shape, dtype=np.int64)
        for k in Alignment._code_tonum.keys():
            self._num_rep[self._text_rep == k] = self._code_tonum[k]

    def gap_frequency(self, threshold=0.2):
        """Computes the gap frequency for a MSA

        Keyword arguments:
        threshold: a floating point number. Positions will be considered overly-gapped if their gap ratio is above this value

        Returns:
        (gap_frequency, background_gap_frequency): a tuple consisting of:
        - gap_frequency: a numpy 1D array of length n_pos with each value
                         set to the gap frequency at the position
                         corresponding to the index
        - background_gap_frequency: a floating point number containing the gap
                                    frequency accross all positions which are
                                    below the desired threshold
        """
        gap_frequency = []

        # per position
        for residue in range(self._text_rep.shape[1]):
            currt_freq = np.count_nonzero(self._text_rep[:,residue] == '-') / self._text_rep.shape[0]
            gap_frequency.append(currt_freq)

        gap_frequency = np.array(gap_frequency)

        working_pos = gap_frequency < threshold
        working_pos = working_pos.reshape(working_pos.shape[0], -1)

        mask = np.ones((self._text_rep.shape[0],1), dtype="bool").dot(working_pos.T)
        bg_freq_proportion = self._text_rep[mask]

        return gap_frequency, np.mean(bg_freq_proportion == '-')

    def get_filtered_alignment(self, threshold=0.2):
        """Returns the filtered alignment
        """
        gf, _ = self.gap_frequency(threshold)
        return self._num_rep[:,gf < threshold].copy()

    def similarity(self, use_filtered=True, threshold=0.2):
        """Computes the similarity between sequences in a MSA

        Keyword arguments:
        use_filtered: make use of the filtered MSA, according to the gap frequency threshold (default True)
        threshold: the threshold to filter overly-gapped positions, unused if use_filtered is False (default 0.2)

        Returns:
        A numpy 2D square matrix of size n_seq * n_seq (with n_seq representing the number of sequences in the alignment)
        """
        src = self.get_filtered_alignment(threshold) if use_filtered else self._num_rep
        n_seqs = src.shape[0]
        res = np.zeros(shape=(n_seqs, n_seqs))
        for seq_1_i in range(n_seqs):
            seq_1 = src[seq_1_i]
            for seq_2_i in range(seq_1_i, n_seqs):
                seq_2 = src[seq_2_i]
                dist = np.mean(np.array([seq_1[i]==seq_2[i] for i in range(len(seq_1))]))
                res[seq_1_i,seq_2_i]=dist
                res[seq_2_i,seq_1_i]=dist
        return res

    def compute_actual_frequency(algnt_col, lbda = 0.03, aa_count=21, weights=None):
        # AA proportion
        bins = np.bincount(algnt_col, weights=weights)

        # Pad with zeroes for AAs not present at currt position
        res = np.pad(bins, (0, aa_count - bins.shape[0]))
        res = np.array(res)

        # Frequency
        if weights is None:
            weights = np.ones(algnt_col.shape)

        res = res / weights.sum()

        # Regularization
        res = (1 - lbda) * res + lbda / aa_count
        return res

    def sca_entropy(self, aa_count=21, weights=None, threshold=0.2):
        """Computes entropy as per SCA guidelines

        Keyword arguments:
        aa_count: the number of aminoacids to consider
        weights: to be applied to the different sequences in the MSA
        threshold: the gap frequency low-pass threshold

        Returns:
        float: the entropy considered for each of the positions in the filtered alignment
        """
        _, gap_bg_freq = self.gap_frequency(threshold)
        working_alignment = self.get_filtered_alignment(threshold)

        mult_fact = 1 - gap_bg_freq
        freq0 = Alignment._freq0 * mult_fact
        freq0g = np.insert(freq0, 0, gap_bg_freq)

        # Prevent MSA without gaps & regularize background frequencies the same way data is regularized

        # freq0g_n = (1 - Alignment._lbda) * freq0g + Alignment._lbda / aa_count
        assert freq0g.sum() == 1

        rel_entropy = []
        for position in range(working_alignment.shape[1]):
            freqs = Alignment.compute_actual_frequency(working_alignment[:,position], weights=weights)
            currt_rel_entropy = np.sum(freqs * np.log(freqs / freq0g))
            rel_entropy.append(currt_rel_entropy)
        return rel_entropy

    def get_similarity_weights(self, similarity, threshold):
        """Get the effective number of sequences after applying weights

        Keyword arguments:
        similarity: a similarity matrix
        threshold: the "width" (delta) of an effective sequence
        """
        return (1 / np.sum(similarity >= threshold, axis=0))

    def weight_bins(self, threshold_increment=.001, low_bin_count=100, weight_vector_count=5, use_filtered=True, threshold=0.2):
        """Computes weight bins using a similarity matrix, in order to regroup "similar-enough" sequences in bins

        Keyword arguments:
        threshold_increment: the rate at which delta should be decremented when trying to build bins
        low_bin_count: the target effective sequence count for the smallest bin
        weight_vector_count: the number of bins to create
        use_filtered: use a filtered version of the alignment according to gap cutoffs
        threshold: the threshold over which positions are considered overly-gapped
        """
        similarity_matrix = self.similarity(use_filtered=use_filtered, threshold=threshold)
        bin_sizes = np.logspace(np.log10(low_bin_count),np.log10(similarity_matrix.shape[0]),weight_vector_count).tolist()
        delta_list = []
        weight_list = []
        currt_delta = np.max(similarity_matrix)
        fnd_values = len(delta_list)
        it = 0
        while fnd_values != weight_vector_count:
            sim_wt = self.get_similarity_weights(similarity_matrix, currt_delta)
            eff_seqs = sum(sim_wt)
            # If the tree is very well structured, there is no point in trying to find values very close to the generated bin sizes
            if eff_seqs < bin_sizes[-1 - fnd_values]:
                delta_list.append(currt_delta)
                weight_list.append(sim_wt)
                fnd_values += 1
            currt_delta -= threshold_increment
            it += 1
        return delta_list, weight_list

    def get_sequence_count(self):
        return self._num_rep.shape[0]

    def compute_most_frequent_aa(self, filtered=True, threshold=.2):
        most_freq_aa = []
        most_freq_aa_freq = []
        src = self.get_filtered_alignment()
        for i in range(src.shape[1]):
            currt_col = src[:,i]
            currt_bins = np.bincount(currt_col)
            currt_max_i = np.argmax(currt_bins)
            most_freq_aa.append(currt_max_i)
            most_freq_aa_freq.append(currt_bins[currt_max_i] / np.sum(currt_bins))
        return np.array(most_freq_aa), np.array(most_freq_aa_freq)