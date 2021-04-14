import numpy as np

freq0 = np.array(
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

def parse_msa(msa_file_path, raw=False, output='numerical', unknown_aa_value='gaps'):
    with open(msa_file_path, 'r') as alg_file:
        alg_lines = alg_file.readlines()
        currt_seq = {
            'seqname': '',
            'seq': ''
        }
        for line in alg_lines:
            # This should be making use of pattern matching (available  )
            if line.startswith('>'):
                currt_seq['seqname'] = line.replace('>','')
            elif line.startswith(';'):
                continue
            else:
                currt_seq['seq'] += line

def split_chars(row):
    return [char for char in row]

def gap_frequency(npdarray, threshold=0.4):
    gap_proportion = []

    # per position
    for residue in range(npdarray.shape[1]):
        currt_freq = np.count_nonzero(npdarray[:,residue] == '-') / npdarray.shape[0]
        gap_proportion.append(currt_freq)

    gap_proportion = np.array(gap_proportion)

    working_pos = gap_proportion < threshold
    working_pos = working_pos.reshape(working_pos.shape[0], -1)

    mask = np.ones((npdarray.shape[0],1), dtype="bool").dot(working_pos.T)
    a = npdarray[mask]
    print(np.mean(a == '-'))

    return gap_proportion, np.mean(a == '-')

code_tonum = {
    '-':0,
    'X':0,
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

def reverse_code(code):
    return {v:k for k, v in code if k != 'X'}

def translate_code(arr, code):
    res = np.copy(arr)
    for k in code.keys():
        res[res == k] = code[k]
    return res.astype('int32')

def compute_freq0(freq0_src=freq0, gap_frequency):
    mult_fact = 1 - gap_frequency
    freq0_src = freq0_src * mult_fact
    freq0g = np.insert(freq0_src, 0, gap_frequency)
    return freq0g

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

def sca_entropy(alignment_array, aa_count=21, weights=None):
    _, gap_bg_freq = gap_frequency(alignment_array)
    num_rep = translate_code(alignment_array, code_tonum)

    freq0=np.array(
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

    mult_fact = 1 - gap_bg_freq
    freq0 = freq0 * mult_fact
    freq0g = np.insert(freq0, 0, gap_bg_freq)

    # Prevent MSA without gaps & regularize background frequencies the same way data is regularized
    lbda = 0.03
    freq0g_n = (1 - lbda) * freq0g + lbda / aa_count
    assert freq0g.sum() == 1

    rel_entropy = []
    for position in range(num_rep.shape[1]):
        freqs = compute_actual_frequency(num_rep[:,position], weights=weights)
        currt_rel_entropy = np.sum(freqs * np.log(freqs / freq0g))
        rel_entropy.append(currt_rel_entropy)
    return rel_entropy

def similarity(alignment_array):
    n_seqs = alignment_array.shape[0]
    res = np.zeros(shape=(n_seqs, n_seqs))
    for seq_1_i in range(n_seqs):
        seq_1 = alignment_array[seq_1_i]
        for seq_2_i in range(seq_1_i, n_seqs):
            seq_2 = alignment_array[seq_2_i]
            dist = np.mean(np.array([seq_1[i]==seq_2[i] for i in range(len(seq_1))]))
            res[seq_1_i,seq_2_i]=dist
            res[seq_2_i,seq_1_i]=dist
    return res

def filter_overly_gapped_positions(alignment, gap_frequency, threshold):
    return alignment[:,gap_frequency < threshold]